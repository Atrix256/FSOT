#pragma once

#include <random>
#include <chrono>

#define DETERMINISTIC() true
#define MULTITHREADED() true

template <typename TVecType>
class FSOTClassBase
{
public:
	FSOTClassBase(float weight = 1.0f)
		: m_weight(weight)
	{
	}
	
	virtual ~FSOTClassBase()
	{
	}

	float m_weight = 1.0f;

	virtual float ICDF(float y, const TVecType& direction) const = 0;
	virtual void Filter(const std::vector<TVecType>& points, float subClassZ, std::vector<int>& memberPoints) const = 0;
	virtual bool GetRandomDirection(std::mt19937& rng, TVecType& direction) const = 0;
};

template <typename TVecType, typename TCDF, typename TFilter>
class FSOTClass : public FSOTClassBase<TVecType>
{
public:

	FSOTClass(float weight = 1.0f)
		: FSOTClassBase<TVecType>(weight)
	{
	}

	float ICDF(float y, const TVecType& direction) const override final
	{
		return m_cdf.ICDF(y, direction);
	}

	void Filter(const std::vector<TVecType>& points, float subClassZ, std::vector<int>& memberPoints) const override final
	{
		return m_filter.Filter(points, subClassZ, memberPoints);
	}

	bool GetRandomDirection(std::mt19937& rng, TVecType& direction) const override final
	{
		return m_cdf.OverrideDirection(direction);
	}

	TCDF m_cdf;
	TFilter m_filter;
};


template <typename TVecType>
class FSOT
{
public:
	FSOT()
	{

	}

	~FSOT()
	{
		for (FSOTClassBase<TVecType>* c : m_classes)
			delete c;
		m_classes.clear();
	}

	std::mt19937 GetRNG(int seed)
	{
	#if DETERMINISTIC()
		std::mt19937 ret(seed);
	#else
		std::random_device rd;
		std::mt19937 ret(rd());
	#endif
		return ret;
	}

	template <typename TCDF, typename TFilter>
	FSOTClass<TVecType, TCDF, TFilter>& AddClass(float weight = 1.0f)
	{
		FSOTClass<TVecType, TCDF, TFilter>* newClass = new FSOTClass<TVecType, TCDF, TFilter>(weight);
		m_classes.push_back(newClass);
		return *newClass;
	}

	std::vector<FSOTClassBase<TVecType>*> m_classes;

	// settings
	int m_batchSize = 64;
	int m_numIterations = 10000;
	bool m_useGradientScalingFactor = true;
	bool m_stratifyLine = false;

	// Debug output settings
	bool m_showProgress = true;
	bool m_saveAvgMovementCSV = true;

	template <typename TDebugOutput, typename TGenerateRandomPointFn, typename TGenerateRandomDirectionFn, typename TProjectPointsFn>
	std::vector<TVecType> GeneratePoints(int numPoints, const char* baseFileName, int numProgressImages, const TGenerateRandomPointFn& GenerateRandomPointFn, const TGenerateRandomDirectionFn& GenerateRandomDirectionFn, const TProjectPointsFn& ProjectPointsFn)
	{
		// get the timestamp of when this started
		std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

		if (m_showProgress)
			printf("==================================\n%s\n==================================\n", baseFileName);

		TDebugOutput debugOuput;

		FILE* avgMovementCSV = nullptr;
		if (m_saveAvgMovementCSV)
		{
			char outputFileNameCSV[1024];
			sprintf_s(outputFileNameCSV, "%s.movement.csv", baseFileName);
			fopen_s(&avgMovementCSV, outputFileNameCSV, "wb");
			fprintf_s(avgMovementCSV, "\"Iteration\",\"Avg. Movement\"\n");
		}

		// Generate starting points
		std::vector<TVecType> points(numPoints);
		{
			std::mt19937 rng = GetRNG(0);
			for (TVecType& p : points)
				GenerateRandomPointFn(rng, p);
		}
		std::vector<TVecType> startingPoints = points;

		// Per batch data
		// Each batch entry has it's own data so the batches can be parallelized
		struct BatchData
		{
			BatchData(int numPoints)
			{
				sorted.resize(numPoints);
				for (int i = 0; i < numPoints; ++i)
					sorted[i] = i;
				projections.resize(numPoints);
				batchDirections.resize(numPoints);
			}

			std::vector<int> sorted;
			std::vector<float> projections;
			std::vector<TVecType> batchDirections;
		};
		std::vector<BatchData> allBatchData(m_batchSize, BatchData(numPoints));

		// FSOT
		int lastPercent = -1;
		for (int iterationIndex = 0; iterationIndex < m_numIterations; ++iterationIndex)
		{
			// Write out progress
			debugOuput.Output(points, baseFileName, iterationIndex, m_numIterations);

			// Do the batches in parallel
			#if MULTITHREADED()
			#pragma omp parallel for
			#endif
			for (int batchIndex = 0; batchIndex < m_batchSize; ++batchIndex)
			{
				BatchData& batchData = allBatchData[batchIndex];

				std::mt19937 rng = GetRNG(iterationIndex * m_batchSize + batchIndex);

				// Select a class randomly by weight
				int selectedClass = 0;
				if (m_classes.size() > 1)
				{
					float totalWeight = 0.0f;
					for (const FSOTClassBase<TVecType>* fsotClass : m_classes)
						totalWeight += fsotClass->m_weight;

					std::uniform_real<float> dist(0.0f, totalWeight);

					float weight = dist(rng);

					weight -= m_classes[0]->m_weight;
					while (weight > 0.0f && selectedClass < m_classes.size() - 1)
					{
						selectedClass++;
						weight -= m_classes[selectedClass]->m_weight;
					}
				}
				const FSOTClassBase<TVecType>& FSOTClass = *m_classes[selectedClass];

				// Select a floating point value [0,1] for sub class selection
				float subClassZ = 0.0f;
				{
					std::uniform_real<float> dist(0.0f, 1.0f);
					subClassZ = dist(rng);
				}

				// Filter the points so that batchData.sorted contains the point indices that are members of the subclass
				FSOTClass.Filter(points, subClassZ, batchData.sorted);

				// if no points selected, don't do anything
				if (batchData.sorted.size() == 0)
					continue;

				// generate a random projection direction
				TVecType direction;
				if (!FSOTClass.GetRandomDirection(rng, direction))
					GenerateRandomDirectionFn(rng, direction);

				// Project the points
				for (int i : batchData.sorted)
					batchData.projections[i] = ProjectPointsFn(direction, points[i]);

				// sort the projections
				std::sort(batchData.sorted.begin(), batchData.sorted.end(),
					[&](uint32_t a, uint32_t b)
					{
						return batchData.projections[a] < batchData.projections[b];
					}
				);

				// update batchDirections
				std::uniform_real_distribution<float> distJitter(0.0f, 1.0f);
				float pointCountScalingFactor = float(batchData.projections.size()) / float(numPoints);
				for (size_t i = 0; i < batchData.sorted.size(); ++i)
				{
					float projection = batchData.projections[batchData.sorted[i]];

					float jitter = 0.5f;
					if (m_stratifyLine)
						jitter = distJitter(rng);
					float projectionTarget = FSOTClass.ICDF((float(i) + jitter) / float(batchData.sorted.size()), direction);

					float projectionDiff = projectionTarget - projection;

					float gradientFactor = 1.0f;
					if (m_useGradientScalingFactor && batchData.sorted.size() > 1)
					{
						int lastIndex = std::max(int(i) - 1, 0);
						int nextIndex = std::min(int(i) + 1, int(batchData.sorted.size()) - 1);
						gradientFactor = (FSOTClass.ICDF((float(nextIndex)) / float(batchData.sorted.size()), direction) - FSOTClass.ICDF((float(lastIndex)) / float(batchData.sorted.size()), direction)) / float(nextIndex - lastIndex);
						gradientFactor *= float(batchData.projections.size());
					}

					float delta = pointCountScalingFactor * projectionDiff / gradientFactor;

					batchData.batchDirections[batchData.sorted[i]] = direction * delta;
				}
			}

			// average all batch directions into batchDirections[0]
			{
				for (int batchIndex = 1; batchIndex < m_batchSize; ++batchIndex)
				{
					float alpha = 1.0f / float(batchIndex + 1);
					for (size_t i = 0; i < numPoints; ++i)
						allBatchData[0].batchDirections[i] = Lerp(allBatchData[0].batchDirections[i], allBatchData[batchIndex].batchDirections[i], alpha);
				}
			}

			// update points
			float totalDistance = 0.0f;
			for (size_t i = 0; i < numPoints; ++i)
			{
				const TVecType& adjust = allBatchData[0].batchDirections[i];
				points[i] = points[i] + adjust;
				totalDistance += Length(adjust);
			}

			int percent = int(100.0f * float(iterationIndex) / float(m_numIterations - 1));
			if (percent != lastPercent)
			{
				lastPercent = percent;
				if (m_showProgress)
					printf("\r[%i%%] %f", percent, totalDistance / float(numPoints));

				if (m_saveAvgMovementCSV)
					fprintf(avgMovementCSV, "\"%i\",\"%f\"\n", iterationIndex, totalDistance / float(numPoints));
			}
		}
		if (m_showProgress)
			printf("\n");

		if (avgMovementCSV)
			fclose(avgMovementCSV);

		// Write out the final results
		debugOuput.Output(points, baseFileName, m_numIterations, m_numIterations);

		// report how long this took
		if (m_showProgress)
		{
			float elpasedSeconds = std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::high_resolution_clock::now() - start).count();
			printf("%0.2f seconds\n\n", elpasedSeconds);
		}

		return points;
	}
};
