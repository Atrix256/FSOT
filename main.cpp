#define _CRT_SECURE_NO_WARNINGS // for stb

// Settings
#define DETERMINISTIC() true
#define MULTITHREADED() true
static const int c_pointImageSize = 256;

static const bool c_writePointSetTxt = false;

static const bool c_drawNonGaussImage = true;

static const bool c_drawGaussImage = false;
static const int c_pointImageGaussSize = 512;
static const float c_pointImageGaussBlobSigma = 1.5f;

#include <random>
#include <vector>
#include <direct.h>
#include <stdio.h>
#include <chrono>

#include "squarecdf.h"
#include "NumericalCDF.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

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

class FSOTClassBase
{
public:
	virtual float ICDF(float y, const float2& direction) const = 0;
	virtual void ProjectPoints(const std::vector<float2>& points, float subclassZ, const float2& direction, std::vector<float>& projections) const = 0;
};

template <typename TICDF, typename TFilter>
class FSOTClass : public FSOTClassBase, public TICDF, public TFilter
{
	float ICDF(float y, const float2& direction) const override final
	{
		return TICDF::ICDF(y, direction);
	}

	void ProjectPoints(const std::vector<float2>& points, float subclassZ, const float2& direction, std::vector<float>& projections) const override final
	{
		return TFilter::ProjectPoints(points, subclassZ, direction, projections);
	}
};

class ICDF_UniformSquare
{
public:
	float ICDF(float y, const float2& direction) const
	{
		// Convert y: square is in [-0.5, 0.5], but y is in [0, 1].
		y = y - 0.5f;

		// Evaluate ICDF
		float x = Square::InverseCDF(y, direction);

		// The CDF is in [-0.5, 0.5], but we want the points to be in [-1,1]
		return x * 2.0f;
	}
};

class Filter_All
{
public:
	void ProjectPoints(const std::vector<float2>& points, float subclassZ, const float2& direction, std::vector<float>& projections) const
	{
		// project all the points
		projections.resize(points.size());
		for (size_t i = 0; i < points.size(); ++i)
			projections[i] = Dot(direction, points[i]);
	}
};

template <int NumChannels>
void PlotGaussian(std::vector<unsigned char>& image, int width, int height, int x, int y, float sigma, unsigned char color[NumChannels])
{
	int kernelRadius = int(std::sqrt(-2.0f * sigma * sigma * std::log(0.005f)));

	int sx = Clamp(x - kernelRadius, 0, width - 1);
	int ex = Clamp(x + kernelRadius, 0, height - 1);
	int sy = Clamp(y - kernelRadius, 0, width - 1);
	int ey = Clamp(y + kernelRadius, 0, height - 1);

	for (int iy = sy; iy <= ey; ++iy)
	{
		unsigned char* pixel = &image[(iy * width + sx) * NumChannels];

		int ky = std::abs(iy - y);
		float kernelY = std::exp(-float(ky * ky) / (2.0f * sigma * sigma));

		for (int ix = sx; ix <= ex; ++ix)
		{
			int kx = std::abs(ix - x);
			float kernelX = std::exp(-float(kx * kx) / (2.0f * sigma * sigma));

			float kernel = kernelX * kernelY;

			for (int i = 0; i < NumChannels; ++i)
			{
				unsigned char oldColor = *pixel;
				unsigned char newColor = (unsigned char)Lerp(float(oldColor), float(color[i]), kernel);
				*pixel = newColor;
				pixel++;
			}
		}
	}
}

void SavePointSet(const std::vector<float2>& points, const char* baseFileName, int index, int total)
{
	// Write out points in text
	if (c_writePointSetTxt)
	{
		char fileName[1024];
		sprintf_s(fileName, "%s_%i_%i.txt", baseFileName, index, total);
		FILE* file = nullptr;
		fopen_s(&file, fileName, "wb");

		fprintf(file, "float points[][2] =\n{\n");

		for (size_t index = 0; index < points.size(); ++index)
			fprintf(file, "    { %ff, %ff },\n", points[index].x, points[index].y);

		fprintf(file, "};\n");

		fclose(file);
	}

	// Draw an image of the points
	{
		std::vector<unsigned char> pixels(c_pointImageSize * c_pointImageSize, 255);
		std::vector<unsigned char> pixelsGauss(c_pointImageGaussSize * c_pointImageGaussSize, 255);

		for (size_t index = 0; index < points.size(); ++index)
		{
			if (c_drawNonGaussImage)
			{
				int x = (int)Clamp((points[index].x * 0.5f + 0.5f) * float(c_pointImageSize - 1), 0.0f, float(c_pointImageSize - 1));
				int y = (int)Clamp((points[index].y * 0.5f + 0.5f) * float(c_pointImageSize - 1), 0.0f, float(c_pointImageSize - 1));
				pixels[y * c_pointImageSize + x] = 0;
			}

			if (c_drawGaussImage)
			{
				unsigned char color[] = { 0 };
				int x = (int)Clamp((points[index].x * 0.5f + 0.5f) * float(c_pointImageGaussSize - 1), 0.0f, float(c_pointImageGaussSize - 1));
				int y = (int)Clamp((points[index].y * 0.5f + 0.5f) * float(c_pointImageGaussSize - 1), 0.0f, float(c_pointImageGaussSize - 1));
				PlotGaussian<1>(pixelsGauss, c_pointImageGaussSize, c_pointImageGaussSize, x, y, c_pointImageGaussBlobSigma, color);
			}
		}

		if (c_drawNonGaussImage)
		{
			char fileName[1024];
			sprintf_s(fileName, "%s_%i_%i.png", baseFileName, index, total);
			stbi_write_png(fileName, c_pointImageSize, c_pointImageSize, 1, pixels.data(), 0);
		}

		if (c_drawGaussImage)
		{
			char fileName[1024];
			sprintf_s(fileName, "%s_%i_%i.gauss.png", baseFileName, index, total);
			stbi_write_png(fileName, c_pointImageGaussSize, c_pointImageGaussSize, 1, pixelsGauss.data(), 0);
		}
	}
}

void GeneratePoints(int numPoints, int batchSize, int numIterations, const char* baseFileName, int numProgressImages, bool useGradientScalingFactor, bool stratifyLine, const std::vector<const FSOTClassBase*>& classes)
{
	// get the timestamp of when this started
	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

	printf("==================================\n%s\n==================================\n", baseFileName);

	FILE* file = nullptr;
	char outputFileNameCSV[1024];
	sprintf(outputFileNameCSV, "%s.csv", baseFileName);
	fopen_s(&file, outputFileNameCSV, "wb");
	fprintf(file, "\"Iteration\",\"Avg. Movement\"\n");

	// Generate starting points
	std::vector<float2> points(numPoints);
	{
		std::mt19937 rng = GetRNG(0);
		std::uniform_real_distribution<float> dist(-1.0f, 1.0f);

		for (float2& p : points)
		{
			p.x = dist(rng);
			p.y = dist(rng);
		}
	}
	std::vector<float2> startingPoints = points;

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
		std::vector<float2> batchDirections;
	};
	std::vector<BatchData> allBatchData(batchSize, BatchData(numPoints));

	// FSOT
	int lastPercent = -1;
	for (int iterationIndex = 0; iterationIndex < numIterations; ++iterationIndex)
	{
		// Write out progress
		if (numProgressImages > 0)
		{
			int progressInterval = numIterations / numProgressImages;
			if (iterationIndex % progressInterval == 0)
				SavePointSet(points, baseFileName, iterationIndex / progressInterval, numProgressImages);
		}

		// Do the batches in parallel
		#if MULTITHREADED()
		#pragma omp parallel for
		#endif
		for (int batchIndex = 0; batchIndex < batchSize; ++batchIndex)
		{
			BatchData& batchData = allBatchData[batchIndex];

			std::mt19937 rng = GetRNG(iterationIndex * batchSize + batchIndex);

			// generate a random projection direction
			float2 direction;
			{
				std::normal_distribution<float> distNormal(0.0f, 1.0f);

				direction.x = distNormal(rng);
				direction.y = distNormal(rng);
				direction = Normalize(direction);
			}

			// Select a class randomly if there is more than 1
			int selectedClass = 0;
			if (classes.size() > 1)
			{
				std::uniform_int_distribution<int> dist(0, (int)classes.size() - 1);
				selectedClass = dist(rng);
			}
			const FSOTClassBase& FSOTClass = *classes[selectedClass];

			// Select a floating point value [0,1] for sub class selection
			float subclassZ = 0.0f;
			{
				std::uniform_real<float> dist(0.0f, 1.0f);
				subclassZ = dist(rng);
			}

			// Project the points
			FSOTClass.ProjectPoints(points, subclassZ, direction, batchData.projections);

			// sort the projections
			batchData.sorted.resize(batchData.projections.size());
			for (int i = 0; i < (int)batchData.sorted.size(); ++i)
				batchData.sorted[i] = i;
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
				if (stratifyLine)
					jitter = distJitter(rng);
				float projectionTarget = FSOTClass.ICDF((float(i) + jitter) / float(batchData.sorted.size()), direction);

				float projectionDiff = projectionTarget - projection;

				float gradientFactor = 1.0f;
				if (useGradientScalingFactor)
				{
					int lastIndex = std::max(int(i) - 1, 0);
					int nextIndex = std::min(int(i) + 1, int(batchData.sorted.size()) - 1);
					gradientFactor = (FSOTClass.ICDF((float(nextIndex)) / float(batchData.sorted.size()), direction) - FSOTClass.ICDF((float(lastIndex)) / float(batchData.sorted.size()), direction)) / float(nextIndex - lastIndex);
					gradientFactor *= float(batchData.projections.size());
				}

				float delta = pointCountScalingFactor * projectionDiff / gradientFactor;

				batchData.batchDirections[batchData.sorted[i]].x = direction.x * delta;
				batchData.batchDirections[batchData.sorted[i]].y = direction.y * delta;
			}
		}

		// average all batch directions into batchDirections[0]
		{
			for (int batchIndex = 1; batchIndex < batchSize; ++batchIndex)
			{
				float alpha = 1.0f / float(batchIndex + 1);
				for (size_t i = 0; i < numPoints; ++i)
				{
					allBatchData[0].batchDirections[i].x = Lerp(allBatchData[0].batchDirections[i].x, allBatchData[batchIndex].batchDirections[i].x, alpha);
					allBatchData[0].batchDirections[i].y = Lerp(allBatchData[0].batchDirections[i].y, allBatchData[batchIndex].batchDirections[i].y, alpha);
				}
			}
		}

		// update points
		float totalDistance = 0.0f;
		for (size_t i = 0; i < numPoints; ++i)
		{
			const float2& adjust = allBatchData[0].batchDirections[i];

			points[i].x += adjust.x;
			points[i].y += adjust.y;

			totalDistance += std::sqrt(adjust.x * adjust.x + adjust.y * adjust.y);
		}

		int percent = int(100.0f * float(iterationIndex) / float(numIterations - 1));
		if (percent != lastPercent)
		{
			lastPercent = percent;
			printf("\r[%i%%] %f", percent, totalDistance / float(numPoints));
			fprintf(file, "\"%i\",\"%f\"\n", iterationIndex, totalDistance / float(numPoints));
		}
	}
	printf("\n");

	fclose(file);

	// Write out the final results
	SavePointSet(points, baseFileName, numProgressImages, numProgressImages);

	// report how long this took
	float elpasedSeconds = std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::high_resolution_clock::now() - start).count();
	printf("%0.2f seconds\n\n", elpasedSeconds);
}

int main(int argc, char** argv)
{
	_mkdir("out");

	// 1) Show the benefit of the gradient scaling fix
	// 2) Compare stratifying the line vs not
	{
		FSOTClass<ICDF_UniformSquare, Filter_All> a;

		// 1) Show the benefit of the gradient scaling fix
		GeneratePoints(1000, 64, 1000, "out/GradFixNoStratifyNo", 1, false, false, { &a });
		GeneratePoints(1000, 64, 1000, "out/GradFixYesStratifyNo", 1, true, false, { &a });

		// 2) Compare stratifying the line vs not
		GeneratePoints(1000, 64, 1000, "out/GradFixYesStratifyYes", 1, true, true, { &a });
	}

	// Multiclass
	{
		// TODO: separate ICDF from membership function, make a templated object take one of each.

		//FSOTClass_UniformSquare a;
		//FSOTClass_UniformSquare b;
	}

	return 0;
}

/*
A class has...
1) A membership lambda. Take a Z value and point count, it returns a vector of points that are in the membership.
2) a target desnity. like an ICDF to project?


FSOT PAPER NOTES:
* FSOT paper has batch size of 64 by default. 4096 iterations. 4096 points.
* they make random points in space and project (not jitter, random!), instead of doing equal area points on the line. getInverseRadonNCube()

BLOG:
* first, compare GradFixNo to GradFixYes, showing how it improves that "overconvergence" thing. should have a better DFT than not letting it go to convergence (compare the 3 point sets, and DFTs)
* Note that you are making the target be the center of each bucket, put through the ICDF, but they are stratifying. Compare the 2, see if there's a difference. note that not timothy lottes mentioned that too, after the last post.

TODO:
* profiling seems to say that GetRNG is a hot spot. maybe generate all RNG in advance?
* could try showing DFTs that are the average of several realizations. may help show if stratification is good or not?
* clean up generate points function
* you need to be able to explain the gradient scaling for your blog post. why does it improve things?
* point sets, and textures
* progressive point set.
* toroidally progressive point set (secret for now? JCGT?)
? what does a continuous membership even mean? maybe show 1 class with a box membership function vs a smooth membership function.
 * with fewer points (higher Z value) in continuous, those points get spread out more. weird.
? how do they support toroidal distance for their points?
* profile again

Toroidal Progressiveness
* compare vs halton, and R2
* projective blue noise and not projective

*/
