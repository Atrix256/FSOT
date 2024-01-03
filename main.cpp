#define _CRT_SECURE_NO_WARNINGS // for stb

// Settings
#define DETERMINISTIC() true
#define MULTITHREADED() true

static const bool c_writePointSetTxt = false;

static const bool c_drawStartingState = false;

static const bool c_drawPointImage = true;
static const int c_pointImageSize = 256;

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
	virtual void Filter(const std::vector<float2>& points, float subClassZ, std::vector<int>& memberPoints) const = 0;
};

template <typename TICDF, typename TFilter>
class FSOTClass : public FSOTClassBase, public TICDF, public TFilter
{
	float ICDF(float y, const float2& direction) const override final
	{
		return TICDF::ICDF(y, direction);
	}

	void Filter(const std::vector<float2>& points, float subClassZ, std::vector<int>& memberPoints) const override final
	{
		return TFilter::Filter(points, subClassZ, memberPoints);
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
	void Filter(const std::vector<float2>& points, float subClassZ, std::vector<int>& memberPoints) const
	{
		// All points
		memberPoints.resize(points.size());
		for (int i = 0; i < (int)points.size(); ++i)
			memberPoints[i] = i;
	}
};

template <float THE_MIN, float THE_MAX>
class Filter_Range
{
public:
	void Filter(const std::vector<float2>& points, float subClassZ, std::vector<int>& memberPoints) const
	{
		int indexMin = std::max(int(THE_MIN * float(points.size())), 0);
		int indexMax = std::min(int(THE_MAX * float(points.size())), int(points.size()));

		memberPoints.resize(indexMax - indexMin);
		for (int i = indexMin; i < indexMax; ++i)
			memberPoints[i - indexMin] = i;
	}
};

class Filter_Progressive
{
public:
	void Filter(const std::vector<float2>& points, float subClassZ, std::vector<int>& memberPoints) const
	{
		int indexMax = Clamp(int(subClassZ * float(points.size())), 0, int(points.size()));
		memberPoints.resize(indexMax);
		for (int i = 0; i < indexMax; ++i)
			memberPoints[i] = i;
	}
};

class PointColoringObjectBase
{
public:
	virtual int NumSets(const std::vector<float2>& points) const = 0;
	virtual std::vector<float2> GetSet(const std::vector<float2>& points, int setIndex) const = 0;
	virtual void GetSetColor(int setIndex, unsigned char& R, unsigned char& G, unsigned char& B) const = 0;
	virtual void GetBackgroundColor(unsigned char& R, unsigned char& G, unsigned char& B) const = 0;
	virtual bool ShowSetPermutations() const = 0;
};

class PointColoringObject_OneSetBlack : public PointColoringObjectBase
{
public:
	int NumSets(const std::vector<float2>& points) const override final
	{
		return 1;
	}

	std::vector<float2> GetSet(const std::vector<float2>& points, int setIndex) const override final
	{
		std::vector<float2> ret(points.size());
		ret = points;
		return ret;
	}

	void GetSetColor(int setIndex, unsigned char& R, unsigned char& G, unsigned char& B) const override final
	{
		R = 0;
		G = 0;
		B = 0;
	}

	void GetBackgroundColor(unsigned char& R, unsigned char& G, unsigned char& B) const override final
	{
		R = 255;
		G = 255;
		B = 255;
	}

	bool ShowSetPermutations() const override final
	{
		return true;
	}
};

template <int NUM_SETS = 0>
class PointColoringObject_Progressive : public PointColoringObjectBase
{
public:
	int NumSets(const std::vector<float2>& points) const override final
	{
		return NUM_SETS ? NUM_SETS : (int)points.size();
	}

	std::vector<float2> GetSet(const std::vector<float2>& points, int setIndex) const override final
	{
		if (NUM_SETS == 0)
		{
			std::vector<float2> ret(setIndex + 1);
			for (int i = 0; i < setIndex; ++i)
				ret[i] = points[i];
			return ret;
		}
		else
		{
			int setEnd = int(float(setIndex + 1) * float(points.size()) / float(NUM_SETS));
			std::vector<float2> ret(setEnd);
			for (int i = 0; i < setEnd; ++i)
				ret[i] = points[i];
			return ret;
		}
	}

	void GetSetColor(int setIndex, unsigned char& R, unsigned char& G, unsigned char& B) const override final
	{
		R = 0;
		G = 0;
		B = 0;
	}

	void GetBackgroundColor(unsigned char& R, unsigned char& G, unsigned char& B) const override final
	{
		R = 255;
		G = 255;
		B = 255;
	}

	bool ShowSetPermutations() const override final
	{
		return false;
	}
};

class PointColoringObject_HalfRedHalfBlue: public PointColoringObjectBase
{
public:
	int NumSets(const std::vector<float2>& points) const override final
	{
		return 2;
	}

	std::vector<float2> GetSet(const std::vector<float2>& points, int setIndex) const override final
	{
		int start = (setIndex == 0) ? 0 : int(points.size()) / 2;
		int end = (setIndex == 0) ? int(points.size()) / 2 : int(points.size());

		std::vector<float2> ret(end - start);
		for (int i = start; i < end; ++i)
			ret[i - start] = points[i];

		return ret;
	}

	void GetSetColor(int setIndex, unsigned char& R, unsigned char& G, unsigned char& B) const override final
	{
		R = (setIndex == 0) ? 255 : 0;
		G = 0;
		B = (setIndex == 1) ? 255 : 0;
	}

	void GetBackgroundColor(unsigned char& R, unsigned char& G, unsigned char& B) const override final
	{
		R = 255;
		G = 255;
		B = 255;
	}

	bool ShowSetPermutations() const override final
	{
		return true;
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

void SavePointSet(const std::vector<float2>& points, const char* baseFileName, int progressIndex, int progressTotal, const PointColoringObjectBase& pointColoringObject)
{
	const int numSets = pointColoringObject.NumSets(points);

	// Write out points in text
	if (c_writePointSetTxt)
	{
		char fileName[1024];
		sprintf_s(fileName, "%s_%i_%i.txt", baseFileName, progressIndex, progressTotal);
		FILE* file = nullptr;
		fopen_s(&file, fileName, "wb");

		if (numSets > 1)
		{
			fprintf(file, "float points[][3] =\n{\n");

			for (int setIndex = 0; setIndex < numSets; ++setIndex)
			{
				fprintf(file, "\n    // Class %i\n", setIndex);
				for (const float2& p : pointColoringObject.GetSet(points, setIndex))
					fprintf(file, "    { %ff, %ff, %ff },\n", p.x, p.y, float(setIndex));
			}
		}
		else
		{
			fprintf(file, "float points[][2] =\n{\n");
			for (size_t index = 0; index < points.size(); ++index)
				fprintf(file, "    { %ff, %ff },\n", points[index].x, points[index].y);
		}

		fprintf(file, "};\n");

		fclose(file);
	}

	// Draw an image of the points
	int imageCount = 0;
	if (c_drawPointImage || c_drawGaussImage)
	{
		unsigned char bg[3];
		pointColoringObject.GetBackgroundColor(bg[0], bg[1], bg[2]);

		std::vector<unsigned char> pixels(c_pointImageSize * c_pointImageSize * 3);
		std::vector<unsigned char> pixelsGauss(c_pointImageGaussSize * c_pointImageGaussSize * 3, 255);

		bool showSetPermutations = pointColoringObject.ShowSetPermutations();
		int loopIndexStart = 0;
		int loopIndexEnd = numSets;
		if (showSetPermutations)
		{
			loopIndexStart = 1;
			loopIndexEnd = 1 << numSets;
		}

		for (int loopIndex = loopIndexStart; loopIndex < loopIndexEnd; ++loopIndex)
		{
			for (size_t index = 0; index < c_pointImageSize * c_pointImageSize; ++index)
				memcpy(&pixels[index * 3], bg, 3);

			for (size_t index = 0; index < c_pointImageGaussSize * c_pointImageGaussSize; ++index)
				memcpy(&pixelsGauss[index * 3], bg, 3);

			int setIndexStart = loopIndex;
			int setIndexEnd = loopIndex + 1;
			if (showSetPermutations)
			{
				setIndexStart = 0;
				setIndexEnd = numSets;
			}

			for (int setIndex = setIndexStart; setIndex < setIndexEnd; ++setIndex)
			{
				if (showSetPermutations && (loopIndex & (1 << setIndex)) == 0)
					continue;

				unsigned char fg[3];
				pointColoringObject.GetSetColor(setIndex, fg[0], fg[1], fg[2]);

				for (const float2& p : pointColoringObject.GetSet(points, setIndex))
				{
					if (c_drawPointImage)
					{
						int x = (int)Clamp((p.x * 0.5f + 0.5f) * float(c_pointImageSize - 1), 0.0f, float(c_pointImageSize - 1));
						int y = (int)Clamp((p.y * 0.5f + 0.5f) * float(c_pointImageSize - 1), 0.0f, float(c_pointImageSize - 1));
						memcpy(&pixels[(y * c_pointImageSize + x) * 3], fg, 3);
					}

					if (c_drawGaussImage)
					{
						int x = (int)Clamp((p.x * 0.5f + 0.5f) * float(c_pointImageGaussSize - 1), 0.0f, float(c_pointImageGaussSize - 1));
						int y = (int)Clamp((p.y * 0.5f + 0.5f) * float(c_pointImageGaussSize - 1), 0.0f, float(c_pointImageGaussSize - 1));
						PlotGaussian<3>(pixelsGauss, c_pointImageGaussSize, c_pointImageGaussSize, x, y, c_pointImageGaussBlobSigma, fg);
					}
				}
			}

			char permutationString[64] = { 0 };
			if (numSets > 1)
			{
				if (!showSetPermutations)
				{
					sprintf_s(permutationString, ".%i", imageCount);
					imageCount++;
				}
				else
				{
					strcat(permutationString, ".");
					for (int setIndex = 0; setIndex < numSets; ++setIndex)
					{
						if ((loopIndex & (1 << setIndex)) == 0)
							strcat(permutationString, "F");
						else
							strcat(permutationString, "T");
					}
				}
			}

			if (c_drawPointImage)
			{
				char fileName[1024];
				sprintf_s(fileName, "%s_%i_%i%s.png", baseFileName, progressIndex, progressTotal, permutationString);
				stbi_write_png(fileName, c_pointImageSize, c_pointImageSize, 3, pixels.data(), 0);
			}

			if (c_drawGaussImage)
			{
				char fileName[1024];
				sprintf_s(fileName, "%s_%i_%i%s.gauss.png", baseFileName, progressIndex, progressTotal, permutationString);
				stbi_write_png(fileName, c_pointImageGaussSize, c_pointImageGaussSize, 3, pixelsGauss.data(), 0);
			}
		}
	}
}

void GeneratePoints(int numPoints, int batchSize, int numIterations, const char* baseFileName, int numProgressImages, bool useGradientScalingFactor, bool stratifyLine, bool toroidalFixup, const std::vector<const FSOTClassBase*>& classes, const PointColoringObjectBase& pointColoringObject)
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
		if (numProgressImages > 0 && c_drawStartingState)
		{
			int progressInterval = numIterations / numProgressImages;
			if (iterationIndex % progressInterval == 0)
				SavePointSet(points, baseFileName, iterationIndex / progressInterval, numProgressImages, pointColoringObject);
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
			float subClassZ = 0.0f;
			{
				std::uniform_real<float> dist(0.0f, 1.0f);
				subClassZ = dist(rng);
			}

			// Filter the points so that batchData.sorted contains the point indices that are members
			FSOTClass.Filter(points, subClassZ, batchData.sorted);

			// if no points selected, don't do anything
			if (batchData.sorted.size() == 0)
				continue;

			// Project the points
			for (int i : batchData.sorted)
				batchData.projections[i] = Dot(direction, points[i]);

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
				if (stratifyLine)
					jitter = distJitter(rng);
				float projectionTarget = FSOTClass.ICDF((float(i) + jitter) / float(batchData.sorted.size()), direction);

				float projectionDiff = projectionTarget - projection;

				float gradientFactor = 1.0f;
				if (useGradientScalingFactor && batchData.sorted.size() > 1)
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

			if (toroidalFixup)
			{
				while (points[i].x < -1.0f)
					points[i].x += 2.0f;

				while (points[i].x > 1.0f)
					points[i].x -= 2.0f;

				while (points[i].y < -1.0f)
					points[i].y += 2.0f;

				while (points[i].y > 1.0f)
					points[i].y -= 2.0f;
			}

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
	SavePointSet(points, baseFileName, numProgressImages, numProgressImages, pointColoringObject);

	// report how long this took
	float elpasedSeconds = std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::high_resolution_clock::now() - start).count();
	printf("%0.2f seconds\n\n", elpasedSeconds);
}

int main(int argc, char** argv)
{
	_mkdir("out");

	// 1) Show the benefit of the gradient scaling fix
	// 2) Compare stratifying the line vs not
	// 3) Compare toroidal fixup vs not
	if (false)
	{
		FSOTClass<ICDF_UniformSquare, Filter_All> a;
		PointColoringObject_OneSetBlack pointColoringObject;

		// 1) Show the benefit of the gradient scaling fix
		GeneratePoints(1000, 64, 1000, "out/GradFixNoStratifyNoToroidalNo", 1, false, false, false, { &a }, pointColoringObject);
		GeneratePoints(1000, 64, 1000, "out/GradFixYesStratifyNoToroidalNo", 1, true, false, false, { &a }, pointColoringObject);

		// 2) Compare stratifying the line vs not
		GeneratePoints(1000, 64, 1000, "out/GradFixYesStratifyYesToroidalNo", 1, true, true, false, { &a }, pointColoringObject);

		// 3) Compare toroidal fixup vs not
		GeneratePoints(1000, 64, 1000, "out/GradFixYesStratifyNoToroidalYes", 1, true, false, true, { &a }, pointColoringObject);
	}

	// Multiclass
	{
		PointColoringObject_HalfRedHalfBlue pointColoringObject;

		FSOTClass<ICDF_UniformSquare, Filter_Range<0.0f, 0.5f>> firstHalf;
		FSOTClass<ICDF_UniformSquare, Filter_Range<0.5f, 1.0f>> secondHalf;
		FSOTClass<ICDF_UniformSquare, Filter_All> all;

		// Same settings as official example from their github, including 50% chance to do the combined set
		GeneratePoints(1024, 256, 4096, "out/Multiclass", 1, true, false, false, { &firstHalf, &secondHalf, &all, &all }, pointColoringObject);
	}

	// Progressive vs Non progressive
	if (false)
	{
		FSOTClass<ICDF_UniformSquare, Filter_Progressive> a;
		FSOTClass<ICDF_UniformSquare, Filter_All> b;

		PointColoringObject_Progressive<4> pointColoringObject;

		// Same settings as official example from their github
		GeneratePoints(1024, 128, 4096, "out/ProgressiveYes", 1, true, false, false, { &a }, pointColoringObject);
		GeneratePoints(1024, 128, 4096, "out/ProgressiveNo", 1, true, false, false, { &b }, pointColoringObject);
	}

	return 0;
}

/*

FSOT PAPER NOTES:
* FSOT paper has batch size of 64 by default. 4096 iterations. 4096 points.
* they make random points in space and project (not jitter, random!), instead of doing equal area points on the line. getInverseRadonNCube()

BLOG:
* first, compare GradFixNo to GradFixYes, showing how it improves that "overconvergence" thing. should have a better DFT than not letting it go to convergence (compare the 3 point sets, and DFTs)
* Note that you are making the target be the center of each bucket, put through the ICDF, but they are stratifying. Compare the 2, see if there's a difference. note that not timothy lottes mentioned that too, after the last post.
* show multiclass
* show progressive

TODO:
* could add a weight for each class, and do "weighted round robin" for selection or something? yeah do this
* projective point sets (and progressive / projective. and toroidally progressive / projective)
 * for projective may need to have the ICDF be able to define the projection direction (axis aligned!)
* your progressive point st isn't the highest quality. why not? They use adam (see slicedOptimalTransportBatchCube_progressive()), maybe that is why? but their code doesn't really even do progressive as far as i can tell...
 * same with multiclass right?
 * could generate and look at their points and see if they are the same
* could try showing DFTs that are the average of several realizations. may help show if stratification is good or not?
* you need to be able to explain the gradient scaling for your blog post. why does it improve things?
* point sets, and noise textures
* progressive point set.
* toroidally progressive point set (secret for now? JCGT?)
? what does a continuous membership even mean? maybe show 1 class with a box membership function vs a smooth membership function.
 * with fewer points (higher Z value) in continuous, those points get spread out more. weird.
? how do they support toroidal distance for their points? i did a toroidal fixup but that doesn't seem to be it?
* profile again

Toroidal Progressiveness
* compare vs halton, and R2
* projective blue noise and not projective
* figure out adam first and use it?

*/
