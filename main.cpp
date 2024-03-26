#include <vector>
#include <direct.h>
#include <stdio.h>

#include "squarecdf.h"
#include "NumericalCDF.h"

#include "FSOT.h"

#include "Filters.h"
#include "CDFs.h"
#include "DebugOutput.h"

void GenerateRandomFloat(std::mt19937& rng, float& p)
{
	std::uniform_real_distribution<float> dist(0.0f, 1.0f);
	p = dist(rng);
}

void GenerateRandomDirectionFloat(std::mt19937& rng, float& direction)
{
	std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
	direction = dist(rng) >= 0.0f ? 1.0f : -1.0f;
}

float ProjectPoints1D(float p, float dir)
{
	return p * dir;
}

void GenerateRandomFloat2(std::mt19937& rng, float2& p)
{
	std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
	p.x = dist(rng);
	p.y = dist(rng);
}

void GenerateRandomDirectionFloat2(std::mt19937& rng, float2& direction)
{
	std::normal_distribution<float> dist(0.0f, 1.0f);
	direction.x = dist(rng);
	direction.y = dist(rng);
	direction = Normalize(direction);
}

int main(int argc, char** argv)
{
	_mkdir("out");

	static const bool c_doGradFixComparison = false;
	static const bool c_doMulticlass = false;
	static const bool c_doProgressive = false;
	static const bool c_doToroidalProgressiveBig = false;
	static const bool c_doToroidalProgressiveSmall = false;
	static const bool c_doProjective = false;

	static const bool c_do1DBNTexture = false;
	static const bool c_do2DBNTexture = true;

	// 1) Show the benefit of the gradient scaling fix
	// 2) Compare stratifying the line vs not
	// 3) Compare toroidal fixup vs not
	if (c_doGradFixComparison)
	{
		// Set up FSOT to make blue noise point sets
		FSOT<float2> fsot;
		fsot.AddClass<CDF_UniformSquare, Filter_All>();
		fsot.m_batchSize = 64;
		fsot.m_numIterations = 1000;

		// 1) Show the benefit of the gradient scaling fix
		fsot.m_useGradientScalingFactor = false;
		fsot.GeneratePoints<DebugOutput_Points_float2>(1024, "out/GradFixNoStratifyNoToroidalNo", 1, GenerateRandomFloat2, GenerateRandomDirectionFloat2, Dot);
		fsot.m_useGradientScalingFactor = true;

		fsot.GeneratePoints<DebugOutput_Points_float2>(1024, "out/GradFixYesStratifyNoToroidalNo", 1, GenerateRandomFloat2, GenerateRandomDirectionFloat2, Dot);

		// 2) Compare stratifying the line vs not
		fsot.m_stratifyLine = true;
		fsot.GeneratePoints<DebugOutput_Points_float2>(1024, "out/GradFixYesStratifyYesToroidalNo", 1, GenerateRandomFloat2, GenerateRandomDirectionFloat2, Dot);
		fsot.m_stratifyLine = false;

		/*
		// 3) Compare toroidal fixup vs not
		fsot.m_toroidalFixup = true;
		fsot.GeneratePoints<PointColoringObject_OneSetBlack>(1024, "out/GradFixYesStratifyNoToroidalYes", 1, GenerateRandomFloat2, GenerateRandomDirectionFloat2, Dot);
		fsot.m_toroidalFixup = false;
		*/
	}

	// Multiclass
	if (c_doMulticlass)
	{
		// Set up FSOT to make 2 class blue noise point sets
		FSOT<float2> fsot;
		fsot.AddClass<CDF_UniformSquare, Filter_Range<0.0f, 0.5f>>();
		fsot.AddClass<CDF_UniformSquare, Filter_Range<0.5f, 1.0f>>();
		fsot.AddClass<CDF_UniformSquare, Filter_All>(2.0f);
		fsot.m_batchSize = 256;
		fsot.m_numIterations = 4096;

		fsot.GeneratePoints<DebugOutput_Points_float2_TwoClass>(1024, "out/Multiclass", 1, GenerateRandomFloat2, GenerateRandomDirectionFloat2, Dot);
	}

	// Progressive vs Non progressive
	if (c_doProgressive)
	{
		FSOT<float2> fsot_progressive;
		fsot_progressive.AddClass<CDF_UniformSquare, Filter_Progressive>();
		fsot_progressive.m_batchSize = 128;
		fsot_progressive.m_numIterations = 4096;

		FSOT<float2> fsot_not;
		fsot_not.AddClass<CDF_UniformSquare, Filter_All>();
		fsot_not.m_batchSize = 128;
		fsot_not.m_numIterations = 4096;

		fsot_progressive.GeneratePoints<DebugOutput_Points_float2_Progressive<4, 1>>(1024, "out/ProgressiveYes", 1, GenerateRandomFloat2, GenerateRandomDirectionFloat2, Dot);
		fsot_not.GeneratePoints<DebugOutput_Points_float2_Progressive<4, 1>>(1024, "out/ProgressiveNo", 1, GenerateRandomFloat2, GenerateRandomDirectionFloat2, Dot);
	}

	// Toroidally Progressive - Many points
	if (c_doToroidalProgressiveBig)
	{
		static const int c_numPoints = 1024;

		FSOT<float2> fsot;
		fsot.m_batchSize = c_numPoints / 8;
		fsot.m_numIterations = 4096;

		for (int i = 0; i < c_numPoints; ++i)
		{
			auto& c = fsot.AddClass<CDF_UniformSquare, Filter_Progressive>();
			c.m_filter.m_firstIndex = i;
		}

		fsot.GeneratePoints<DebugOutput_Points_float2_Progressive<4, 4>>(c_numPoints, "out/TProgressiveBig", 1, GenerateRandomFloat2, GenerateRandomDirectionFloat2, Dot);
	}

	// Toroidally Progressive - Fewer points
	if (c_doToroidalProgressiveSmall)
	{
		static const int c_numPoints = 64;

		FSOT<float2> fsot;
		fsot.m_batchSize = c_numPoints / 8;
		fsot.m_numIterations = 4096;

		for (int i = 0; i < c_numPoints; ++i)
		{
			auto& c = fsot.AddClass<CDF_UniformSquare, Filter_Progressive>();
			c.m_filter.m_firstIndex = i;
		}

		fsot.GeneratePoints<DebugOutput_Points_float2_Progressive<4, 4>>(c_numPoints, "out/TProgressiveSmall", 1, GenerateRandomFloat2, GenerateRandomDirectionFloat2, Dot);
	}

	// Projective
	if (c_doProjective)
	{
		static const int c_numPoints = 1024;

		// Use weighting from the projective blue noise paper
		{
			static const float c_MaxRadius1D = MaxPackedSphereRadiusDimension1(c_numPoints);
			static const float c_MaxRadius2D = MaxPackedSphereRadiusDimension2(c_numPoints);

			static const float c_totalWeight = c_MaxRadius1D + c_MaxRadius1D + c_MaxRadius2D;

			FSOT<float2> fsot;
			fsot.AddClass<CDF_UniformXAxis, Filter_All>(c_MaxRadius1D / c_totalWeight);
			fsot.AddClass<CDF_UniformYAxis, Filter_All>(c_MaxRadius1D / c_totalWeight);
			fsot.AddClass<CDF_UniformSquare, Filter_All>(c_MaxRadius2D / c_totalWeight);
			fsot.m_batchSize = 64;
			fsot.m_numIterations = 10000;

			// For 1000 points, these weights are: 0.028, 0.028, 0.94.
			// if we do the other way in the projective blue noise paper where it's 1D/2D and 2D/2D respectively, it's nearly the same:
			// 0.029, 0.029, 1. 

			fsot.GeneratePoints<DebugOutput_Points_float2>(c_numPoints, "out/Projective_Paper", 1, GenerateRandomFloat2, GenerateRandomDirectionFloat2, Dot);
		}

		// Use equal weighting for the projections, like the sliced optimal transport paper uses
		{
			FSOT<float2> fsot;
			fsot.AddClass<CDF_UniformXAxis, Filter_All>();
			fsot.AddClass<CDF_UniformYAxis, Filter_All>();
			fsot.AddClass<CDF_UniformSquare, Filter_All>();
			fsot.m_batchSize = 64;
			fsot.m_numIterations = 10000;

			fsot.GeneratePoints<DebugOutput_Points_float2>(c_numPoints, "out/Projective_Equal", 1, GenerateRandomFloat2, GenerateRandomDirectionFloat2, Dot);
		}

		// use double weighting for the combined result
		{
			FSOT<float2> fsot;
			fsot.AddClass<CDF_UniformXAxis, Filter_All>();
			fsot.AddClass<CDF_UniformYAxis, Filter_All>();
			fsot.AddClass<CDF_UniformSquare, Filter_All>(2.0f);
			fsot.m_batchSize = 64;
			fsot.m_numIterations = 10000;

			fsot.GeneratePoints<DebugOutput_Points_float2>(c_numPoints, "out/Projective_SemiEqual", 1, GenerateRandomFloat2, GenerateRandomDirectionFloat2, Dot);
		}
	}

	if (c_do1DBNTexture)
	{
		static const int c_numPixels = 1024;
		static const int c_filterSize = 4;

		FSOT<float> fsot;
		fsot.m_batchSize = 512;
		fsot.m_numIterations = 10000;

		// This promotes a flat histogram across the texture
		fsot.AddClass<CDF_Uniform1D, Filter_All>();

		// These promote diversity in small regions of pixels
		for (int i = 0; i < c_numPixels; ++i)
		{
			auto& c = fsot.AddClass<CDF_Uniform1D, Filter_Range1D>();
			c.m_filter.m_index = i;
			c.m_filter.m_radius = c_filterSize;
		}

		fsot.GeneratePoints<DebugOutput_Texture1D_float>(c_numPixels, "out/Texture1D", 1, GenerateRandomFloat, GenerateRandomDirectionFloat, ProjectPoints1D);

		// TODO: tune the settings
		// TODO: should we jitter?
		// TODO: compare with flat histogram class, and not. flat histogram thing makes variety i think, instead of everything going flat
		// TODO: also, box filter vs gauss!
		// TODO: should this be simplified for 1d textures? don't need random projections, only random selections and no repeats? or maybe ok as is.
		//   Or maybe we need an option for shuffling then choosing the first N instead of randomly selecting.
	}

	if (c_do2DBNTexture)
	{
		static const int c_imageSize = 64;
		static const int c_numPixels = c_imageSize * c_imageSize;

		static const int c_filterSizes[] = { 1, 2, 3 };

		for (int filterSize : c_filterSizes)
		{
			FSOT<float> fsot;
			fsot.m_batchSize = 512;
			fsot.m_numIterations = 10000;

			// This promotes a flat histogram across the texture
			fsot.AddClass<CDF_Uniform1D, Filter_All>();

			// These promote diversity in small regions of pixels
			for (int i = 0; i < c_numPixels; ++i)
			{
				auto& c = fsot.AddClass<CDF_Uniform1D, Filter_Range2D>();
				c.m_filter.m_index = i;
				c.m_filter.m_radius = filterSize;
				c.m_filter.m_width = c_imageSize;
			}

			char fileNameBase[256];
			sprintf_s(fileNameBase, "out/Texture2D_f%i", filterSize);

			fsot.GeneratePoints<DebugOutput_Texture2D_float>(c_numPixels, fileNameBase, 1, GenerateRandomFloat, GenerateRandomDirectionFloat, ProjectPoints1D);
		}

		// TODO: finish this
	}


	return 0;
}

/*
TODO: try Adam!

TODO:
* make "DebugOutput_Texture2D_float::Output" spit out 2d arrays for txt. and csv?
* maybe combine DebugOutput_Texture1D_float and DebugOutput_Texture2D_float into one thing
* the debug output objects need to let the user change values on them, so pass an object to generate points instead of a type.
 * like DebugOutput_Texture2D_float.m_width!
* for toroidally progressive, we could also try and add a "all points are blue together" constraint?
* are you biasing towards the first class in random selection? your toroidally progressive point sets seem like that might be true?
* the code saving progressive point sets as text like for a header file just saves the same points over and over again.
* clean up code and make a "2d" version for textures.
 * ! start by making a 1d texture. all it is, is "positional membership", with a bunch of overlapping sets of membership.
  * it makes diversity in small regions.
 * could actually make textures with this same 1d code. unwrap the image from 2d to 1d. generalizes to ND
 * look at DFT
 * also histogram
 * and threshold
 * compare box vs gauss -> both for membership and for ICDF?
 ? can we make spatiotemporal blue noise? i do think so!
 ? can we make projective (scalar) blue noise textures?
* toroidally progressive point set (secret for now? JCGT?)
 * and projective while toroidally progressives, after
* noise textures
* your progressive point st isn't the highest quality. why not? They use adam (see slicedOptimalTransportBatchCube_progressive()), maybe that is why? but their code doesn't really even do progressive as far as i can tell...
 * same with multiclass right?
 * could generate and look at their points and see if they are the same
 * perhaps same with projective noise. I don't know that you have the right weighting yet
* could try showing DFTs that are the average of several realizations. may help show if stratification is good or not? do it as one off special work in that specialized python script we haven't pulled over
* you need to be able to explain the gradient scaling for your blog post. why does it improve things?
? how do they support toroidal distance for their points? i did a toroidal fixup but that doesn't seem to be it?
* profile again
* clean up this code, like how you define classes and stuff.
* need to figure out adam?
? OT tends to make blue noise. is there a way to make it red etc? maybe using LPF (etc) filters as a pdf?
 * could try a box filter instead of a gaussian? oh wait that's just a uniform distribution.
* your projective noise with balanced weights isn't giving as good a power spectrum as the SOT paper does.
 * see page 9: http://www.geometry.caltech.edu/pubs/PBCIW+20.pdf
 * maybe try their code if they have any
? should we do something with numerical CDF (images?). might look better after adam is sorted out!

Toroidal Progressiveness
* compare vs halton, and R2
* projective blue noise and not projective
* figure out adam first and use it?




FSOT PAPER NOTES:
* FSOT paper has batch size of 64 by default. 4096 iterations. 4096 points.
* they make random points in space and project (not jitter, random!), instead of doing equal area points on the line. getInverseRadonNCube()

BLOG:
* first, compare GradFixNo to GradFixYes, showing how it improves that "overconvergence" thing. should have a better DFT than not letting it go to convergence (compare the 3 point sets, and DFTs)
* Note that you are making the target be the center of each bucket, put through the ICDF, but they are stratifying. Compare the 2, see if there's a difference. note that not timothy lottes mentioned that too, after the last post.
* show multiclass
 * mention that you added a weighting to each class
* show progressive
* talk about classes and subclasses, and selection of them
 * continuous subclasses -> the points move slower by being part of a smaller group (m/n) but they also move as part of the larger subset. BUT they move to a different location. kind of strange, not sure if that's good behavior
* show projective blue noise
 ! projective blue noise paper: https://resources.mpi-inf.mpg.de/ProjectiveBlueNoise/ProjectiveBlueNoise.pdf
 * talk about how you balanced the weights, like the paper says to.
 * section 3.6 says they just use equal weight for all of em, and didn't do much analysis otherwise. http://www.geometry.caltech.edu/pubs/PBCIW+20.pdf
! this is a great framework.  Hard to know how many batches to do.  Hard to know how many iterations to do, or a learning rate, but Adam helps.

*/
