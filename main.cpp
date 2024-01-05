#include <vector>
#include <direct.h>
#include <stdio.h>

#include "squarecdf.h"
#include "NumericalCDF.h"

#include "FSOT.h"

#include "Filters.h"
#include "CDFs.h"
#include "DebugOutput.h"


void GenerateRandomFloat2(std::mt19937& rng, float2& p)
{
	std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
	p.x = dist(rng);
	p.y = dist(rng);
}

int main(int argc, char** argv)
{
	_mkdir("out");

	static const bool c_doGradFixComparison = false;
	static const bool c_doMulticlass = false;
	static const bool c_doProgressive = false;
	static const bool c_doProjective = true;

	static const bool c_do1DBNTexture = true;

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
		fsot.m_useGradientScalingFactor = false;
		fsot.m_stratifyLine = false;

		// 1) Show the benefit of the gradient scaling fix
		fsot.GeneratePoints<DebugOutput_Points_float2>(1024, "out/GradFixNoStratifyNoToroidalNo", 1, GenerateRandomFloat2);

		fsot.m_useGradientScalingFactor = true;
		fsot.GeneratePoints<DebugOutput_Points_float2>(1024, "out/GradFixYesStratifyNoToroidalNo", 1, GenerateRandomFloat2);

		// 2) Compare stratifying the line vs not
		fsot.m_stratifyLine = true;
		fsot.GeneratePoints<DebugOutput_Points_float2>(1024, "out/GradFixYesStratifyYesToroidalNo", 1, GenerateRandomFloat2);
		fsot.m_stratifyLine = false;

		/*
		// 3) Compare toroidal fixup vs not
		fsot.m_toroidalFixup = true;
		fsot.GeneratePoints<PointColoringObject_OneSetBlack>(1024, "out/GradFixYesStratifyNoToroidalYes", 1, GenerateRandomFloat2);
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
		fsot.m_useGradientScalingFactor = true;
		fsot.m_stratifyLine = false;

		fsot.GeneratePoints<DebugOutput_Points_float2_TwoClass>(1024, "out/Multiclass", 1, GenerateRandomFloat2);
	}

	// Progressive vs Non progressive
	if (c_doProgressive)
	{
		FSOT<float2> fsot_progressive;
		fsot_progressive.AddClass<CDF_UniformSquare, Filter_Progressive>();
		fsot_progressive.m_batchSize = 128;
		fsot_progressive.m_numIterations = 4096;
		fsot_progressive.m_useGradientScalingFactor = true;
		fsot_progressive.m_stratifyLine = false;

		FSOT<float2> fsot_not;
		fsot_not.AddClass<CDF_UniformSquare, Filter_All>();
		fsot_not.m_batchSize = 128;
		fsot_not.m_numIterations = 4096;
		fsot_not.m_useGradientScalingFactor = true;
		fsot_not.m_stratifyLine = false;

		fsot_progressive.GeneratePoints<DebugOutput_Points_float2_Progressive<4>>(1024, "out/ProgressiveYes", 1, GenerateRandomFloat2);
		fsot_not.GeneratePoints<DebugOutput_Points_float2_Progressive<4>>(1024, "out/ProgressiveNo", 1, GenerateRandomFloat2);
	}

#if 0

	// Projective
	if (c_doProjective)
	{
		static const int c_numPoints = 1000;

		PointColoringObject_OneSetBlack pointColoringObject;

		// Use weighting from the projective blue noise paper
		{
			static const float c_MaxRadius1D = MaxPackedSphereRadiusDimension1(c_numPoints);
			static const float c_MaxRadius2D = MaxPackedSphereRadiusDimension2(c_numPoints);

			static const float c_totalWeight = c_MaxRadius1D + c_MaxRadius1D + c_MaxRadius2D;

			FSOTClass<ICDF_UniformXAxis, Filter_All> a(c_MaxRadius1D / c_totalWeight);
			FSOTClass<ICDF_UniformYAxis, Filter_All> b(c_MaxRadius1D / c_totalWeight);
			FSOTClass<ICDF_UniformSquare, Filter_All> c(c_MaxRadius2D / c_totalWeight);

			// For 1000 points, these weights are: 0.028, 0.028, 0.94.
			// if we do the other way in the projective blue noise paper where it's 1D/2D and 2D/2D respectively, it's nearly the same:
			// 0.029, 0.029, 1. 

			GeneratePoints(c_numPoints, 64, 10000, "out/Projective_Paper", 1, true, false, false, { &a, &b, &c }, pointColoringObject);
		}

		// Use equal weighting for the projections, like the sliced optimal transport paper uses
		{
			FSOTClass<ICDF_UniformXAxis, Filter_All> a;
			FSOTClass<ICDF_UniformYAxis, Filter_All> b;
			FSOTClass<ICDF_UniformSquare, Filter_All> c;

			GeneratePoints(c_numPoints, 64, 10000, "out/Projective_Equal", 1, true, false, false, { &a, &b, &c }, pointColoringObject);
		}

		// use double weighting for the combined result
		{
			FSOTClass<ICDF_UniformXAxis, Filter_All> a;
			FSOTClass<ICDF_UniformYAxis, Filter_All> b;
			FSOTClass<ICDF_UniformSquare, Filter_All> c(2.0f);

			GeneratePoints(c_numPoints, 64, 10000, "out/Projective_SemiEqual", 1, true, false, false, { &a, &b, &c }, pointColoringObject);
		}
	}

	if (c_do1DBNTexture)
	{
		//FSOT<1> plan;

		//FSOTClass<ICDF_UniformSquare, Filter_Range<0.5f, 1.0f>> secondHalf;
	}

#endif

	return 0;
}

/*

TODO:
* finish cleaning up all this code
* filters could be a lambda. ICDF is more complex though
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
