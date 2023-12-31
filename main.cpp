#include "FSOTClass.h"
#include <random>
#include <vector>

// Settings
#define DETERMINISTIC() true

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

int main(int argc, char** argv)
{
	static const int c_numPoints = 32;
	static const int c_batchSize = 64;
	static const int c_iterationCount = 1024;

	for (int iterationIndex = 0; iterationIndex < c_iterationCount; ++iterationIndex)
	{
		for (int batchIndex = 0; batchIndex < c_batchSize; ++batchIndex)
		{
			std::mt19937 rng = GetRNG(iterationIndex * c_batchSize + batchIndex);


		}
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

TODO:
! you need to figure out what the gradient scaling is all about. it improves quality, fixes that "overconvergence" thing.
 * there are 2 scaling factors and the code may not match the paper in labels.
 * Code says gradient scaling is the number of points in the subclass, divided by the total number of points. This is multiplied into the adjustment.
 * Code says other scaling is the next (+1) location minus the last (-1) location divided by 2. Multiplied by number of points in subclass. The adjustment is divided by this amount.
 * Paper says that gradient scaling is the average bin length divided by the current bin length. Seems to disagree with code for naming!
* omp across batches
* point sets, and textures
* progressive point set.
* toroidally progressive point set (secret for now? JCGT)
* write out csv of convergence, then can graph it
? what does a continuous membership even mean? maybe show 1 class with a box membership function vs a smooth membership function.
 * with fewer points (higher Z value) in continuous, those points get spread out more. weird.

Toroidal Progressiveness
* compare vs halton, and R2
* projective blue noise and not projective

*/
