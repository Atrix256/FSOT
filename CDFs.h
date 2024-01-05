#pragma once

#include "maths.h"

class CDF_UniformSquare
{
public:
	static bool OverrideDirection(float2& direction)
	{
		return false;
	}

	static float ICDF(float y, const float2& direction)
	{
		// Convert y: square is in [-0.5, 0.5], but y is in [0, 1].
		y = y - 0.5f;

		// Evaluate ICDF
		float x = Square::InverseCDF(y, direction);

		// The CDF is in [-0.5, 0.5], but we want the points to be in [-1,1]
		return x * 2.0f;
	}
};

class CDF_UniformXAxis
{
public:
	static bool OverrideDirection(float2& direction)
	{
		direction.x = 1.0f;
		direction.y = 0.0f;
		return true;
	}

	static float ICDF(float y, const float2& direction)
	{
		// Convert y: square is in [-0.5, 0.5], but y is in [0, 1].
		y = y - 0.5f;

		// Evaluate ICDF
		float x = y;

		// The CDF is in [-0.5, 0.5], but we want the points to be in [-1,1]
		return x * 2.0f;
	}
};

class CDF_UniformYAxis
{
public:
	static bool OverrideDirection(float2& direction)
	{
		direction.x = 0.0f;
		direction.y = 1.0f;
		return true;
	}

	static float ICDF(float y, const float2& direction)
	{
		// Convert y: square is in [-0.5, 0.5], but y is in [0, 1].
		y = y - 0.5f;

		// Evaluate ICDF
		float x = y;

		// The CDF is in [-0.5, 0.5], but we want the points to be in [-1,1]
		return x * 2.0f;
	}
};
