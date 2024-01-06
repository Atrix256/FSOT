#pragma once

#include "maths.h"

class CDF_UniformSquare
{
public:
	bool OverrideDirection(float2& direction) const
	{
		return false;
	}

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

class CDF_UniformXAxis
{
public:
	bool OverrideDirection(float2& direction) const
	{
		direction.x = 1.0f;
		direction.y = 0.0f;
		return true;
	}

	float ICDF(float y, const float2& direction) const
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
	bool OverrideDirection(float2& direction) const
	{
		direction.x = 0.0f;
		direction.y = 1.0f;
		return true;
	}

	float ICDF(float y, const float2& direction) const
	{
		// Convert y: square is in [-0.5, 0.5], but y is in [0, 1].
		y = y - 0.5f;

		// Evaluate ICDF
		float x = y;

		// The CDF is in [-0.5, 0.5], but we want the points to be in [-1,1]
		return x * 2.0f;
	}
};

class CDF_Uniform1D
{
public:
	bool OverrideDirection(float& direction) const
	{
		direction = 1.0f;
		return true;
	}

	float ICDF(float y, const float& direction) const
	{
		return y;
	}
};
