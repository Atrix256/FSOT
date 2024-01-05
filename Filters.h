#pragma once

class Filter_All
{
public:
	template <typename T>
	static void Filter(const std::vector<T>& points, float subClassZ, std::vector<int>& memberPoints)
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
	template <typename T>
	static void Filter(const std::vector<T>& points, float subClassZ, std::vector<int>& memberPoints)
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
	template <typename T>
	static void Filter(const std::vector<T>& points, float subClassZ, std::vector<int>& memberPoints)
	{
		int indexMax = Clamp(int(subClassZ * float(points.size())), 0, int(points.size()));
		memberPoints.resize(indexMax);
		for (int i = 0; i < indexMax; ++i)
			memberPoints[i] = i;
	}
};
