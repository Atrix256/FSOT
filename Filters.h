#pragma once

class Filter_All
{
public:
	template <typename T>
	void Filter(const std::vector<T>& points, float subClassZ, std::vector<int>& memberPoints) const
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
	void Filter(const std::vector<T>& points, float subClassZ, std::vector<int>& memberPoints) const
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
	void Filter(const std::vector<T>& points, float subClassZ, std::vector<int>& memberPoints) const
	{
		int indexMax = Clamp(int(subClassZ * float(points.size())), 0, int(points.size()));
		memberPoints.resize(indexMax);
		for (int i = 0; i < indexMax; ++i)
			memberPoints[i] = i;// int((i + m_firstIndex) % points.size());
	}

	int m_firstIndex = 0;
};

class Filter_Range1D
{
public:
	template <typename T>
	void Filter(const std::vector<T>& points, float subClassZ, std::vector<int>& memberPoints) const
	{
		int centerIndex = m_index;

		memberPoints.clear();
		for (int i = -m_radius; i <= m_radius; ++i)
			memberPoints.push_back(int((centerIndex + i + points.size()) % points.size()));
	}

	int m_index = 0;
	int m_radius = 0;
};

class Filter_Range2D
{
public:
	template <typename T>
	void Filter(const std::vector<T>& points, float subClassZ, std::vector<int>& memberPoints) const
	{
		int centerIndex = m_index;
		int height = int(points.size()) / m_width;

		int centerX = centerIndex % m_width;
		int centerY = centerIndex / m_width;

		memberPoints.clear();

		for (int iy = -m_radius; iy <= m_radius; ++iy)
		{
			int py = (centerY + iy + height) % height;

			for (int ix = -m_radius; ix <= m_radius; ++ix)
			{
				int px = (centerX + ix + m_width) % m_width;

				memberPoints.push_back(py * m_width + px);
			}
		}
	}

	int m_index = 0;
	int m_radius = 0;
	int m_width = 1;
};
