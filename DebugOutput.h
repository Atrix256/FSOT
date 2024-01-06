#pragma once

#include "maths.h"
#include <vector>

class DebugOutput_Points_float2
{
public:
	void Output(const std::vector<float2>& points, const char* baseFileName, int iteration, int iterationCount);

	bool m_reportStartingState = false;
	int m_numProgressResults = 1;

	bool m_writeTxt = true;
	bool m_writeCsv = true;

	bool m_makePointImage = true;
	int m_pointImageSize = 256;

	bool m_makeGaussImage = true;
	int m_gaussImageSize = 512;
	float m_gaussImageSigma = 1.5f;

protected:
	virtual int NumSets(const std::vector<float2>& points) const
	{
		return 1;
	}

	virtual std::vector<float2> GetSet(const std::vector<float2>& points, int setIndex) const
	{
		return points;
	}

	virtual void GetSetColor(int setIndex, unsigned char& R, unsigned char& G, unsigned char& B) const
	{
		R = 0;
		G = 0;
		B = 0;
	}

	virtual bool ShowSetPermutations() const
	{
		return true;
	}

private:
	void MakePointsCSV(const std::vector<float2>& points, const char* fileName);
	void MakePointsText(const std::vector<float2>& points, const char* fileName);
	void MakePointImage(const std::vector<float2>& points, const char* baseFileName);
	void MakeGaussImage(const std::vector<float2>& points, const char* baseFileName);
};

class DebugOutput_Points_float2_TwoClass : public DebugOutput_Points_float2
{
protected:
	virtual int NumSets(const std::vector<float2>& points) const
	{
		return 2;
	}

	virtual std::vector<float2> GetSet(const std::vector<float2>& points, int setIndex) const override final
	{
		int start = (setIndex == 0) ? 0 : int(points.size()) / 2;
		int end = (setIndex == 0) ? int(points.size()) / 2 : int(points.size());

		std::vector<float2> ret(end - start);
		for (int i = start; i < end; ++i)
			ret[i - start] = points[i];

		return ret;
	}

	virtual void GetSetColor(int setIndex, unsigned char& R, unsigned char& G, unsigned char& B) const override final
	{
		R = (setIndex == 0) ? 255 : 0;
		G = 0;
		B = (setIndex == 1) ? 255 : 0;
	}

	virtual bool ShowSetPermutations() const override final
	{
		return true;
	}
};

template <int PROGRESSION_COUNT, int NUM_PROGRESSION_OFFSETS>
class DebugOutput_Points_float2_Progressive : public DebugOutput_Points_float2
{
protected:
	virtual int NumSets(const std::vector<float2>& points) const
	{
		return PROGRESSION_COUNT * NUM_PROGRESSION_OFFSETS;
	}

	virtual std::vector<float2> GetSet(const std::vector<float2>& points, int setIndex) const override final
	{
		int progressionIndex = setIndex % PROGRESSION_COUNT;
		int progressionOffsetIndex = setIndex / PROGRESSION_COUNT;

		int setOffset = int(float(points.size()) * float(progressionOffsetIndex) / float(NUM_PROGRESSION_OFFSETS));

		int setEnd = int(float(progressionIndex + 1) * float(points.size()) / float(PROGRESSION_COUNT));
		std::vector<float2> ret(setEnd);
		for (int i = 0; i < setEnd; ++i)
			ret[i] = points[(i + setOffset) % points.size()];
		return ret;
	}

	virtual void GetSetColor(int setIndex, unsigned char& R, unsigned char& G, unsigned char& B) const override final
	{
		R = 0;
		G = 0;
		B = 0;
	}

	virtual bool ShowSetPermutations() const override final
	{
		return false;
	}
};

class DebugOutput_Texture1D_float
{
public:
	void Output(const std::vector<float>& points, const char* baseFileName, int iteration, int iterationCount);

	bool m_reportStartingState = true;
	int m_numProgressResults = 1;

	bool m_writeTxt = true;
	bool m_writeCsv = true;

	bool m_makeImage = true;
	int m_imageSize = 256;
};

class DebugOutput_Texture2D_float
{
public:
	void Output(const std::vector<float>& points, const char* baseFileName, int iteration, int iterationCount);

	bool m_reportStartingState = true;
	int m_numProgressResults = 1;

	bool m_writeTxt = true;
	bool m_writeCsv = true;

	bool m_makeImage = true;
	int m_imageSize = 256;

	// TODO: this needs to be able to be set by the people calling into FSOT. maybe take an object instead of a type.
	int m_width = 64;
};
