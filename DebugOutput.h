#pragma once

#include "maths.h"
#include <vector>

class DebugOutput_SetInfoBase
{
public:
private:
};

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

	DebugOutput_SetInfoBase* m_setInfo = nullptr;

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