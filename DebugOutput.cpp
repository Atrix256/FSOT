#define _CRT_SECURE_NO_WARNINGS // for stb

#include "DebugOutput.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

template <int NumChannels>
static void PlotGaussian(std::vector<unsigned char>& image, int width, int height, int x, int y, float sigma, unsigned char color[NumChannels])
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

void DebugOutput_Points_float2::Output(const std::vector<float2>& points, const char* baseFileName, int iteration, int iterationCount)
{
	// Don't report progress if we shouldn't
	if (m_numProgressResults == 0)
		return;

	int progressInterval = iterationCount / m_numProgressResults;
	int progressIndex = (iteration == iterationCount) ? m_numProgressResults : iteration / progressInterval;
	if (iteration != iterationCount && (iteration % progressInterval) != 0)
		return;

	if (iteration == 0 && !m_reportStartingState)
		return;

	char fileName[1024];

	if (m_writeTxt)
	{
		sprintf_s(fileName, "%s_%i_%i.txt", baseFileName, progressIndex, m_numProgressResults);
		MakePointsText(points, fileName);
	}

	if (m_writeCsv)
	{
		sprintf_s(fileName, "%s_%i_%i.csv", baseFileName, progressIndex, m_numProgressResults);
		MakePointsCSV(points, fileName);

	}

	if (m_makePointImage)
	{
		sprintf_s(fileName, "%s_%i_%i", baseFileName, progressIndex, m_numProgressResults);
		MakePointImage(points, fileName);
	}

	if (m_makeGaussImage)
	{
		sprintf_s(fileName, "%s_%i_%i", baseFileName, progressIndex, m_numProgressResults);
		MakeGaussImage(points, fileName);
	}
}

void DebugOutput_Points_float2::MakePointsCSV(const std::vector<float2>& points, const char* fileName)
{
	const int numSets = NumSets(points);

	FILE* file = nullptr;
	fopen_s(&file, fileName, "wb");

	fprintf(file, "\"X\",\"Y\",\"Class\"\n");

	for (int setIndex = 0; setIndex < numSets; ++setIndex)
	{
		for (const float2& p : GetSet(points, setIndex))
			fprintf(file, "\"%f\",\"%f\",\"%i\"\n", p.x, p.y, setIndex);
	}

	fclose(file);
}

void DebugOutput_Points_float2::MakePointsText(const std::vector<float2>& points, const char* fileName)
{
	const int numSets = NumSets(points);

	FILE* file = nullptr;
	fopen_s(&file, fileName, "wb");

	if (numSets > 1)
	{
		fprintf(file, "float points[][3] =\n{\n");

		for (int setIndex = 0; setIndex < numSets; ++setIndex)
		{
			fprintf(file, "\n    // Class %i\n", setIndex);
			for (const float2& p : GetSet(points, setIndex))
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

void DebugOutput_Points_float2::MakePointImage(const std::vector<float2>& points, const char* baseFileName)
{
	const int numSets = NumSets(points);

	unsigned char bg[3] = { 255, 255, 255 };

	std::vector<unsigned char> pixels(m_pointImageSize * m_pointImageSize * 3);

	bool showSetPermutations = ShowSetPermutations();
	int loopIndexStart = 0;
	int loopIndexEnd = numSets;
	if (showSetPermutations)
	{
		loopIndexStart = 1;
		loopIndexEnd = 1 << numSets;
	}

	int imageCount = 0;
	for (int loopIndex = loopIndexStart; loopIndex < loopIndexEnd; ++loopIndex)
	{
		for (size_t index = 0; index < m_pointImageSize * m_pointImageSize; ++index)
			memcpy(&pixels[index * 3], bg, 3);

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
			GetSetColor(setIndex, fg[0], fg[1], fg[2]);

			for (const float2& p : GetSet(points, setIndex))
			{
				int x = (int)Clamp((p.x * 0.5f + 0.5f) * float(m_pointImageSize - 1), 0.0f, float(m_pointImageSize - 1));
				int y = (int)Clamp((p.y * 0.5f + 0.5f) * float(m_pointImageSize - 1), 0.0f, float(m_pointImageSize - 1));
				memcpy(&pixels[(y * m_pointImageSize + x) * 3], fg, 3);
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

		char fileName[1024];
		sprintf_s(fileName, "%s%s.points.png", baseFileName, permutationString);
		stbi_write_png(fileName, m_pointImageSize, m_pointImageSize, 3, pixels.data(), 0);
	}
}

void DebugOutput_Points_float2::MakeGaussImage(const std::vector<float2>& points, const char* baseFileName)
{
	const int numSets = NumSets(points);

	unsigned char bg[3] = { 255, 255, 255 };

	std::vector<unsigned char> pixels(m_gaussImageSize * m_gaussImageSize * 3);

	bool showSetPermutations = ShowSetPermutations();
	int loopIndexStart = 0;
	int loopIndexEnd = numSets;
	if (showSetPermutations)
	{
		loopIndexStart = 1;
		loopIndexEnd = 1 << numSets;
	}

	int imageCount = 0;
	for (int loopIndex = loopIndexStart; loopIndex < loopIndexEnd; ++loopIndex)
	{
		for (size_t index = 0; index < m_gaussImageSize * m_gaussImageSize; ++index)
			memcpy(&pixels[index * 3], bg, 3);

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
			GetSetColor(setIndex, fg[0], fg[1], fg[2]);

			for (const float2& p : GetSet(points, setIndex))
			{
				int x = (int)Clamp((p.x * 0.5f + 0.5f) * float(m_gaussImageSize - 1), 0.0f, float(m_gaussImageSize - 1));
				int y = (int)Clamp((p.y * 0.5f + 0.5f) * float(m_gaussImageSize - 1), 0.0f, float(m_gaussImageSize - 1));
				PlotGaussian<3>(pixels, m_gaussImageSize, m_gaussImageSize, x, y, m_gaussImageSigma, fg);
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

		char fileName[1024];
		sprintf_s(fileName, "%s%s.gauss.png", baseFileName, permutationString);
		stbi_write_png(fileName, m_gaussImageSize, m_gaussImageSize, 3, pixels.data(), 0);
	}
}
