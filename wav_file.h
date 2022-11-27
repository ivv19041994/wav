#pragma once

#include <string>
#include <map>
#include <vector>
#include <fstream>

class WavFile {
public:
	WavFile(const std::string& path);

	WavFile intersection(const WavFile& other) const;

	double SimilarIntersection(const WavFile& other) const;

	void save(const std::string& path) const;
	void SetNumChannels(size_t numChannels);

private:

	enum audio_format {
		pcm = 1
	};

#pragma pack(push, 1)
	struct fmt {
		audio_format audioFormat:16;
		uint16_t numChannels;
		uint32_t sampleRate;
		uint32_t byteRate;
		uint16_t blockAlign;
		uint16_t bitsPerSample;
	};
#pragma pack(pop)

	std::string chunkId;
	std::string format;

	std::map<std::string, std::vector<unsigned char>> subchanks;
	struct fmt fmt;
	std::vector<std::vector<double>> data;
	std::vector<int16_t> data_debug;

	std::pair<long long, double> CorrelationMax(const WavFile& other) const;

	double correlation(const std::vector<double>& f1, const std::vector<double>& f2, long long x) const;
	size_t CalculateChankSize() const;

	void FillSubchanksByFile(const std::string& path);

	uint16_t ReadWord(std::vector<unsigned char>::iterator it);
	uint32_t ReadDWord(std::vector<unsigned char>::iterator it);

	void PutWord(std::basic_fstream <unsigned char>& fs, uint16_t word) const;
	void PutWord(std::basic_fstream <unsigned char>& fs, int16_t word) const;
	void PutDWord(std::basic_fstream <unsigned char>& fs, uint32_t dword) const;
	void PutDWord(std::basic_fstream <unsigned char>& fs, size_t dword) const;

	void WriteChankHeader(std::basic_fstream <unsigned char>& fs) const;
	void WriteFmtSubchank(std::basic_fstream <unsigned char>& fs) const;
	void WriteDataSubchank(std::basic_fstream <unsigned char>& fs) const;
	void WriteOtherSubchanks(std::basic_fstream <unsigned char>& fs) const;

	void FillFmtSubchank();
	void FillDataSubchank();

	void NormalizeDataChannels(double max_coefficient);




};

