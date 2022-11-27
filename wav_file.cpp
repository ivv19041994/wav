#include "wav_file.h"

#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <execution>
#include <cmath>

#include <iostream>

using namespace std;

void WavFile::FillSubchanksByFile(const std::string& path) {
	std::basic_fstream <unsigned char> file(path, ios_base::binary | ios_base::in);

	if (!file) {
		throw runtime_error("can't open "s + path);
	}

	vector<unsigned char> readed;
	while (!file.eof()) {
		array<unsigned char, 1024> temp;
		file.read(temp.data(), 1024);
		readed.insert(readed.end(), temp.begin(), temp.begin() + file.gcount());
	}

	//cout << "readed size = " << readed.size() << endl;

	chunkId = string(readed.begin(), readed.begin() + 4);
	//cout << chunkId << endl;

	if (chunkId != "RIFF") {
		throw runtime_error("Invalid chank format: "s + chunkId);
	}

	{
		size_t chankSize = readed[4] | readed[5] << 8 | readed[6] << 16 | readed[7] << 24;
		//cout << chankSize << endl;
		if (readed.size() != chankSize + 8) {
			throw runtime_error("Invalid chank size format: "s + to_string(readed.size()) + " != "s + to_string(chankSize) + " + 8"s);
		}
	}

	format = string(readed.begin() + 8, readed.begin() + 12);
	//cout << format << endl;
	if (format != "WAVE") {
		throw runtime_error("Invalid format in header: "s + format);
	}

	size_t offset = 12;
	string sunchank_name;
	for (sunchank_name = { readed.begin() + offset, readed.begin() + offset + 4 };
		/*sunchank_name != "data" && (offset + 8) < readed.size()*/;
		sunchank_name = { readed.begin() + offset, readed.begin() + offset + 4 }) {
		offset += 4;
		size_t subchank_size = readed[offset] | readed[offset + 1] << 8 | readed[offset + 2] << 16 | readed[offset + 3] << 24;
		offset += 4;
		if (subchank_size > readed.size() - offset) {
			throw runtime_error("Invalid subchank size in header: "s + sunchank_name);
		}
		subchanks[sunchank_name] = { readed.begin() + offset, readed.begin() + offset + subchank_size };
		offset += subchank_size;
		if (offset + 8 >= readed.size()) {
			break;
		}
	}

}

uint16_t WavFile::ReadWord(vector<unsigned char>::iterator it) {
	return it[0] | it[1] << 8;
}

uint32_t WavFile::ReadDWord(vector<unsigned char>::iterator it) {
	return it[0] | it[1] << 8 | it[2] << 16 | it[3] << 24;
}

void WavFile::FillFmtSubchank() {
	if (subchanks.count("fmt ") == 0) {
		throw runtime_error("Invalid header: subchank \"fmt \" is not exist"s);
	}

	vector<unsigned char>& fmt_subchank = subchanks["fmt "];

	if (fmt_subchank.size() != 16) {
		if (fmt_subchank.size() >= 2) {
			throw runtime_error("Invalid header: unsupported format = "s + to_string(ReadWord(fmt_subchank.begin())));
		}
		else {
			throw runtime_error("Invalid header: unsupported format"s);
		}
	}

	fmt.audioFormat = static_cast<enum audio_format>(ReadWord(fmt_subchank.begin()));
	fmt.numChannels = ReadWord(fmt_subchank.begin() + 2);
	fmt.sampleRate = ReadDWord(fmt_subchank.begin() + 4);
	fmt.byteRate = ReadDWord(fmt_subchank.begin() + 8);
	fmt.blockAlign = ReadWord(fmt_subchank.begin() + 12);
	fmt.bitsPerSample = ReadWord(fmt_subchank.begin() + 14);

	subchanks.erase("fmt ");

	//cout
	//	<< "audioFormat     = " << fmt.audioFormat << endl
	//	<< "numChannels     = " << fmt.numChannels << endl
	//	<< "sampleRate      = " << fmt.sampleRate << endl
	//	<< "byteRate        = " << fmt.byteRate << endl
	//	<< "blockAlign      = " << fmt.blockAlign << endl
	//	<< "bitsPerSample   = " << fmt.bitsPerSample << endl;

	if (fmt.audioFormat != pcm) {
		throw runtime_error("Invalid header: unsupported format = "s + to_string(fmt.audioFormat));
	}
}

void WavFile::FillDataSubchank() {
	if (subchanks.count("data") == 0) {
		throw runtime_error("Invalid file: haven't data subchank "s);
	}

	size_t data_size = subchanks["data"].size();
	//cout << "data_size = " << data_size << endl;

	data.resize(fmt.numChannels);

	{
		if (fmt.bitsPerSample % 8) {
			throw runtime_error("Unsuppotred file: bits per sample not briefly 8. bits per sample = "s + to_string(fmt.bitsPerSample));
		}

		size_t byte_per_sample = fmt.bitsPerSample / 8;
		if (byte_per_sample > 2) {
			throw runtime_error("Unsuppotred file: bits per sample = "s + to_string(fmt.bitsPerSample));
		}
		for (auto it = subchanks["data"].begin(); it != subchanks["data"].end(); it += fmt.blockAlign) {
			auto cur = it;
			for (int channel = 0; channel < fmt.numChannels; ++channel) {
				int16_t sample = 0;
				int offset = 0;
				for (size_t byte = 0; byte < byte_per_sample; ++byte, ++cur) {
					sample |= (*cur) << offset;
					offset += 8;
				}
				data[channel].push_back(sample);
			}
		}
	}
	subchanks.erase("data");

	//cout << "data.size() = " << data.size() << " data[0].size() = " << data[0].size() << endl;
}

void WavFile::NormalizeDataChannels(double max_coefficient) {
	vector<double> max_vector(data.size());
	transform(execution::seq, data.begin(), data.end(), max_vector.begin(), [](const auto& sub_vec) {
		return *max_element(sub_vec.begin(), sub_vec.end());
		}
	);
	double m = *max_element(max_vector.begin(), max_vector.end());
	//cout << m << endl;
	m *= max_coefficient;

	for (auto& v : data) {
		transform(execution::seq, v.begin(), v.end(), v.begin(), [m](const auto& el) {
			return el / m;
			}
		);
	}
}


WavFile::WavFile(const std::string& path)
{
	FillSubchanksByFile(path);

	FillFmtSubchank();

	FillDataSubchank();

	NormalizeDataChannels(1.2);
}

double WavFile::SimilarIntersection(const WavFile& other) const {
	if (fmt.sampleRate != other.fmt.sampleRate) {
		throw runtime_error("Unsuppotred intersection for dif sample rates"s);
	}

	if (fmt.audioFormat != other.fmt.audioFormat) {
		throw runtime_error("Unsuppotred intersection for dif audio format"s);
	}

	auto [corr_offset, corr] = CorrelationMax(other);//на сколько нужно сдвинуть this объект, что бы получить пересечение с other

	long long common_size_wav;

	
	vector<double>::const_iterator it1, it2;
	if (corr_offset < 0) {
		//если отрицательный, то двигать начало надо у этого файла
		common_size_wav = min(static_cast<long long>(data[0].size()) + corr_offset, static_cast<long long>(other.data[0].size()));
		long long start_index = -corr_offset;
		it1 = data[0].begin() + start_index;
		it2 = other.data[0].begin();
	}
	else {
		common_size_wav = min(static_cast<long long>(other.data[0].size()) - corr_offset, static_cast<long long>(data[0].size()));
		it1 = data[0].begin();
		it2 = other.data[0].begin() + corr_offset;
	}

	double ret = transform_reduce(execution::seq,
		it1, it1 + common_size_wav,
		it2,
		0.0,
		std::plus<>{},
		[](double f, double s) {
			return std::abs(f - s);
		});

	return ret / common_size_wav / corr;
}

WavFile WavFile::intersection(const WavFile& other) const {

	if (fmt.sampleRate != other.fmt.sampleRate) {
		throw runtime_error("Unsuppotred intersection for dif sample rates"s);
	}

	if (fmt.audioFormat != other.fmt.audioFormat) {
		throw runtime_error("Unsuppotred intersection for dif audio format"s);
	}

	auto [corr_offset, corr] = CorrelationMax(other);//на сколько нужно сдвинуть this объект, что бы получить пересечение с other

	long long new_size_wav;
	long long start_index;
	if (corr_offset < 0) {
		//если отрицательный, то двигать начало надо у этого файла
		new_size_wav = min(static_cast<long long>(data[0].size()) + corr_offset, static_cast<long long>(other.data[0].size()));
		start_index = -corr_offset;
	} else {
		new_size_wav = min(static_cast<long long>(other.data[0].size()) - corr_offset, static_cast<long long>(data[0].size()));
		start_index = 0;
	}

	//cout << new_size_wav << " " << start_index << endl;

	WavFile ret = *this;

	for (auto& channel_data : ret.data) {
		channel_data = vector<double>(channel_data.begin() + start_index, channel_data.begin() + start_index + new_size_wav);
	}

	return ret;
}

pair<long long, double> WavFile::CorrelationMax(const WavFile& other) const {

	long long best_correlation = 0;
	double best_correlation_val;

	long long cnt = data[0].size();
	long long end = other.data[0].size() + 1;
	cnt = -cnt;
	best_correlation = cnt;
	best_correlation_val = correlation(other.data[0], data[0], cnt);

	for (++cnt; cnt < end; ++cnt) {
		double cur = std::abs(correlation(other.data[0], data[0], cnt));
		if (cur > best_correlation_val) {
			best_correlation_val = cur;
			best_correlation = cnt;
		} else if (cur < best_correlation_val / 10) {
			cnt += 1;//прорядим вычисления, correlation очень большая и толстая фнкция
		} else if (cur < best_correlation_val / 100 || cur < 0.1) {
			cnt += 100;//прорядим вычисления, correlation очень большая и толстая фнкция
		}
	}

	//cout << "best_correlation = " << best_correlation << " " << best_correlation_val << endl;
	return { best_correlation, best_correlation_val };
}

static double getVal(const vector<double>& v, long long x) {
	if (x < 0) {
		return 0;
	}
	if (x < v.size()) {
		return v[x];
	}
	return 0;
}

double WavFile::correlation(const vector<double>& f1, const vector<double>& f2, long long x) const {

	long long ret = 0;

	long long from = max(0ll, x);
	long long to = min(static_cast<long long>(f1.size()), static_cast<long long>(f2.size()) + x);
	if (to <= from) {
		return 0;
	}

	ret = transform_reduce(execution::par,
		f1.begin() + from, f1.begin() + to,
		f2.begin() + from - x,
		0ll,
		std::plus<>{},
		[](double f, double s) {
			return static_cast<long long>(f * s * 1000000000);
		});

	return static_cast<double>(ret) / 1000000000;
}

size_t WavFile::CalculateChankSize() const {
	size_t chankSize =
		4 //WAVE
		+ (subchanks.size() + 2) * 8 //subname + size
		;

	for (const auto& [key, val] : subchanks) {
		chankSize += val.size();
	}

	chankSize += 16;//fmt_size
	chankSize += data[0].size() * fmt.blockAlign;
	return chankSize;
}

void WavFile::PutWord(std::basic_fstream <unsigned char>& fs, uint16_t word) const {
	fs.put(word & 0xFF);
	fs.put((word >> 8) & 0xFF);
}

void WavFile::PutWord(std::basic_fstream <unsigned char>& fs, int16_t word) const {
	fs.put(word & 0xFF);
	fs.put((word >> 8) & 0xFF);
}

void WavFile::PutDWord(std::basic_fstream <unsigned char>& fs, uint32_t dword) const {
	fs.put(dword & 0xFF);
	fs.put((dword >> 8) & 0xFF);
	fs.put((dword >> 16) & 0xFF);
	fs.put((dword >> 24) & 0xFF);
}

void WavFile::PutDWord(std::basic_fstream <unsigned char>& fs, size_t dword) const {
	PutDWord(fs, static_cast<uint32_t>(dword));
}

void WavFile::save(const std::string& path) const {
	std::basic_fstream <unsigned char> file(path, ios_base::binary | ios_base::out | ios_base::trunc);
	if (!file) {
		throw runtime_error("can't create "s + path);
	}

	WriteChankHeader(file);

	WriteFmtSubchank(file);

	WriteOtherSubchanks(file);

	WriteDataSubchank(file);
}

void WavFile::WriteChankHeader(std::basic_fstream <unsigned char>& file) const {
	size_t chankSize = CalculateChankSize();

	file.write(reinterpret_cast<const unsigned char*>(chunkId.c_str()), chunkId.size());
	PutDWord(file, chankSize);
	file.write(reinterpret_cast<const unsigned char*>(format.c_str()), format.size());
}

void WavFile::WriteOtherSubchanks(std::basic_fstream <unsigned char>& file) const {
	for (const auto& [key, data] : subchanks) {
		file.write(reinterpret_cast<const unsigned char*>(key.c_str()), 4);
		PutDWord(file, data.size());
		file.write(data.data(), data.size());
	}
}

void WavFile::WriteFmtSubchank(std::basic_fstream <unsigned char>& file) const {
	file.write(reinterpret_cast<const unsigned char*>("fmt "), 4);
	PutDWord(file, 16u);

	PutWord(file, static_cast<uint16_t>(fmt.audioFormat));
	PutWord(file, fmt.numChannels);
	PutDWord(file, fmt.sampleRate);
	PutDWord(file, fmt.byteRate);
	PutWord(file, fmt.blockAlign);
	PutWord(file, fmt.bitsPerSample);
}

void WavFile::WriteDataSubchank(std::basic_fstream <unsigned char>& file) const {
	file.write(reinterpret_cast<const unsigned char*>("data"), 4);
	PutDWord(file, fmt.blockAlign * data[0].size());

	switch (fmt.bitsPerSample) {
	case 8:
		for (size_t i = 0; i < data[0].size(); ++i) {
			for (size_t channel = 0; channel < data.size(); ++channel) {
				int8_t val = numeric_limits<int8_t>::max() * data[channel][i];
				file.put(val);
			}
		}
		break;
	case 16:
		for (size_t i = 0; i < data[0].size(); ++i) {
			for (size_t channel = 0; channel < data.size(); ++channel) {
				double dval = static_cast<double>(numeric_limits<int16_t>::max() - 100) * data[channel][i];
				int16_t val = dval;
				PutWord(file, val);
			}
		}
		break;
	default:
		throw runtime_error("Not supported: do yourself");
	}
}

void WavFile::SetNumChannels(size_t numChannels) {
	if (numChannels == 0 || fmt.numChannels == numChannels) {
		return;
	}

	if (fmt.numChannels > numChannels) {
		int cnt = fmt.numChannels - numChannels + 1;
		size_t last = numChannels - 1;
		for (size_t i = 0; i < data[0].size(); ++i) {
			double n = transform_reduce(
				execution::seq,
				data.begin() + last, 
				data.end(), 0.0, 
				std::plus<>{},
				[i](const vector<double>& l) {
					return l[i];
				}
			);
			data[last][i] = n / cnt;
		}
		data.resize(numChannels);
	} else {
		vector<vector<double>> temp;
		data.reserve(numChannels);
		int cnt = numChannels - fmt.numChannels;
		temp.reserve(cnt);
		for (int i = 0; i < cnt; ++i) {
			temp.push_back(data[0]);
		}
		for (auto& v : temp) {
			data.push_back(move(v));
		}
	}
	fmt.blockAlign /= fmt.numChannels;
	fmt.byteRate /= fmt.numChannels;
	fmt.numChannels = numChannels;
	fmt.blockAlign *= fmt.numChannels;
	fmt.byteRate *= fmt.numChannels;
}

