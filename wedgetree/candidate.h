#ifndef HEADER_3DC82AAF72A4411C887EB04B07B09032
#define HEADER_3DC82AAF72A4411C887EB04B07B09032

#include <vector>
#include <stdexcept>
#include <iostream>
#include <limits>
#include <cmath>

struct Candidate {
	size_t timeseries_pos;
	size_t length;
	std::vector<double> const& timeseries;
	std::vector<double> series_normalized;
	double mean;
	double stddev;
	double max;
	double min;
	double range;
	Candidate() = delete;
	Candidate(const std::vector<double>& timeseries, size_t position, size_t length): 
		timeseries_pos(position),
		length(length),
		timeseries(timeseries),
		series_normalized(length, 0.0),
		mean(0.0),
		stddev(0.0),
		range(0.0),
		min(std::numeric_limits<double>::max()),
		max(std::numeric_limits<double>::min())
	{
		double ex = 0.0, ex2 = 0.0;
		size_t timeseries_len = timeseries.size();
		if (timeseries_pos + length > timeseries_len)
			throw new std::invalid_argument("Candidate length exceeds timeseries end position");
		for (size_t i = 0; i != length; ++i) {
			double d = timeseries[i + timeseries_pos];
			ex += d;
			ex2 += d*d;
		}
		mean = ex / length;
		stddev = ex2 / length;
		stddev = sqrt(stddev - mean * mean);
		for (size_t i = 0; i != length; i++) {
			double d = timeseries[i + timeseries_pos];
			series_normalized[i] = (d - mean) / stddev;
			if (d > max) max = d;
			if (d < min) min = d;
		}
		//lemire_envelope = LemireEnvelope(series_normalized, 0, WARPING_r);
		range = max - min;
	}
	Candidate(Candidate const& src) :
		timeseries_pos(src.timeseries_pos),
		length(src.length),
		timeseries(src.timeseries),
		series_normalized(src.series_normalized),
		mean(src.mean),
		stddev(src.stddev),
		max(src.max),
		min(src.min),
		range(src.range)
	{ 
		//std::cout << "copy constructor" << std::endl;
	}
	//visual studio does not yet support default move constructors as of VS 2013
	Candidate(Candidate&& src) :
		timeseries_pos(src.timeseries_pos),
		length(src.length),
		timeseries(src.timeseries),
		series_normalized(std::move(src.series_normalized)),
		mean(src.mean),
		stddev(src.stddev),
		max(src.max),
		min(src.min),
		range(src.range)
	{ }
};

#endif