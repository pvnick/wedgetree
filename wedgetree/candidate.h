#ifndef _CANDIDATE_H_
#define _CANDIDATE_H_

#include <vector>
#include <stdexcept>

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
};

#endif