/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#include "Util.h"
#include <algorithm>

namespace mtlib {

	std::string Util::space2underscore(std::string text)
	{
		for (std::string::iterator it = text.begin(); it != text.end(); ++it) {
			if (*it == ' ') {
				*it = '_';
			}
		}
		return text;
	}

	double Util::getCost(const std::vector<double>& values, const ObjectiveFunction& f)
	{
		switch (f) {
		case MAX: return getMax(values);
		case SUM: return getSum(values);
		case AVERAGE: return getAverage(values);
		case SQUARED_AVERAGE: return getSquaredAverage(values);
		case STD_DEVIATION: return getStdDeviation(values);
		case TEMPORAL_VARIANCE: return getTemporalVariance(values);
		default:return 0;
		}
	}

	double Util::getCostAll(const std::vector<double>& values, const std::vector<ObjectiveFunction>& funs)
	{
		double ret = 0;
		for (const auto& f : funs) {
			auto cost = Util::getCost(values, f);
			ret += cost;
		}
		return ret / funs.size();
	}

	std::string Util::getCurrentTimestamp()
	{
		
		time_t t = time(0);   // get time now
		struct tm * now = localtime(&t);

		std::stringstream ss;

		ss << (now->tm_year + 1900) << '-'
			<< (now->tm_mon + 1) << '-'
			<< now->tm_mday << '_'
			<< now->tm_hour << '-'
			<< now->tm_min << '-'
			<< now->tm_sec;

		return ss.str();

	}

	double Util::getMax(const std::vector<double>& values)
	{
		double ret = 0;
		for (const auto& v : values) {
			ret = std::max(ret, std::fabs(v));
		}
		return ret;
	}

	double Util::getSum(const std::vector<double>& values)
	{
		double ret = 0.0;
		for (const auto& v : values) {
			ret += std::fabs(v);
		}
		return ret;
	}

	double Util::getAverage(const std::vector<double>& values)
	{
		return getSum(values) / values.size();
	}

	double Util::getSquaredAverage(const std::vector<double>& values)
	{
		double ret = 0.0;
		for (const auto& v : values) {
			ret += (v * v);
		}
		return std::sqrt(ret / values.size());
	}

	double Util::getStdDeviation(const std::vector<double>& values)
	{
		double avg = 0;

		for (const auto& v : values) {
			avg += v;
		}
		
		avg /= values.size();

		double ret = 0.0;

		for (const auto& v : values) {
			ret += ((v - avg) * (v - avg));
		}
		
		return std::sqrt(ret / values.size());
	}
	double Util::getTemporalVariance(const std::vector<double>& values)
	{

		if (values.size() < 3) {
			return 0.0;
		}

		// interpolate between neighbouring values for the current value
		// and return the sum of differences
		double sum = 0;

		// right now just look at the two neighbours - left and right
		for (int i = 1; i < values.size() - 1; i++) {
			auto expected = (values[i - 1] + values[i + 1]) / 2;
			sum += pow((expected - values[i]), 2.0);
		}

		return sum / (values.size() - 2);

	}
}