/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/


#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <string>

namespace mtlib {

	struct Vertex3D {
		size_t x;
		size_t y;
		size_t z;

		Vertex3D() {}

		Vertex3D(size_t x_, size_t y_, size_t z_) {
			x = x_;
			y = y_;
			z = z_;
		}

	};

	struct Vertex3Df {
		double x;
		double y;
		double z;

		Vertex3Df() {
			x = 0; y = 0; z = 0;
		}

		Vertex3Df(double x_, double y_, double z_) {
			x = x_;
			y = y_;
			z = z_;
		}

		Vertex3Df operator*(const double& s) {
			Vertex3Df ret;
			ret.x = x * s;
			ret.y = y * s;
			ret.z = z * s;
			return ret;
		}

		Vertex3Df operator/(const double& s) {
			Vertex3Df ret;
			ret.x = x / s;
			ret.y = y / s;
			ret.z = z / s;
			return ret;
		}

		Vertex3Df operator+(const Vertex3Df& v1) {
			Vertex3Df ret;
			ret.x = x + v1.x;
			ret.y = y + v1.y;
			ret.z = z + v1.z;
			return ret;
		}

		/*void operator=(const Vertex3Df& v1) {
			x = v1.x;
			y = v1.y;
			z = v1.z;
		}*/

		static double EucDist(const Vertex3Df& v1, const Vertex3Df& v2) {
			return sqrt((v1.x - v2.x)*(v1.x - v2.x) + (v1.y - v2.y)*(v1.y - v2.y) + (v1.z - v2.z)*(v1.z - v2.z));
		}

	};

	struct Edge {
		int n1;
		int n2;

		Edge() {

		}

		Edge(int n1_, int n2_) {
			n1 = n1_;
			n2 = n2_;
		}

	};

	typedef struct Dimensions {
		size_t nx;
		size_t ny;
		size_t nz;

		Dimensions() {

		}

		Dimensions(size_t x_, size_t y_, size_t z_) {
			nx = x_;
			ny = y_;
			nz = z_;
		}

		void init(size_t x_, size_t y_, size_t z_) {
			nx = x_;
			ny = y_;
			nz = z_;
		}

		size_t getNumVoxels() const {
			return nx * ny * nz;
		}

	} Dimensions;

	typedef struct Voxel {
		double val;
		size_t idx;
		size_t sortedIdx;
	} Voxel;

	//Histograms

	struct Hist {
		std::vector<int> bins;
		int vol;
		double totmass;
		Vertex3Df cm;
		double minVal;
		double maxVal;
		double hdf;
		int buckets;

		long long getMemory() const {
			return bins.size() * sizeof(int) + 2 * sizeof(int) + 4 * sizeof(double) + sizeof(Vertex3Df);
		}

		Hist() {
			buckets = 1;
			bins.resize(1, 0);
			vol = 0;
			minVal = std::numeric_limits<double>::max();
			maxVal = std::numeric_limits<double>::lowest();
			totmass = 0;
			cm.x = 0;
			cm.y = 0;
			cm.z = 0;
		}

		void init(int buckets_) {
			buckets = buckets_;
			bins.resize(buckets, 0);
		}

		void add_elem(const double& minValue, const double& maxValue, const double& x, Vertex3Df& v) {

			minVal = std::min(minVal, x);
			maxVal = std::max(maxVal, x);

			double normVal = (x - minValue) / (maxValue - minValue);

			cm = (cm * totmass + v * normVal) / (totmass + normVal);
			totmass += normVal;

			int idx = static_cast<int>(buckets * normVal);
			if (idx == buckets) {
				idx--;
			}
			bins[idx]++;
			vol++;
		}

		void add_hist(Hist& h) {
			minVal = std::min(minVal, h.minVal);
			maxVal = std::max(maxVal, h.maxVal);
			for (size_t i = 0; i < buckets; i++) {
				bins[i] += h.bins[i];
			}

			cm = (cm * totmass + h.cm * h.totmass) / (totmass + h.totmass);
			totmass += h.totmass;
			vol += h.vol;

		}

	};

	//Tracking Regions

	struct Region {
		int timestep;
		int index; // from i_0 - i_n
		int vindex; // index of i_0 in volume - index of i_n in volume
		Hist hist;

		Region() {

		}

		Region(int t, int i) {
			timestep = t;
			index = i;
		}

		static bool compareRegionTopo(const Region& R1, const Region& R2) {
			return R1.timestep < R2.timestep;
		}

		long long getMemory() const {
			return 3 * sizeof(int) + hist.getMemory();
		}

	};
}