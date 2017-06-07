/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#pragma once

#include "globals.h"
#include <vector>
#include <functional>
#include <algorithm>

namespace mtlib {

	class ScalarField {
	public:
		ScalarField();

		~ScalarField();

		void init(const Dimensions& d);

		void setDataAt(const size_t& idx, const double val);

		std::vector<double> getData() const { return m_data; }

		//const

		double F(const size_t& idx) const;

		void sortMap(std::function<bool(Voxel, Voxel)> comp, std::vector<size_t>& vertexOrder) const;

		double getValueRange() const;

		double getMinVal() const;

		double getMaxVal() const;

		double getHeightDifference(const size_t& idx1, const size_t& idx2) const;

		double getNormalizedHeightDifference(const size_t& idx1, const size_t& idx2) const;

		Dimensions getDims() const;

		//statics

		static size_t getVertex(const Vertex3D& v, const Dimensions& dims);

		static Vertex3D getVertex3(const size_t& idx, const Dimensions& dims);

		static void getNeighbors(const size_t& idx, const Dimensions& dims, std::vector<size_t>& neighbors);

		static double normEuclideanDist2(const size_t& idx1, const size_t& idx2, const Dimensions& dims); // returns the squared Euclidean Distance between the two volume indices

	private:
		double minValue_;
		double maxValue_;
		std::vector<double> m_data;
		Dimensions m_dims;
	};

}