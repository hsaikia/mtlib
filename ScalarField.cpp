/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#include "ScalarField.h"

namespace mtlib {

	ScalarField::ScalarField() {}

	ScalarField::~ScalarField() {}

	void ScalarField::init(const Dimensions& d) {
		m_dims.nx = d.nx;
		m_dims.ny = d.ny;
		m_dims.nz = d.nz;

		m_data.clear();
		m_data.resize(m_dims.getNumVoxels());

		minValue_ = std::numeric_limits<double>::max();
		maxValue_ = std::numeric_limits<double>::lowest();
	}

	void ScalarField::setDataAt(const size_t& idx, const double val) {
		m_data[idx] = val;
		minValue_ = std::min(minValue_, val);
		maxValue_ = std::max(maxValue_, val);
	}

	double ScalarField::F(const size_t& idx) const {
		return m_data[idx];
	}

	Dimensions ScalarField::getDims() const {
		return m_dims;
	}

	Vertex3D ScalarField::getVertex3(const size_t& idx, const Dimensions& dims) {
		return Vertex3D(idx % dims.nx, (idx / dims.nx) % dims.ny, (idx / (dims.nx * dims.ny)) % dims.nz);
	}

	size_t ScalarField::getVertex(const Vertex3D& v, const Dimensions& dims) {
		return dims.nx * (dims.ny * v.z + v.y) + v.x;
	}

	void ScalarField::getNeighbors(const size_t& idx, const Dimensions& dims, std::vector<size_t>& neighbors) {

		neighbors.clear();

		const int d[3] = { -1, 0, 1 };

		Vertex3D v = getVertex3(idx, dims);

		int X = (int)v.x;
		int Y = (int)v.y;
		int Z = (int)v.z;

		for (size_t i = 0; i < 3; ++i)
		{
			if ((X + d[i] < 0) || (X + d[i] > dims.nx - 1)) continue;
			for (size_t j = 0; j < 3; ++j)
			{
				if ((Y + d[j] < 0) || (Y + d[j] > dims.ny - 1)) continue;
				for (size_t k = 0; k < 3; ++k)
				{
					if ((Z + d[k] < 0) || (Z + d[k] > dims.nz - 1)) continue;
					if (d[i] == 0 && d[j] == 0 && d[k] == 0) continue;
					//printf("[x y z] = [%d %d %d]\n", dx[i], dy[j], dz[k]);

					Vertex3D nv(v.x + d[i], v.y + d[j], v.z + d[k]);
					neighbors.push_back(getVertex(nv, dims));

				}
			}
		}
	}

	double ScalarField::normEuclideanDist2(const size_t & idx1, const size_t & idx2, const Dimensions& dims)
	{
		Vertex3D v1 = getVertex3(idx1, dims);
		Vertex3D v2 = getVertex3(idx2, dims);

		Vertex3Df v1f(static_cast<double>(v1.x) / dims.nx, static_cast<double>(v1.y) / dims.ny, static_cast<double>(v1.z) / dims.nz);
		Vertex3Df v2f(static_cast<double>(v2.x) / dims.nx, static_cast<double>(v2.y) / dims.ny, static_cast<double>(v2.z) / dims.nz);

		return (v1f.x - v2f.x)*(v1f.x - v2f.x) + (v1f.y - v2f.y)*(v1f.y - v2f.y) + (v1f.z - v2f.z)*(v1f.z - v2f.z);
	}

	void ScalarField::sortMap(std::function<bool(Voxel, Voxel)> comp, std::vector<size_t>& vertexOrder) const {

		size_t numVox = m_dims.getNumVoxels();

		vertexOrder.clear();
		vertexOrder.resize(numVox);

		std::vector<Voxel> Map(numVox);

		for (size_t i = 0; i < numVox; i++) {
			Map[i].idx = i;
			Map[i].val = F(i);
		}

		sort(Map.begin(), Map.end(), comp);

		for (size_t i = 0; i < Map.size(); i++) {
			vertexOrder[i] = Map[i].idx;
		}

	}

	double ScalarField::getValueRange() const
	{
		return maxValue_ - minValue_;
	}

	double ScalarField::getMinVal() const
	{
		return minValue_;
	}

	double ScalarField::getMaxVal() const
	{
		return maxValue_;
	}

	double ScalarField::getHeightDifference(const size_t & idx1, const size_t & idx2) const
	{
		return std::fabs(F(idx1) - F(idx2));
	}

	double ScalarField::getNormalizedHeightDifference(const size_t & idx1, const size_t & idx2) const
	{
		return std::fabs(F(idx1) - F(idx2)) / (maxValue_ - minValue_);
	}

}