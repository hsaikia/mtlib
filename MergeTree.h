/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#pragma once

#include "globals.h"
#include "ScalarField.h"
#include <map>
#include <set>
#include <functional>
#include <fstream>
#include <istream>
#include <string>

namespace mtlib {

	typedef size_t VertexIdx;
	typedef int FeatureIdx;

	enum CPType { EXTREMUM, SADDLE, GBL_EXT, REGULAR };
	enum SimplificationMethod { SIM_THRESHOLD, SIM_NUMBER };
	enum HierarchyType { PERSISTENCE, INTEGRATED_INTENSITY_DIFFERENCE };

	struct CriticalPoint {
		VertexIdx idx;
		CPType typ;
	};

	struct Branch {
		size_t brReprExtIdx;	// branch representative extremum index
		size_t saddleIdx;		// saddle index at which the branch is plucked off from
		double weight;			// normalized height difference of this branch with root branch being 1.0 in case persistence simplication is used
								// sum of f(x) - f(s) where x is a voxel in this branch, and s is the saddle in case integrated intensity simplification is used.
		size_t level;
		int vol;				// cumulative volume of this branch (branch volume + all child branches' volume)
	};

	struct ConciseSMT {
		std::map<int, int> unAugIdxTree;						// map {cp_i => i --> cp_j => j}, -1 is global extremum
		std::vector<CriticalPoint> cps;							// remaining subtrees
	};

	struct SimplifiedMergeTree {
		std::map<FeatureIdx, FeatureIdx> unAugIdxTree;						// map {cp_i => i --> cp_j => j}, -1 is global extremum
		std::vector<CriticalPoint> cps;										// remaining critical points
		std::vector<FeatureIdx> levelSetReprIdxMap;							// map {v_i -> cp_j => j}
		std::vector<double> simplifiedScalarField;							// All voxels of simplified edges get flattened out
		std::map <FeatureIdx, std::vector<FeatureIdx> > iUnAugIdxTree;		// map {cp_i => i --> child(cp_i) == <cp_j> => <j> }
		FeatureIdx rootIdxTree;												// root of inverted tree. Last saddle.
		Dimensions dims;

		SimplifiedMergeTree() {

		}

		~SimplifiedMergeTree() {

		}

		void getIsoLevelSets(std::vector<FeatureIdx>& belongsTo, FeatureIdx& maxFeatureIdx, const double isoValue);

		ConciseSMT getConciseSMT() {
			ConciseSMT c;
			c.cps = cps;
			c.unAugIdxTree = unAugIdxTree;
			return c;
		}
		void writeUnAugTree(const std::string& filename) {
			std::ofstream f;
			f.open(filename);

			for (const auto& e : unAugIdxTree) {
				f << cps[e.first].idx << "->" << cps[e.second].idx << "\n";
			}

			f.close();
		}
	};

	class MergeTree {

	public:
		struct doCompare
		{
			doCompare(const MergeTree& info) : m_info(info) { } // only if you really need the object state
			const MergeTree& m_info;

			bool operator()(const size_t& idx1, const size_t& idx2)
			{
				size_t rep_idx1 = idx1;
				size_t rep_idx2 = idx2;
				if (m_info.brIdxMap_.find(idx1) != m_info.brIdxMap_.end()) {
					rep_idx1 = m_info.brs_[m_info.brIdxMap_.find(idx1)->second].brReprExtIdx;
				}
				if (m_info.brIdxMap_.find(idx2) != m_info.brIdxMap_.end()) {
					rep_idx2 = m_info.brs_[m_info.brIdxMap_.find(idx2)->second].brReprExtIdx;
				}
				return m_info.sortedOrder_[rep_idx1] < m_info.sortedOrder_[rep_idx2];

			}
		};
		MergeTree(ScalarField& sf, std::vector<VertexIdx>& vertexOrder, const HierarchyType& htyp);

		void createTree();

		void createAugTree();

		static void drawUnAugTree(std::ofstream& file, const std::map<int, int>& unAugTree, const std::string& prefix);

		void drawUnAugTree(std::ofstream& file, const std::string& prefix);

		void writeUnAugTree(const std::string& filename);

		void writeSimpMap(const std::string& filename);

		void drawBDTree(std::ofstream& file, const std::string& prefix);

		void writeStatistics(std::string filename);

		void createSimplifiedMergeTree(double threshold, size_t numToRemain, SimplifiedMergeTree& smt, SimplificationMethod method) const;

		static bool isAncestor(const int & nodeIdx, const int & ancestorIdx, const std::map<int, int>& unAugTree);

		size_t getTreeSize(); // num of critical points

		size_t getBDTSize(); // number of branches

		std::vector<size_t> augTree_;								// size = entire dataset

		std::vector<CriticalPoint> cps_; // making public for testing, should be private

	private:
		ScalarField& sf_;
		HierarchyType hierarchyType_;
		std::vector<VertexIdx>& vertexOrder_;	//ordered list of vertices
		std::vector<VertexIdx> sortedOrder_;	// vertex to position in sorted array
		size_t numVoxels_;
		//Union Find functions

		VertexIdx rootCpIdx_;											// index of global min (max) of join(split) tree

		std::vector<VertexIdx> activeLevelSetReprMap_;					// updates using the find function to the bottom-most CP

		std::map<VertexIdx, VertexIdx> unAugTree_;						// map to indicate child CP to parent CP idx in volume, size = size of CPs
		//std::vector<CriticalPoint> cps_;
		std::map<VertexIdx, FeatureIdx> cpIdxMap_;						// map of cp's volume idx to index in 'cps' array

		//Branch Decomposition
		std::map<size_t, size_t> BDTree_;								// map to indicate child CP to parent CP in 'brs' array, size = size of BDT
		std::vector<Branch> brs_;
		std::map<VertexIdx, FeatureIdx> brIdxMap_;						// map of cp's volume idx to index in 'brs' array

		//simplification map, based on the BDT
		std::map<FeatureIdx, FeatureIdx> simpMap_;						// x simplifies to simpMap_[x], so any member of x points to simpMap_[x] after x is simplified

		std::vector<VertexIdx> levelSetReprMap_;						// map to say which voxel belongs to which CP, size = entire dataset

		size_t findComponent(const size_t& idx);
		void createNewComponent(const size_t& idx, const CPType& typ);
		void mergeWithOneComponent(const size_t& idx, const size_t& belongsTo);
		void mergeWithManyComponents(const size_t& idx, const std::vector<size_t>& saddleOf);
		//void createSimplificationOrder();
		void doConsistencyCheck();
		int createSimplificationOrderForPath(const int& startCpIdx, const std::map<size_t, bool>& survivedCPs, std::map<int, int>& simpMap) const;
		void setSimpMap(int idx, int master);
	};
}