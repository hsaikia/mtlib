/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#include "FileIO.h"

namespace mtlib {

	void FileIO::readRegionsFromCSV(const std::string & filename, std::map<dagkey, Region>& dagNodes)
	{
		dagNodes.clear();
		std::ifstream f;
		f.open(filename);

		std::string line;

		bool first = true;

		int count = 0;

		while (std::getline(f, line))
		{
			//std::cout << "Reading Line " << count << "\n";

			std::istringstream iss(line);

			if (first) {
				first = false;

				std::string token;

				while (std::getline(iss, token, ',')) {
					count++;
				}

				std::cout << "Histogram Size is " << count - 11 << "\n";

				continue;
			}

			Region R;

			R.hist.init(count - 11);

			int sttype, par_id;
			dagkey key;

			char c; // comma

			iss >> R.timestep >> c >> R.index >> c >> R.vindex >> c >> sttype >> c >> par_id >> c >> key >> c >> R.hist.hdf
				>> c >> R.hist.vol >> c >> R.hist.cm.x >> c >> R.hist.cm.y >> c >> R.hist.cm.z;

			for (auto& bin : R.hist.bins) {
				iss >> c >> bin;
			}

			dagNodes[key] = R;

		}

		f.close();
	}


	void FileIO::readMergeTreesFromCSV(const std::string & filename, std::vector<std::map<int, int>>& unAugTrees)
	{
		unAugTrees.clear();

		std::ifstream f;
		f.open(filename);

		std::string line;

		bool first = true;

		int count = 0;

		while (std::getline(f, line))
		{
			if (first) {
				first = false;
				continue;
			}

			count++;

			std::istringstream iss(line);
			int t, id, vid, sttype, par_id, seq_id, st_vol;
			float st_hdf;
			Vertex3Df cm;
			char c; // comma
			iss >> t >> c >> id >> c >> vid >> c >> sttype >> c >> par_id >> c >> seq_id >> c >> st_hdf >> c >> st_vol >> c >> cm.x >> c >> cm.y >> c >> cm.z;

			if (t >= unAugTrees.size()) {
				unAugTrees.resize(t + 1);
			}
			unAugTrees[t].insert(std::pair<int, int>(id, par_id));

		}

		std::cout << "[Tracking::readMergeTreesFromCSV] Number of Merge Trees " << unAugTrees.size() << ". Num lines read " << count << "\n";

		f.close();
	}

	void FileIO::writeCostMatrixToDisk(const std::string & filename, const std::vector<std::vector<double>>& costMatrix)
	{
		std::ofstream f;
		f.open(filename);

		f << costMatrix.size() << "\n";

		for (size_t i = 0; i < costMatrix.size(); i++) {
			for (size_t j = 0; j < costMatrix[i].size(); j++) {
				f << costMatrix[i][j] << "\t";
			}
			f << "\n";
		}
		f.close();
	}

	void FileIO::readCostMatrixFromDisk(const std::string & filename, std::vector<std::vector<double>>& costMatrix)
	{
		std::ifstream f;
		f.open(filename);

		size_t n;
		f >> n;

		costMatrix.clear();
		costMatrix.resize(n, std::vector<double>(n));

		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n; j++) {
				f >> costMatrix[i][j];
			}
		}

		f.close();
	}

	void FileIO::drawMatches(const std::map<size_t, size_t>& matches, const std::string& prefix1, const std::string& prefix2, std::ofstream & file)
	{
		for (const auto& match : matches) {
			file << prefix1 << match.first << "->" << prefix2 << match.second << "[color=\"#00b300\",style=\"setlinewidth(2)\"];\n";
		}
	}

	void FileIO::addHeaderToStatisticsFile(std::ofstream & fs, const int hist_buckets)
	{
		fs << "Timestep,Merge_Tree_Id,Volume_Vertex,Subtree_Repr_Type,Merge_Tree_Parent_Id,DAG_Key,Subtree_Height_Difference,Subtree_Volume,Subtree_Center_of_Mass_x,Subtree_Center_of_Mass_y,Subtree_Center_of_Mass_z";
		for (auto i = 0; i < hist_buckets; i++) {
			fs << ",Hist_" << i;
		}
		fs << "\n";
	}

	void FileIO::addToStatistics(std::ofstream & fs, const int & timestep, const ConciseSMT & smt, const std::vector<Hist>& hists, const double globalRange)
	{
		for (auto i = 0; i < smt.cps.size() - 1; i++) {

			int mtParentId = i;
			if (smt.unAugIdxTree.find(i) != smt.unAugIdxTree.end()) {
				mtParentId = smt.unAugIdxTree.find(i)->second;
			}
			else {
				continue;
			}

			if (mtParentId == smt.unAugIdxTree.size()) {
				mtParentId = -1; // global minimum
			}

			auto key = DagNode::makeKey(timestep, i);

			fs << timestep << ","
				<< i << ","
				<< smt.cps[i].idx << ","
				<< smt.cps[i].typ << ","
				<< mtParentId << ","
				<< key << ","
				<< (hists[i].maxVal - hists[i].minVal) * globalRange << ","
				<< hists[i].vol << ","
				<< hists[i].cm.x << ","
				<< hists[i].cm.y << ","
				<< hists[i].cm.z;

			for (const auto& bin : hists[i].bins) {
				fs << "," << bin;
			}

			fs << "\n";
		}
	}

}
