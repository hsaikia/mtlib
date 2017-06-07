/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#pragma once
#include <vector>
#include <sstream>

namespace mtlib {

	enum ObjectiveFunction { MAX, SUM, AVERAGE, SQUARED_AVERAGE, STD_DEVIATION, TEMPORAL_VARIANCE, NUM_OBJECTIVE_FUNCTIONS };

	class Util {
	public:
		// The following two functions are taken from http://www.geeksforgeeks.org/counting-inversions/
		// and modified by Himangshu Saikia
		template <typename T>
		static int mergeSort(std::vector<T>& arr);
		template <typename T>
		static int _mergeSort(std::vector<T>& arr, std::vector<T>& temp, int left, int right);
		template <typename T>
		static int merge(std::vector<T>& arr, std::vector<T>& temp, int left, int mid, int right);
		static std::string space2underscore(std::string text);
		
		static double getCost(const std::vector<double>& values, const ObjectiveFunction& f);

		static double getCostAll(const std::vector<double>& values, const std::vector<ObjectiveFunction>& funs);

		static std::string getCurrentTimestamp();

		static double getMax(const std::vector<double>& values);
		static double getSum(const std::vector<double>& values);
		static double getAverage(const std::vector<double>& values);
		static double getSquaredAverage(const std::vector<double>& values);
		static double getStdDeviation(const std::vector<double>& values);
		static double getTemporalVariance(const std::vector<double>& values);

	};

	template<typename T>
	inline int Util::mergeSort(std::vector<T>& arr)
	{
		std::vector<T> temp(arr.size());
		return _mergeSort(arr, temp, 0, arr.size() - 1);
	}

	template<typename T>
	inline int Util::_mergeSort(std::vector<T>& arr, std::vector<T>& temp, int left, int right)
	{
		int mid, inv_count = 0;
		if (right > left)
		{
			/* Divide the array into two parts and call _mergeSortAndCountInv()
			for each of the parts */
			mid = (right + left) / 2;

			/* Inversion count will be sum of inversions in left-part, right-part
			and number of inversions in merging */
			inv_count += _mergeSort(arr, temp, left, mid);
			inv_count += _mergeSort(arr, temp, mid + 1, right);

			/*Merge the two parts*/
			inv_count += merge(arr, temp, left, mid + 1, right);
		}
		return inv_count;
	}

	template<typename T>
	inline int Util::merge(std::vector<T>& arr, std::vector<T>& temp, int left, int mid, int right)
	{
		int i, j, k;
		int inv_count = 0;

		i = left; /* i is index for left subarray*/
		j = mid; /* j is index for right subarray*/
		k = left; /* k is index for resultant merged subarray*/
		while ((i <= mid - 1) && (j <= right))
		{
			if (arr[i] <= arr[j])
			{
				temp[k++] = arr[i++];
			}
			else
			{
				temp[k++] = arr[j++];

				/*this is tricky -- see above explanation/diagram for merge()*/
				inv_count = inv_count + (mid - i);
			}
		}

		/* Copy the remaining elements of left subarray
		(if there are any) to temp*/
		while (i <= mid - 1)
			temp[k++] = arr[i++];

		/* Copy the remaining elements of right subarray
		(if there are any) to temp*/
		while (j <= right)
			temp[k++] = arr[j++];

		/*Copy back the merged elements to original array*/
		for (i = left; i <= right; i++)
			arr[i] = temp[i];

		return inv_count;

	}

}