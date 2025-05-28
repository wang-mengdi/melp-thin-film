//////////////////////////////////////////////////////////////////////////
// Auxillary functions relevent to Array(i.e., std::vector) operators
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __ArrayFunc_h__
#define __ArrayFunc_h__
#include "Common.h"
#include "AuxFunc.h"
#include <iostream>

namespace ArrayFunc {
	////Definitions
	//check if an Array contains invalid values, like nan, inf
	template<class T> bool Numerical_Check(const Array<T>& arr, const std::string& name = "", bool crash_on_fail = true);
	template<class T> int Largest_Norm_Element(const Array<T>& arr);//return -1 on empty
	template<class T> real Largest_Norm(const Array<T>& arr);

	//check deviation
	real Max_Abs_Error(const real x0, const Array<real>& arr);
	real Upper_Rel_Error(const real x0, const Array<real>& arr);//only count for arr[i]>x0. ignore all arr[i]<x0.
	real Max_Rel_Error(const real x0, const Array<real>& arr);
	real L2_Abs_Error(const real x0, const Array<real>& arr);
	real L2_Rel_Error(const real x0, const Array<real>& arr);

	//Sequential operations
	template<class T> void Resize_To(Array<T>& a, const int n, const T val);

	//Parallel Array Operations
	template<class T> void Array_Add(Array<T>& a, const Array<T>& b, const real c = 1);//a+=b*c
	template<class T> void Copy(Array<T>& a, const Array<T>& b);
	template<class Fvoid, class T> void Exec_Each(Fvoid f, Array<T>& a);
	template<class F1int> void Calc_Each(F1int f, Array<decltype(f(0))>& arr);
}

namespace ArrayFunc{
	////Implementations of template functions
	template<class T>
	bool Numerical_Check(const Array<T>& arr, const std::string& name, bool crash_on_fail)
	{
		int invalid_cnt = 0;
		for (int i = 0; i < arr.size(); i++) {
			//std::cout << "i: " << i << std::endl;
			//std::cout << "thing: " << arr[i] << std::endl;
			if (!Is_Valid_Number(arr[i])) {
				Info("ArrayFunc::Numerical_Check fails for {} at index {}: {}", name, i, arr[i]);
				invalid_cnt++;
				if (invalid_cnt >= 10) {
					Info(".....");
					break;
				}
			}
		}
		if (invalid_cnt) {
			if (crash_on_fail) exit(1);
			else return false;
		}
		return true;
	}
	template<class T>
	int Largest_Norm_Element(const Array<T>& arr)
	{
		int idx = -1; real max_norm = -1.0;
		for (int i = 0; i < arr.size(); i++) {
			real v = AuxFunc::Norm<T>(arr[i]);
			if (v >= max_norm) {
				idx = i;
				max_norm = v;
			}
		}
		return idx;
	}
	template<class T>
	real Largest_Norm(const Array<T>& arr)
	{
		int idx = Largest_Norm_Element(arr);
		if (idx < 0) return 0;
		else return AuxFunc::Norm<T>(arr[idx]);
	}
	template<class T>
	void Resize_To(Array<T>& a, const int n, const T val)
	{
		a.resize(n);
		std::fill(a.begin(), a.end(), val);
	}
	template<class T>
	void Array_Add(Array<T>& a, const Array<T>& b, const real c)
	{
		AuxFunc::Assert(a.size() == b.size(), "ArrayFunc::Array_Sum: size unmatch");
#pragma omp parallel for
		for (int i = 0; i < a.size(); i++) {
			a[i] += b[i] * c;
		}
	}
	template<class T>
	void Copy(Array<T>& a, const Array<T>& b)
	{
		AuxFunc::Assert(a.size() == b.size(), "ArrayFunc::Copy: size unmatch");
		int n = a.size();
#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			a[i] = b[i];
		}
	}
	template<class Fvoid, class T> 
	void Exec_Each(Fvoid f, Array<T>& a)
	{
		int N = a.size();
#pragma omp parallel for
		for (int i = 0; i < N; i++) f(i);
	}
	template<class F1int> 
	void Calc_Each(F1int f, Array<decltype(f(0))>& arr)
	{
		int N = arr.size();
		arr.resize(N);
#pragma omp parallel for
		for (int i = 0; i < N; i++) arr[i] = f(i);
	}
}

#endif
