//////////////////////////////////////////////////////////////////////////
// ArrayIO, here provides more convenient I/O, based on basic functions provided by File.h
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __ArrayIO_h__
#define __ArrayIO_h__


#include "File.h"
#include "Common.h"

namespace BinaryDataIO {
	//// Declaration of functions

	////Write Vector array to binary, padding to 3D. Memory layout:
	//Bytes [0,4): an integer n, indicating number of vectors (NOT number of bytes)
	//Bytes [4,4+n*3*sizeof(T)): binary value of n 3D vectors.
	//Note: always store 3D vectors in binary file. If the original vectors are less than 3D, we fill it with 0.
	template<class T, int d> bool Write_Vector_Array_3D(const std::string& file_name, const Array<Vector<T, d> >& arr);
	template<class T, int d, class F> bool Write_Vector_Array_3D(const std::string& file_name, F& f, std::uint32_t n);//i-th vector is f(i)
	template<class T, int d> bool Read_Vector_Array_3D(const std::string& file_name, Array<Vector<T, d> >& arr);

	////No padding, just write memory.
	//Bytes [0,4): an integer n, indicating number of scalars (NOT number of bytes)
	//Bytes [4,4+n*sizeof(T)): binary value of n elements
	//Note that for Eigen::Vector or Eigen::Matrix, they're just the flattened element array.

	//Stream Style
	void Write_Array_Stream_Content(std::ostream& output, const Array<bool>& arr);
	template<class T> void Write_Array_Stream_Content(std::ostream& output, const Array<T>& arr) {
		const T* data = (const T*)(arr.data());
		File::Write_Binary_Array<T>(output, data, arr.size());
	}
	template<class T> void Write_Array_Stream(std::ostream& output, const Array<T>& arr) {
		std::uint32_t n = (std::uint32_t)arr.size();
		File::Write_Binary<std::uint32_t>(output, n);
		Write_Array_Stream_Content(output, arr);
	}
	template<class T, class F> void Write_Array_Stream(std::ostream& output, F& f, std::uint32_t n) {
		Array<T> arr(n);
		for (int i = 0; i < n; i++) arr[i] = f(i);
		Write_Array_Stream(output, arr);
	}
	void Read_Array_Stream_Content(std::istream& input, Array<bool>& arr, const std::uint32_t& n);
	template<class T> void Read_Array_Stream_Content(std::istream& input, Array<T>& arr, const std::uint32_t& n) {
		arr.resize(n);
		T* data = (T*)(arr.data());
		File::Read_Binary_Array<T>(input, data, n);
	}
	template<class T> void Read_Array_Stream(std::istream& input, Array<T>& arr) {
		std::uint32_t n;
		File::Read_Binary<std::uint32_t>(input, n);
		Read_Array_Stream_Content(input, arr, n);
	}

	//File Style
	template<class T> bool Write_Array(const std::string& file_name, const Array<T>& arr) {
		std::ofstream output(file_name, std::ios::binary); 
		if (!output) return false; 
		Write_Array_Stream(output, arr); 
		output.close(); 
		return true;
	}
	template<class T, class F> bool Write_Array(const std::string& file_name, F& f, std::uint32_t n) {
		Array<T> arr(n);
		for (int i = 0; i < n; i++) arr[i] = f(i);
		return Write_Array(file_name, arr);
	}
	template<class T> bool Read_Array(const std::string& file_name, Array<T>& arr) {
		std::ifstream input(file_name, std::ios::binary); 
		if (!input) return false; 
		Read_Array_Stream(input, arr); 
		input.close(); 
		return true;
	}


	//// Implementation of functions


	template<class T, int d>
	bool Write_Vector_Array_3D(const std::string& file_name, const Array<Vector<T, d>>& arr)
	{
		assert(1 <= d && d <= 3);
		std::ofstream output(file_name, std::ios::binary);
		if (!output) return false;
		std::uint32_t n = (std::uint32_t)arr.size();
		File::Write_Binary(output, n);
		T* data = new T[n * 3];
		memset(data, 0, n * 3 * sizeof(T));
#pragma omp parallel for
		for (int i = 0; i < (int)n; i++) {
			for (int axis = 0; axis < d; axis++) data[i * 3 + axis] = arr[i](axis);
		}
		File::Write_Binary_Array(output, data, n * 3);
		delete[] data;
		output.close();
		return true;
	}
	template<class T, int d, class F>
	bool Write_Vector_Array_3D(const std::string& file_name, F& f, std::uint32_t n)
	{
		assert(1 <= d && d <= 3);
		std::ofstream output(file_name, std::ios::binary);
		if (!output) return false;
		File::Write_Binary(output, n);
		T* data = new T[n * 3];
		memset(data, 0, n * 3 * sizeof(T));
#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			Vector<T, d> vec = f(i);
			for (int axis = 0; axis < d; axis++) data[i * 3 + axis] = vec(axis);
		}
		File::Write_Binary_Array(output, data, n * 3);
		delete[] data;
		output.close();
		return true;
	}
	template<class T, int d>
	bool Read_Vector_Array_3D(const std::string& file_name, Array<Vector<T, d>>& arr)
	{
		assert(1 <= d && d <= 3);
		std::ifstream input(file_name, std::ios::binary);
		if (!input) return false;
		std::uint32_t n;
		File::Read_Binary(input, n);
		arr.resize(n);
		T* data = new T[n * 3];
		File::Read_Binary_Array(input, data, n * 3);
#pragma omp parallel for
		for (int i = 0; i < (int)n; i++) {
			for (int axis = 0; axis < d; axis++) arr[i](axis) = data[i * 3 + axis];
		}
		delete[] data;
		input.close();
		return true;
	}
}

#endif
