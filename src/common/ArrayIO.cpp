#include "ArrayIO.h"

void BinaryDataIO::Write_Array_Stream_Content(std::ostream& output, const Array<bool>& arr)
{
	std::uint32_t n = (std::uint32_t)arr.size();
	bool* data = new bool[n];
#pragma omp parallel for
	for (int i = 0; i < (int)n; i++) data[i] = arr[i];
	File::Write_Binary_Array<bool>(output, data, n);
	delete[] data;
}

void BinaryDataIO::Read_Array_Stream_Content(std::istream& input, Array<bool>& arr, const std::uint32_t& n) {
	arr.resize(n);
	bool* data = new bool[n];
	File::Read_Binary_Array<bool>(input, data, (int)n);
#pragma omp parallel for
	for (int i = 0; i < (int)n; i++) arr[i] = data[i];
	delete[] data;
}