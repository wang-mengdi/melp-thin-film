//////////////////////////////////////////////////////////////////////////
// Particles
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Particles.h"
#include "File.h"
#include "AuxFunc.h"

//////////////////////////////////////////////////////////////////////////
////points

template<int d, typename T> void Points<d, T>::Print_Attributes(void) {
	std::cout << "Total " << att_map.size() << " attributes:\n";
	for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
		std::cout << (iter->first) << " ";
	}
	std::cout << "\n";
}

template<int d, typename T> void Points<d, T>::Resize(const int size) {
	for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
		(iter->second)->Resize(size);
	}
}

template<int d, typename T> void Points<d, T>::Reserve(const int size) {
	for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
		(iter->second)->Reserve(size);
	}
}

template<int d, typename T> int Points<d, T>::Add_Element() {
	int idx = -1;
	for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
		int now_idx = (iter->second)->Add_Element();
		if (idx != -1 && idx != now_idx) {
			std::cerr << "Points::Add_Element error: unexpected size of " << iter->first << "\n";
			exit(1);
		}
		idx = now_idx;
	}
	return idx;
}

template<int d, typename T> int Points<d, T>::Add_Position_Only(const VectorD& pos) {
	int idx = Add_Element();
	X(idx) = pos;
	return idx;
}

template<int d, typename T> int Points<d, T>::Add_Duplicate_Element(int src_idx) {
	int idx = -1;
	for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
		int new_idx = (iter->second)->Add_Elements(1);
		if (idx != -1 && idx != new_idx) {
			std::cerr << "Points::Add_Duplicate_Element error: unexpected size of " << iter->first << "\n";
			exit(1);
		}
		idx = new_idx;
		iter->second->Copy_Element(new_idx, iter->second, src_idx);
	}
	return idx;
}

template<int d, typename T> int Points<d, T>::Add_Elements(const int n) {
	int idx = -1;
	for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
		int now_idx = (iter->second)->Add_Elements(n);
		if (idx != -1 && idx != now_idx) {
			std::cerr << "Points::Add_Element error: unexpected size of " << iter->first << "\n";
			exit(0);
		}
		idx = now_idx;
	}
	return idx;
}

template<int d, typename T> int Points<d, T>::Join(const Points<d, T>& src) {
	int n_old = Size(), n_append = src.Size();
	Resize(n_old + n_append);
	for (int i = 0; i < n_append; i++) {
		Copy_Element_From(n_old + i, src, i);
	}
	return n_old;
}

template<int d, typename T> int Points<d, T>::Delete_Elements(const Array<int>& is_deleted) {
	int deleted_size = 0;
	for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
		deleted_size = (iter->second)->Delete_Elements(is_deleted);
	}
	return deleted_size;
}

template<int d, typename T> int Points<d, T>::Delete_Elements(const std::set<int>& to_delete)
{
	int deleted_size = 0;
	for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
		deleted_size = (iter->second)->Delete_Elements(to_delete);
	}
	return deleted_size;
}

template<int d, typename T> int Points<d, T>::Delete_Elements_Safe(Array<int> is_deleted) {
	return Delete_Elements(is_deleted);
}

//we call filter_func for all n elements once at the beginning, so you can safely use the class itself in filter_func.
//return number of remaining elements
template<int d, typename T> int Points<d, T>::Filter_Elements(std::function<bool(const int)> filter_func)
{
	int n = Size();
	Array<int> is_deleted(n);
	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		is_deleted[i] = !filter_func(i);
	}
	return Delete_Elements(is_deleted);
}

template<int d, typename T> void Points<d, T>::Copy_Element_From(const int idx, const Points<d, T>& src, const int src_idx) {
	for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
		auto src_iter = src.att_map.find(iter->first); assert(src_iter != src.att_map.end());
		iter->second->Copy_Element(idx, src_iter->second, src_idx);
	}
}

template<int d, typename T>
void Points<d, T>::Copy_Element_From(const int idx, const Points<d, T>& src, const int src_idx, const Array<std::string> attrs)
{
	for (int i = 0; i < attrs.size(); i++) {
		const std::string& att_name = attrs[i];
		auto dst_iter = att_map.find(att_name); assert(dst_iter != att_map.end());
		auto src_iter = src.att_map.find(att_name); assert(src_iter != src.att_map.end());
		dst_iter->second->Copy_Element(idx, src_iter->second, src_idx);
	}
}

template<int d, typename T> void Points<d, T>::Write_Binary(std::ostream& output) const
{
	File::Write_Binary(output, Size());
	for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
		(iter->second)->Write_Binary_Content(output);
	}
}

template<int d, typename T> void Points<d, T>::Read_Binary(std::istream& input)
{
	int n = 0; File::Read_Binary(input, n);
	for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
		(iter->second)->Read_Binary_Content(input, n);
	}
}

template<int d, typename T>
bool Points<d, T>::Save_Snapshot(const std::string& file_name)const
{
	std::ofstream output(file_name, std::ios::binary);
	if (!output) {
		std::cerr << "Points::Save_Snapshot error: cannot write file " << file_name << "\n";
		exit(0);
		return false;
	}
	Write_Binary(output);
	output.close();
	return true;
}

template<int d, typename T>
bool Points<d, T>::Load_Snapshot(const std::string& file_name)
{
	std::ifstream input(file_name, std::ios::binary);
	if (!input) {
		std::cerr << "Points::Load_Snapshot error: cannot read file " << file_name << "\n";
		exit(0);
		return false;
	}
	Read_Binary(input);
	input.close();
	return true;
}

//template<int d, typename T> int Points<d, T>::Add_Point(VectorD pos) {
//	X()->push_back(pos);
//	int n = Size();
//	this->Resize(n);//Critical: use this to access child class
//	return n - 1;
//}
//
//template<int d, typename T> int Points<d, T>::Add_Points(const Array<VectorD>& vertices) {
//	int n=(int)vertices.size();
//	int idx = this->Add_Elements(n);
//	for (int i = 0; i < n; i++) {
//		this->X(i + idx) = vertices[i];}
//	return idx;
//}
//
template<int d, typename T> void Points<d, T>::Write_To_File_3d(const std::string& file_name) const
{
	if constexpr (d == 3) {
		File::Write_Binary_To_File(file_name, *this);
	}
	else {
		Points<3, T> p3; p3.Resize(Size());
		AuxFunc::Dim_Conversion_Array<T, d, 3>(*X(), *p3.X(), (T)0);
		File::Write_Binary_To_File(file_name, p3);
	}
}

//template<int d, typename T>
//void Points<d, T>::Write_To_File_3d_Fast(const std::string& file_name) const
//{
//	RenderFunc::Write_Points<d, T>(file_name, XRef());
//}

template class Points<2,double>;
template class Points<3,double>;
template class Points<2,float>;
template class Points<3,float>;

//////////////////////////////////////////////////////////////////////////
////tracker points

template<int d,typename T> void TrackerPoints<d,T>::Write_To_File_3d(const std::string& file_name) const
{
	if constexpr (d==3){
		File::Write_Binary_To_File(file_name,*this);}
	else{
		TrackerPoints<3,T> p3;p3.Resize(Size());
		AuxFunc::Dim_Conversion_Array<T,d,3>(*X(),*p3.X(),(T)0);
		AuxFunc::Dim_Conversion_Array<T,d,3>(*V(),*p3.V(),(T)0);
		*p3.I()=*I();
		File::Write_Binary_To_File(file_name,p3);}
}

template class TrackerPoints<2,double>;
template class TrackerPoints<3,double>;
template class TrackerPoints<2,float>;
template class TrackerPoints<3,float>;

//////////////////////////////////////////////////////////////////////////
////particles

template<int d,typename T> void Particles<d,T>::Write_To_File_3d(const std::string& file_name) const
{
	if constexpr (d==3){
		File::Write_Binary_To_File(file_name,*this);}
	else {
		Particles<3, T> p3; p3.Resize(Size());
		AuxFunc::Dim_Conversion_Array<T, d, 3>(*X(), *p3.X(), (T)0);
		AuxFunc::Dim_Conversion_Array<T, d, 3>(*V(), *p3.V(), (T)0);
		AuxFunc::Dim_Conversion_Array<T, d, 3>(*F(), *p3.F(), (T)0);
		*p3.M() = *M();
		*p3.C() = *C();
		*p3.I() = *I();
		File::Write_Binary_To_File(file_name, p3);
	}
}

template class Particles<2,double>;
template class Particles<3,double>;
template class Particles<2,float>;
template class Particles<3,float>;
