//////////////////////////////////////////////////////////////////////////
// Particles
// Copyright (c) (2018-),Bo Zhu, Mengdi Wang
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Particles_h__
#define __Particles_h__
#include <memory>
#include <map>
#include <set>
#include <fstream>
#include "Common.h"
#include "ArrayPointer.h"
#include "AuxFunc.h"

//////////////////////////////////////////////////////////////////////////
////macros to define helper functions for manipulating particle attributes
#define Declare_Attribute(T,A,a)																						\
	public:T& A(const int i){return (*a)[i];}																			\
	const T& A(const int i) const {return (*a)[i];}																		\
	Array<T>* A(){return a.get();}																						\
	const Array<T>* A() const {return a.get();}																			\
	ArrayPtr<T> A##Ptr(){return a;}																						\
	const ArrayPtr<T> A##Ptr() const {return a;}																		\
	Array<T>& A##Ref(){return *a;}																						\
	const Array<T>& A##Ref() const {return *a;}				  															\
	void Rebind_##A(Array<T>* _a){AuxFunc::Assert(_a->size()==this->Size(),"Points class rebine mismatch");a.reset(_a);}\
	std::function<T(const int )> A##_callable=[&](const int i){return A(i);};																\
    protected:ArrayPtr<T> a;public:

////parent functions
#define Declare_Attribute_Base_Func(...)											\
virtual void New_Attributes(){New_Att(__VA_ARGS__);auto vec=AuxFunc::Split_String(#__VA_ARGS__,", ");std::reverse(vec.begin(),vec.end());Register_Att(att_map,vec,__VA_ARGS__);}

////call the parent virtual function first
#define Declare_Attribute_Inherent_Func(...)															\
virtual void New_Attributes(){Base::New_Attributes();New_Att(__VA_ARGS__);auto vec=AuxFunc::Split_String(#__VA_ARGS__,", ");std::reverse(vec.begin(),vec.end());Register_Att(this->att_map,vec,__VA_ARGS__);}

//////////////////////////////////////////////////////////////////////////
////CODESAMPLE: variadic template functions for manipulating attributes
////Here is an example of of the usage of variadic tempalte
//template <typename T> void print(int n,const T &t){std::cout<<"n"<<n<<": "<<t;}
//template <typename T,typename...Args> void print(int n,const T &t,const Args&...rest)
//{std::cout<<n<<": "<<t<<" ";print(10,rest...);}

template<typename T> void New_Att(ArrayPtr<T> & att)
{if(att==nullptr)att=std::make_shared<Array<T> >();}
template<typename T,typename...Args>
void New_Att(ArrayPtr<T> & att,Args & ...rest)
{New_Att<T>(att);New_Att(rest...);}

template<class T> void Register_Att(std::map<std::string, std::shared_ptr<ArrayPointerBase>>& att_map, Array<std::string>& reverse_names, ArrayPtr<T>& att)
{
	att_map[reverse_names.back()] = std::make_shared<ArrayPointerDerived<T>>(att);
	reverse_names.pop_back();
}
template<typename T, typename... Args>
void Register_Att(std::map<std::string, std::shared_ptr<ArrayPointerBase>>& att_map, Array<std::string>& reverse_names, ArrayPtr<T>& att, Args & ...rest)
{
	Register_Att<T>(att_map, reverse_names, att); 
	Register_Att(att_map, reverse_names, rest...);
}

//////////////////////////////////////////////////////////////////////////
////particle classes
//////////////////////////////////////////////////////////////////////////

template<int d, typename T = real> class Points
{
	using Base = Points<d, T>;
	using VectorD = Vector<T, d>;
public:
	std::map<std::string, std::shared_ptr<ArrayPointerBase>> att_map;
	Points() { New_Attributes(); }
	virtual ~Points() {}

	Declare_Attribute(VectorD, X, x);						////position
	Declare_Attribute_Base_Func(x);
	void Print_Attributes(void);

	//NOTE: all insertion operations set all new elements to 0 as default.
	//If you are modifying, please maintain this assertion. - Mengdi
	int Size() const { return (int)(*x).size(); }
	void Resize(const int size);
	void Reserve(const int size);
	int Add_Element(void);
	int Add_Position_Only(const VectorD& pos);
	int Add_Duplicate_Element(int idx);//add an element which is the duplicate of idx. return the newly added index (equal to the number before adding)
	int Add_Elements(const int n);//return the size before adding
	int Join(const Points<d, T>& src);//join src to it. return the size before append
	int Delete_Elements(const Array<int>& is_deleted);//return the number of remaining elements, is_deleted[i]!=0 means a deletion.
	int Delete_Elements(const std::set<int>& to_delete);//return the number of remaining elements.
	int Delete_Elements_Safe(Array<int> is_deleted);//delete elements however use pass-value style, so you can use it to delete itself.
	int Filter_Elements(std::function<bool(const int)> filter_func);//only keep elements such that filter_func(i)==true. filter_func can call the Points itself.
	void Copy_Element_From(const int idx, const Points<d, T>& src, const int src_idx);
	void Copy_Element_From(const int idx, const Points<d, T>& src, const int src_idx, const Array<std::string> attrs);//only copy attributes in attrs

	////I/O functions
	virtual void Write_Binary(std::ostream& output) const;
	virtual void Read_Binary(std::istream& input);
	bool Save_Snapshot(const std::string& file_name = "points.bin")const;
	bool Load_Snapshot(const std::string& file_name = "points.bin");

	////lambda functions
	template<class Fvoid>
	void Exec_Each(Fvoid f)
	{
		int N = this->Size();
#pragma omp parallel for
		for (int i = 0; i < N; i++) f(i);
	}
	template<class F1int>
	void Calc_Each(F1int f, Array<decltype(f(0))>& arr)
	{
		int N = this->Size();
		arr.resize(N);
#pragma omp parallel for
		for (int i = 0; i < N; i++) arr[i] = f(i);
	}

public:
	////Old functions, may not extended to all derived classes
	//Add
	//virtual int Add_Point(VectorD pos);//return the index of new point
	//virtual int Add_Points(const Array<VectorD>& vertices);

	//IO
	virtual void Write_To_File_3d(const std::string& file_name) const;
	//virtual void Write_To_File_3d_Fast(const std::string& file_name) const;	////Fast write function writes positions to files with type of float
};

template<int d, typename T = real> class TrackerPoints : public Points<d, T>
{
	using VectorD = Vector<T, d>; using Base = Points<d, T>;
public:
	using Base::Size; using Base::Resize; using Base::X; using Base::x;
	TrackerPoints() { New_Attributes(); }

	////attributes
	Declare_Attribute(VectorD, V, v);			////velocity
	Declare_Attribute(short, I, idx);			////flag
	Declare_Attribute_Inherent_Func(v, idx);

	////IO
	virtual void Write_To_File_3d(const std::string& file_name) const;
};

template<int d, typename T = real> class Particles : public Points<d, T>
{
	using VectorD = Vector<T, d>; using Base = Points<d, T>;
public:
	using Base::Size; using Base::Resize; using Base::X;
	Particles() { New_Attributes(); }

	Declare_Attribute(VectorD, V, v);			////velocity
	Declare_Attribute(VectorD, F, f);			////force
	Declare_Attribute(T, M, m);					////mass
	Declare_Attribute(T, C, c);					////color
	Declare_Attribute(int, I, idx);				////index
	Declare_Attribute_Inherent_Func(v, f, m, c, idx);

	////IO
	virtual void Write_To_File_3d(const std::string& file_name) const;
};

#endif
