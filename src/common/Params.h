//////////////////////////////////////////////////////////////////////////
// A universal parameter class that can be used in all simulation systems
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#pragma once
#include <iostream>

class AnyObjectBase {
public:
};

template<class T>
class AnyObject :public AnyObjectBase {
public:
	std::shared_ptr<T> value = nullptr;
	AnyObject(){}
	AnyObject(const T& a) { value = std::make_shared<T>(a); }
	void operator = (const T &a) { value = std::make_shared<T>(a); }
	operator T() const { AuxFunc::Assert(value != nullptr, "AnyObject error: nullptr"); return *value; }
	T operator () (void) const { AuxFunc::Assert(value != nullptr, "AnyObject error: nullptr"); return *value; }
};

#define Declare_Param(T,a)              \
public: AnyObject<T> a;                 \

////call the parent virtual function first
#define Register_Params(...)															\
virtual void Register_Attributes(){auto vec=AuxFunc::Split_String(#__VA_ARGS__,", ");std::reverse(vec.begin(),vec.end());this->Register_Att(vec,__VA_ARGS__);}

class Params {
public:
	Params() {
	}
	std::map<std::string, std::shared_ptr<AnyObjectBase>> params_map;
	template<class T> void Register_Att(Array<std::string>& reverse_names, AnyObject<T>& att)
	{
		this->params_map[reverse_names.back()] = std::make_shared<AnyObject<T>>(att);
		reverse_names.pop_back();
	}
	template<typename T, typename... Args>
	void Register_Att(Array<std::string>& reverse_names, AnyObject<T>& att, Args & ...rest)
	{
		//std::cout << "register att: "; for (int i = 0; i < reverse_names.size(); i++) std::cout << reverse_names[i] << " "; std::cout << "\n";
		this->Register_Att<T>(reverse_names, att);
		this->Register_Att(reverse_names, rest...);
	}
	virtual void Register_Attributes() {}
};

class ParamsExample : public Params {
public:
	Declare_Param(real, cfl);
	Declare_Param(Vector3, g);
	Declare_Param(int, type);
	Register_Params(cfl, g, type);
	ParamsExample(real _cfl, Vector3 _g, int _type)
		:cfl(_cfl), g(_g)
	{
		type = _type;
		Register_Attributes();
	}
};