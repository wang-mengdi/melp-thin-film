//#####################################################################
// Functional Programming tools for particle system
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################
#pragma once

#include "Common.h"

//I for index
template<class T> using PairIFunc = std::function<T(const int i, const int j)>;
template<class T> using SingleIFunc = std::function<T(const int i)>;
//S for scalar
using PairIFuncS = PairIFunc<real>;
using SingleIFuncS = SingleIFunc<real>;
//V for vector
template<int d> using PairIFuncV = PairIFunc<Vector<real, d> >;
template<int d> using SingleIFuncV = SingleIFunc<Vector<real, d> >;

template<typename T1, typename T2>
decltype(auto) operator * (SingleIFunc<T1> func, const T2 a) {
	return [=](const int i){
		return func(i) * a;
	};
}
template<typename T1, typename T2> decltype(auto) operator * (const T1 a, SingleIFunc<T2> func) { return func * a; }
template<typename T1, typename T2>
decltype(auto) operator * (PairIFunc<T1> func, const T2 a) {
	return [=](const int i, const int j){
		return func(i, j) * a;
	};
}
template<typename T1, typename T2> decltype(auto) operator * (const T1 a, PairIFunc<T2> func) { return func * a; }