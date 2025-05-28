//////////////////////////////////////////////////////////////////////////
// Particle Systems that have further physics meaning other than just positions
// Please refer to docs/sph-dev-zh.md
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "Common.h"
#include "Particles.h"
#include "NeighborSearcher.h"
#include "Kernels.h"
#include "Params.h"
#include "ParticleLambda.h"

////IMPORTANT NOTE: there should be another class as the parent class of them.

class ParamsSPHParticles : public Params {
public:
	Declare_Param(real, h);
	Declare_Param(KernelType, mode);
	Register_Params(h, mode);
	ParamsSPHParticles(const real &_h, const KernelType& _mode) {
		h = _h;
		mode = _mode;
		Register_Attributes();
	}
};

template<int d>
class SPHParticles : public Points<d>
{
	Typedef_VectorDii(d);
public:
	using Base = Points<d>;
	using Base::X;
	SPHParticles(const ParamsSPHParticles &params) {
		New_Attributes();
		kernel = KernelSPH((real)params.h, (KernelType)params.mode);
		//nbs_searcher = std::make_shared<NeighborKDTree<d>>();
		nbs_searcher = std::make_shared<NeighborHashing<d, 50>>(params.h);
	}
	Declare_Attribute(real, M, m);//mass
	Declare_Attribute(VectorD, V, v);//velocity
	Declare_Attribute(real, Rho, rho);//mass density

	Declare_Attribute_Inherent_Func(m, v, rho);

	//Mass density is calculated exactly by this kernel (including weight function)
	KernelSPH kernel;
	std::shared_ptr<NeighborSearcher<d>> nbs_searcher;

	virtual void Update(void);
	//The VectorD that f takes is r_ij. Here we follow the notation that r_ij=r_i-r_j
	template<class T> T Neighbor_Sum(const int i, PairIFunc<T> f)const;
	template<class T> SingleIFunc<T> Neighbor_Sum_Wrapper(PairIFunc<T> term)const;
	//Parallel
	template<class Func> void Exec_Each(Func f);
	template<class T> void Calc_Each(SingleIFunc<T> f, Array<T>& arr);
	template<class T> void Add_Each(SingleIFunc<T> f, Array<T>& arr);
	//template<class T> void Calculate_All_Pair(PairIFunc<T> f, Array<T>& arr);
	//template<class T> void Add_All_Pair(PairIFunc<T> &f, Array<T>& arr);
	//template<class T> T Weighted_Sum(const VectorD& pos, const Array<T>& arr, KernelType mode = KernelType::NONE);
};


template<int d>
template<class T>
inline T SPHParticles<d>::Neighbor_Sum(const int i, PairIFunc<T> f)const
{
	const Array<int>& nbs = nbs_searcher->Neighbors(i);
	T sum = Zero<T>();
	for (const int &j : nbs) {
		T val = f(i, j);
		sum += val;
	}
	return sum;
}

template<int d>
template<class T>
inline SingleIFunc<T> SPHParticles<d>::Neighbor_Sum_Wrapper(PairIFunc<T> term) const {
	AuxFunc::Assert(term != nullptr, "SPHParticles<d>::Neighbor_Sum_Wrapper: passed a null function");
	SingleIFunc<T> func = [&, term](const int i) {
		const Array<int>& nbs = nbs_searcher->Neighbors(i);
		T sum = Zero<T>();
		for (const int &j : nbs) {
			T val = term(i, j);
			sum += val;
		}
		return sum;
	};
	return func;
}

template<int d>
template<class Func>
inline void SPHParticles<d>::Exec_Each(Func f)
{
	int N = this->Size();
#pragma omp parallel for
	for (int i = 0; i < N; i++) f(i);
}

template<int d>
template<class T>
inline void SPHParticles<d>::Calc_Each(SingleIFunc<T> f, Array<T>& arr)
{
	AuxFunc::Assert(f != nullptr, "Calc_Each: passed a null function");
	int N = this->Size();
#pragma omp parallel for
	for (int i = 0; i < N; i++) arr[i] = f(i);
}

template<int d>
template<class T>
inline void SPHParticles<d>::Add_Each(SingleIFunc<T> f, Array<T>& arr)
{
	AuxFunc::Assert(f != nullptr, "Add_Each: passed a null function");
	int N = this->Size();
#pragma omp parallel for
	for (int i = 0; i < N; i++) arr[i] += f(i);
}

template<int d> class NDSPHParticles : public Particles<d>
{
	Typedef_VectorDii(d); 
	Typedef_MatrixD(d); 
	using Base = Particles<d>;
public:

	using Base::Size; using Base::Resize;
	using Base::X; using Base::V; using Base::F; using Base::M;
	using Base::C;	///curvature

	NDSPHParticles() { New_Attributes(); }

	Declare_Attribute(real, ND, nrho);						////number density
	Declare_Attribute(real, Vol, vol);						////volume
	Declare_Attribute(real, P, p);							////pressure
	Declare_Attribute(VectorD, SF, sf);						////surface force
	Declare_Attribute(VectorD, SN, sn);						////surface normal

	Declare_Attribute_Inherent_Func(nrho, vol, p, sf, sn);	////does not need the attributes from the base class
};

