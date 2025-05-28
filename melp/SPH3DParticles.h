//////////////////////////////////////////////////////////////////////////
// 3D Particles
// Copyright (c) (2018-), Yitong Deng, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __SPH3DParticles_h__
#define __SPH3DParticles_h__
#include "NAParticles.h"
#include "EulerParticles.h"
#include "ArrayFunc.h"



template<int d> class SPH3DParticles: public NAParticles<d>
{
	Typedef_VectorDii(d);
	Typedef_MatrixD(d);
	using Base= NAParticles<d>;
	using VectorT = Vector<real, d - 1>;
public:
	using Base::nbs_searcher;
	using Base::X;
	using Base::XRef;
	using Base::M;
	using Base::Vol;

	using Base::Size;
	using Base::World_Relative_Vec;

	////Particle attributes

	SPH3DParticles() {
		New_Attributes();
		nbs_searcher = std::make_shared<NeighborKDTree<d> >();
		nbs_searcher->Update_Points(XRef());
		//std::cout << "LagrangeParticles initialized\n";
		//Points<d>::Print_Attributes();
	}

	Declare_Attribute(VectorD, SN, sn);
	Declare_Attribute(real, NDen, nden);
	Declare_Attribute(real, Rho, rho);

	Declare_Attribute_Inherent_Func(sn, nden, rho);

	void Initialize() {
		Update();
	}

	virtual void Update(void) {
		nbs_searcher->Build_Data(XRef());
	}

	void Set_Values(real mass) {
		int n = Size();
#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			M(i) = mass / n;
		}
	}

	template<class D> D Laplacian(const int i, const D f_i, const Array<int>& nbs, const Array<D>& nbv, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const;
	template<class D> D Laplacian(const int i, const Array<D>& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const;
	template<class D> D Laplacian(const int i, const std::function<D(const int)>& f, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const;

	template<class ArrayD> VectorD Gradient_Difference(const int i, const ArrayD& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const;
	VectorD Gradient_Difference(const int i, const std::function<real(const int)>& f, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const { return Gradient_Difference(i, Function_Indexer<real>(f), kernel, radius, kernel_type); }
	template<class ArrayD> VectorD Gradient(const int i, const ArrayD& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const;
	VectorD Gradient(const int i, const std::function<real(const int)>& f, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const { return Gradient(i, Function_Indexer<real>(f), kernel, radius, kernel_type); }
};

template<int d>
template<class ArrayD>
inline Vector<real, d> SPH3DParticles<d>::Gradient_Difference(const int i, const ArrayD& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	const auto& nbs = nbs_searcher->Find_Neighbors(X(i), radius);
	VectorD grad_f = VectorD::Zero();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		real S_j = Vol(j);
		VectorD lr_ji = World_Relative_Vec(j, i);
		VectorD grad_W = kernel.Grad<d>(lr_ji, radius, kernel_type);
		grad_f += S_j * (f_arr[j] - f_arr[i]) * grad_W;
	}

	return grad_f;
}

template<int d>
template<class ArrayD>
inline Vector<real, d> SPH3DParticles<d>::Gradient(const int i, const ArrayD& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	const auto& nbs = nbs_searcher->Find_Neighbors(X(i), radius);
	VectorD grad_f = VectorD::Zero();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		real S_j = Vol(j);
		VectorD lr_ji = World_Relative_Vec(j, i);
		VectorD grad_W = kernel.Grad<d>(lr_ji, radius, kernel_type);
		grad_f += S_j * (f_arr[j]) * grad_W;
	}

	return grad_f;
}

template<int d>
template<class D>
inline D SPH3DParticles<d>::Laplacian(const int i, const D f_i, const Array<int>& nbs, const Array<D>& nbv, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	D lap = Zero<D>();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		if (i == j) continue;
		real S_j = Vol(j);
		VectorD lr_ji = -World_Relative_Vec(i, j);
		real norm_ji = std::max(lr_ji.norm(), 1e-8);
		VectorD grad_W = kernel.Grad<d>(lr_ji, radius, kernel_type);
		D f_j = nbv[k];
		lap += S_j * (f_j - f_i) * 2 * grad_W.norm() / norm_ji;
	}
	return lap;
}

template<int d>
template<class D>
inline D SPH3DParticles<d>::Laplacian(const int i, const Array<D>& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	const auto& nbs = nbs_searcher->Find_Neighbors(X(i), radius);
	Array<D> nbv(nbs.size());
	for (int k = 0; k < nbs.size(); k++) { nbv[k] = f_arr[nbs[k]]; }
	return Laplacian(i, f_arr[i], nbs, nbv, kernel, radius, kernel_type);
}

template<int d>
template<class D>
inline D SPH3DParticles<d>::Laplacian(const int i, const std::function<D(const int)>& f, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	const auto& nbs = nbs_searcher->Find_Neighbors(X(i), radius);
	Array<D> nbv(nbs.size());
	for (int k = 0; k < nbs.size(); k++) { nbv[k] = f(nbs[k]); }
	return Laplacian(i, f(i), nbs, nbv, kernel, radius, kernel_type);
}

#endif
