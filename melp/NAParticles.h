//////////////////////////////////////////////////////////////////////////
// Common base class for LagrangeParticles and EulerParticles
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __SSPHParticles_h__
#define __SSPHParticles_h__
#include "GeometryParticles.h"
#include "Kernels.h"
#include "NeighborSearcher.h"
#include "ImplicitShape.h"
#include "PointSet.h"
#include "AuxFunc.h"

class NeighborInfo {
private:
	Array<int> data;
public:
	real radius;
	NeighborInfo();
	void assign(real r, const Array<int>& nbs);
	int size(void)const;
	int operator [] (int k)const;
};

template<typename T>
class Function_Indexer {
public:
	std::function<T(const int)> f;
	Function_Indexer(std::function<T(const int)> _f) { f = _f; }
	T operator [] (int idx) const {
		return f(idx);
	}
};

template<int d> 
class NAParticles : public GeometryParticles<d> {
	Typedef_VectorDii(d);
	Typedef_MatrixD(d);
public:
	using Base = GeometryParticles<d>;
	using VectorT = Vector<real, d - 1>;
	using Base::Size;
	using Base::X;
	using Base::F;
	using Base::V;
	using Base::E;
	using Base::SA;
public:
	Declare_Attribute(real, NDen, nden);
	Declare_Attribute(int, Status, status);//0:normal, 1:to delete
	Declare_Attribute(real, Vol, vol);
	Declare_Attribute(real, H, h);
	Declare_Attribute(real, Soap, soap);
	Declare_Attribute(real, Gamma, gamma);
	Declare_Attribute(real, Den, den);
	Declare_Attribute(real, Temp, temp);
	Declare_Attribute(real, WM, wm);

	NAParticles() {
		New_Attributes();
	}

	Declare_Attribute_Inherent_Func(nden, status, vol, h, soap, gamma, den, temp, wm);

	const real EPS = 1e-8;//threshold used by avg weight
	std::shared_ptr<NeighborSearcher<d>> nbs_searcher;
	Array<real> phis;
	Array<VectorD> bnd_normals;

	////Basic operators
	VectorD World_Relative_Vec(int p, int q) const { return X(q) - X(p); }
	
	//store neighborhood results
	void Save_Neighbor_Infos(const Array<VectorD>& pos, Array<NeighborInfo>& infos, real radius);

	//Weighted AVG or SUM with kernels
	//Fill nbs_searcher for corresponding private functions
	template<class T> T World_Avg_Value(const VectorD& pos, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type, bool& has_nbs);
	template<class T> T World_Avg_Value_KNN(const VectorD& pos, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, const int k, KernelType kernel_type, bool& has_nbs);
	template<class T> T World_Sum_Value(const VectorD& pos, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type)const;
	MatrixD World_Avg_Value_Affine(const VectorD& pos, const Array<VectorD>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type)const;
	template<class T> T Surface_Sum_Value(const VectorD& pos, const MatrixD& frame, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type)const;
	VectorT Surface_Sum_Grad(const VectorD& pos, const MatrixD& frame, const Array<real>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type)const;
	template<class T> T Surface_Convolution_Value(const VectorD& pos, const MatrixD& frame, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type)const;
	template<class T> T World_Convolution_Value(const VectorD& pos, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type)const;

	template<typename T> T Surface_Laplacian(const VectorD& pos, const MatrixD& frame, const T& value, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type)const;
	template<typename T> T Surface_Laplacian(const VectorD& pos, const MatrixD& frame, const T& value, const std::function<T(const int)>& f, const KernelSPH& kernel, const real radius, KernelType kernel_type)const;
	real Surface_Density_Derivative(const VectorD& pos, const MatrixD& frame, const VectorD& value, const Array<VectorD>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type)const;
	real Surface_Density_Derivative_Weighted(const VectorD& pos, const MatrixD& frame, const VectorD& value, const Array<VectorD>& f_arr, const Array<real>& w_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type)const;
	real Surface_Number_Density_Derivative(const VectorD& pos, const MatrixD& frame, const VectorD& value, const Array<VectorD>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type)const;
	real Surface_Divergence(const VectorD& pos, const MatrixD& frame, const VectorD& value, const Array<VectorD>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type)const;
	template<typename T> VectorT Surface_Gradient_Symmetric(const VectorD& pos, const MatrixD& frame, const T& mirrored_f, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const;


	////Updates
	virtual void Update(void) = 0;
	
	//Boundary treatment
	void Update_Phis(const ImplicitShape<d>& analytical_boundary);
	//It will call Update_Phis
	void Correct_Velocity_With_Analytical_Boundary(const real dt, const real scale, const ImplicitShape<d>& analytical_boundary);
	void Correct_Position_With_Analytical_Boundary(const real dt, const real scale, const ImplicitShape<d>& analytical_boundary);

	////Diagnosis
	void Print_Dynamics(int i){ std::cout << "particle " << i << ",X=(" << X(i).transpose() << "),F=(" << F(i).transpose() << "),V=(" << V(i).transpose() << "),vel rate=" << V(i).norm() << "\n"; }
	void Print_Fastest_Particle(void);

	bool Is_Trivial(void) { return Size() < 1; }

	void Apply_Rotation(int idx, MatrixD R) {
		X(idx) = R * X(idx);
		V(idx) = R * V(idx);
		F(idx) = R * F(idx);
		E(idx) = R * E(idx);
	}

	void Apply_Translation(int idx, VectorD T) {
		X(idx) = X(idx) + T;
	}
};

template<int d>
void Check_Farest_Pairs(const NAParticles<d>& a, const NAParticles<d>& b, std::string pair_name) {
	real farest_dis = -1;
	int a_idx = -1, b_idx = -1;
	for (int i = 0; i < a.Size(); i++) {
		int j = b.nbs_searcher->Find_Nearest_Nb(a.X(i));
		real r = (a.X(i) - b.X(j)).norm();
		if (r > farest_dis) {
			farest_dis = r;
			a_idx = i;
			b_idx = j;
		}
	}
	std::cout << "Farest pair " << pair_name << " is: (" << a_idx << "," << b_idx << ") with distance " << farest_dis << "\n";
}

template<int d>
void Check_Farest_Pairs(const NAParticles<d>& a, std::string name) {
	real farest_dis = -1;
	int idx1 = -1, idx2 = -1;
	for (int i = 0; i < a.Size(); i++) {
		Array<int> nbs;
		a.nbs_searcher->Find_K_Nearest_Nb(a.X(i), 2, nbs);
		int j = nbs[0]; if (j == i) j = nbs[1];
		real r = (a.X(i) - a.X(j)).norm();
		if (r > farest_dis) {
			farest_dis = r;
			idx1 = i, idx2 = j;
		}
	}
	std::cout << "Farest pair " << name << " is: (" << idx1 << "," << idx2 << ") with distance " << farest_dis << "\n";
}


template<int d>
template<class T>
inline T NAParticles<d>::World_Avg_Value(const VectorD& pos, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type, bool& has_nbs)
{
	Array<int> nbs = nbs_searcher->Find_Neighbors(pos, radius);
	T sum = Zero<T>();
	real total_weight = 0;
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		VectorD wr_pj = X(j) - pos;
		real w = kernel.Weight<d>(wr_pj, radius, kernel_type);
		sum += f_arr[j] * w;
		total_weight += w;
	}
	if (total_weight > EPS) {
		has_nbs = true;
		return sum / total_weight;
	}
	else {
		has_nbs = false;
		return Zero<T>();
	}
}

template<int d>
template<class T>
inline T NAParticles<d>::World_Avg_Value_KNN(const VectorD& pos, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, const int k, KernelType kernel_type, bool& has_nbs)
{
	Array<int> nbs;
	int num_found = nbs_searcher->Find_K_Nearest_Nb(pos, k, nbs);
	//std::cout << "num found: " << num_found << std::endl;
	T sum = Zero<T>();
	real total_weight = 0;
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		VectorD wr_pj = X(j) - pos;
		//std::cout << "nearest dist: " << wr_pj.norm() << std::endl;
		real w = kernel.Weight<d>(wr_pj, radius, kernel_type);
		//std::cout << "w: " << w << std::endl;
		sum += f_arr[j] * w;
		total_weight += w;
		//std::cout << "f_arr j: " << f_arr[j] << std::endl;
	}
	if (total_weight > EPS) {
		has_nbs = true;
		return sum / total_weight;
	}
	else {
		has_nbs = false;
		return Zero<T>();
	}
}

template<int d>
template<class T>
inline T NAParticles<d>::World_Sum_Value(const VectorD& pos, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type) const
{
	Array<int> nbs = nbs_searcher->Find_Neighbors(pos, radius);
	T sum = Zero<T>();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		VectorD wr_pj = X(j) - pos;
		real w = kernel.Weight<d>(wr_pj, radius, kernel_type);
		sum += f_arr[j] * w;
	}
	return sum;
}


template<int d>
Matrix<real, d> NAParticles<d>::World_Avg_Value_Affine(const VectorD& pos, const Array<VectorD>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type)const
{
	Array<int> nbs = nbs_searcher->Find_Neighbors(pos, radius);
	MatrixD sum = MatrixD::Zero();
	real total_weight = 0;
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		VectorD wr_pj = X(j) - pos;
		real w = kernel.Weight<d>(wr_pj, radius, kernel_type);
		sum += f_arr[j] * (wr_pj.transpose()) * w;
		total_weight += w;
	}
	if (total_weight > EPS) {
		sum *= 1. / total_weight;
		return sum;
	}
	else {
		return MatrixD::Zero();
	}
}

template<int d>
template<class T>
inline T NAParticles<d>::Surface_Sum_Value(const VectorD& pos, const MatrixD& frame, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type) const
{
	Array<int> nbs = nbs_searcher->Find_Neighbors(pos, radius);
	T sum = Zero<T>();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		VectorD wr_pj = X(j) - pos;
		VectorT sr_pj = PointSet<d>::Rotate_To_TPlane(wr_pj, frame);
		real w = kernel.Weight<d-1>(sr_pj, radius, kernel_type);
		sum += f_arr[j] * w;
	}
	return sum;
}

template<int d>
Vector<real, d-1> NAParticles<d>::Surface_Sum_Grad(const VectorD& pos, const MatrixD& frame, const Array<real>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type) const
{
	Array<int> nbs = nbs_searcher->Find_Neighbors(pos, radius);
	VectorT sum = VectorT::Zero();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		VectorD wr_jp = pos - X(j);
		VectorT sr_jp = PointSet<d>::Rotate_To_TPlane(wr_jp, frame);
		VectorT grad_W = kernel.Grad<d - 1>(sr_jp, radius, kernel_type);
		sum += f_arr[j] * grad_W;
	}
	return sum;
}

template<int d>
template<class T>
inline T NAParticles<d>::Surface_Convolution_Value(const VectorD& pos, const MatrixD& frame, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type) const
{
	Array<int> nbs = nbs_searcher->Find_Neighbors(pos, radius);
	T sum = Zero<T>();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		VectorD wr_pj = X(j) - pos;
		VectorT sr_pj = PointSet<d>::Rotate_To_TPlane(wr_pj, frame);
		real w = kernel.Weight<d - 1>(sr_pj, radius, kernel_type);
		sum += f_arr[j] * w * SA(j);
	}
	return sum;
}

template<int d>
template<class T>
inline T NAParticles<d>::World_Convolution_Value(const VectorD& pos, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type) const
{
	Array<int> nbs = nbs_searcher->Find_Neighbors(pos, radius);
	T sum = Zero<T>();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		VectorD wr_pj = X(j) - pos;
		real w = kernel.Weight<d>(wr_pj, radius, kernel_type);
		sum += f_arr[j] * w * Vol(j);
	}
	return sum;
}


template<int d>
template<class T>
inline T NAParticles<d>::Surface_Laplacian(const VectorD& pos, const MatrixD& frame, const T& value, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type) const
{
	Array<int> nbs = nbs_searcher->Find_Neighbors(pos, radius);
	T lap = Zero<T>();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		real S_j = SA(j);
		VectorD wr_pj = X(j) - pos;
		VectorT sr_pj = PointSet<d>::Rotate_To_TPlane(wr_pj, frame);
		real norm_ji = std::max(sr_pj.norm(), radius * 0.001);
		VectorT grad_W = kernel.Grad<d - 1>(sr_pj, radius, kernel_type);
		T f_j = f_arr[j];
		lap += S_j * (f_j - value) * 2 * grad_W.norm() / norm_ji;
	}
	return lap;
}

template<int d>
template<class T>
inline T NAParticles<d>::Surface_Laplacian(const VectorD& pos, const MatrixD& frame, const T& value, const std::function<T(const int)>& f, const KernelSPH& kernel, const real radius, KernelType kernel_type) const
{
	Array<int> nbs = nbs_searcher->Find_Neighbors(pos, radius);
	T lap = Zero<T>();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		real S_j = SA(j);
		VectorD wr_pj = X(j) - pos;
		VectorT sr_pj = PointSet<d>::Rotate_To_TPlane(wr_pj, frame);
		real norm_ji = std::max(sr_pj.norm(), radius * 0.001);
		VectorT grad_W = kernel.Grad<d - 1>(sr_pj, radius, kernel_type);
		T f_j = f(j);
		lap += S_j * (f_j - value) * 2 * grad_W.norm() / norm_ji;
	}
	return lap;
}

template<int d>
inline real NAParticles<d>::Surface_Density_Derivative(const VectorD& pos, const MatrixD& frame, const VectorD& value, const Array<VectorD>& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	const auto& nbs = nbs_searcher->Find_Neighbors(pos, radius);
	real div = 0.;
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		real S_j = SA(j);
		real Vol_j = H(j) * S_j; // that is LVol, not Vol
		VectorD wr_jp = pos-X(j);
		VectorT sr_jp = PointSet<d>::Rotate_To_TPlane(wr_jp, frame);
		VectorT grad_W = kernel.Grad<d - 1>(sr_jp, radius, kernel_type);
		//VectorT fr_jp = PointSet<d>::Rotate_To_TPlane(value, frame) - PointSet<d>::Rotate_To_TPlane(f_arr[j], frame);
		VectorT fr_jp = PointSet<d>::Project_To_TPlane(value-f_arr[j], frame);
		div += Vol_j * fr_jp.dot(grad_W);
	}
	return div;
}

template<int d>
inline real NAParticles<d>::Surface_Density_Derivative_Weighted(const VectorD& pos, const MatrixD& frame, const VectorD& value, const Array<VectorD>& f_arr, const Array<real>& w_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	const auto& nbs = nbs_searcher->Find_Neighbors(pos, radius);
	real div = 0.;
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		real S_j = SA(j);
		real Vol_j = H(j) * S_j; // that is LVol, not Vol
		VectorD wr_jp = pos - X(j);
		VectorT sr_jp = PointSet<d>::Rotate_To_TPlane(wr_jp, frame);
		VectorT grad_W = kernel.Grad<d - 1>(sr_jp, radius, kernel_type);
		//VectorT fr_jp = PointSet<d>::Rotate_To_TPlane(value, frame) - PointSet<d>::Rotate_To_TPlane(f_arr[j], frame);
		VectorT fr_jp = PointSet<d>::Project_To_TPlane(value - f_arr[j], frame);
		div += Vol_j * fr_jp.dot(grad_W) * std::max<real>(0., 1.-w_arr[j]);
	}
	return div;
}

template<int d>
inline real NAParticles<d>::Surface_Number_Density_Derivative(const VectorD& pos, const MatrixD& frame, const VectorD& value, const Array<VectorD>& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	const auto& nbs = nbs_searcher->Find_Neighbors(pos, radius);
	real div = 0.;
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		VectorD wr_jp = pos - X(j);
		VectorT sr_jp = PointSet<d>::Rotate_To_TPlane(wr_jp, frame);
		VectorT grad_W = kernel.Grad<d - 1>(sr_jp, radius, kernel_type);
		VectorT fr_jp = PointSet<d>::Project_To_TPlane(value - f_arr[j], frame);
		//VectorT fr_jp = PointSet<d>::Rotate_To_TPlane(value, frame) - PointSet<d>::Rotate_To_TPlane(f_arr[j], frame);
		div += fr_jp.dot(grad_W);
	}
	return div;
}

template<int d>
inline real NAParticles<d>::Surface_Divergence(const VectorD& pos, const MatrixD& frame, const VectorD& value, const Array<VectorD>& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	const auto& nbs = nbs_searcher->Find_Neighbors(pos, radius);
	real div = 0.;
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		real S_j = SA(j);
		VectorD wr_jp = pos - X(j);
		VectorT sr_jp = PointSet<d>::Rotate_To_TPlane(wr_jp, frame);
		VectorT grad_W = kernel.Grad<d - 1>(sr_jp, radius, kernel_type);
		VectorT fr_jp = PointSet<d>::Project_To_TPlane(f_arr[j]-value, frame);
		//VectorT fr_jp = PointSet<d>::Rotate_To_TPlane(f_arr[j], frame) - PointSet<d>::Rotate_To_TPlane(value, frame);
		div += S_j * fr_jp.dot(grad_W);
	}
	return div;
}

template<int d>
template<typename T>
inline Vector<real, d - 1> NAParticles<d>::Surface_Gradient_Symmetric(const VectorD& pos, const MatrixD& frame, const T& mirrored_f, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	//return grad_f * rho_i;
	const auto& nbs = nbs_searcher->Find_Neighbors(pos, radius);
	VectorT grad_f = VectorT::Zero();

	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		VectorD wr_jp = pos - X(j);
		VectorT sr_jp = PointSet<d>::Rotate_To_TPlane(wr_jp, frame);
		VectorT grad_W = kernel.Grad<d - 1>(sr_jp, radius, kernel_type);
		grad_f += SA(j) * mirrored_f * grad_W;
	}

	return grad_f;
}

#endif

