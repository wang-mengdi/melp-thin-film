#ifndef __EulerParticles_h__
#define __EulerParticles_h__
#include "NAParticles.h"

template<int d> class EulerParticles : public NAParticles<d>
{
	Typedef_VectorDii(d);
	Typedef_MatrixD(d);
public:
	using Base = NAParticles<d>;
	using VectorT = Vector<real, d - 1>;					////tangential vector type
	using MatrixT = Matrix<real, d - 1>;					////tangential vector type

	using Base::nbs_searcher;
	
	using Base::X;
	using Base::V;
	using Base::SA;

	using Base::Normal;
	using Base::Size;
	//using Base::Surface_Laplacian;

	EulerParticles() {
		New_Attributes();
		//std::cout << "EulerParticles initialized\n";
		//Points<d>::Print_Attributes();
	}

	////Particle attributes
	Declare_Attribute(int, B, b); ////1 for boundary particles

	Declare_Attribute(real, LM, lm); ////1 for boundary particles
	Declare_Attribute(real, LVol, lvol);
	Declare_Attribute(VectorD, LMom, lmom);
	Declare_Attribute(VectorD, LAMom, lamom); // note that the affine momentum is in-plane only

	Declare_Attribute(VectorD, LV_T, lv_t);
	Declare_Attribute(VectorD, LV, lv);

	Declare_Attribute(real, KH, kh);     ////laplacian of height field
	Declare_Attribute(real, KS, ks);     ////laplacian of surface
	Declare_Attribute(real, ST, st); //surface tension coefficient
	Declare_Attribute(real, P, p);	
	Declare_Attribute(VectorD, FauxV, fauxv);
	Declare_Attribute(VectorD, KSN, ksn);
	Declare_Attribute(int, Ext, ext); // if is external
	Declare_Attribute(int, Inside_Another, inside_another); // if is external
	Declare_Attribute(real, BH_Ratio, bh_ratio); // if is external
	Declare_Attribute(VectorD, BH_LV, bh_lv); // if is external
	

	Declare_Attribute_Inherent_Func(b, lm, lvol, lmom, lamom, lv_t, lv, kh, ks, st, p, fauxv, ksn, ext, inside_another, bh_ratio, bh_lv);	////does not need the attributes from the base class

	PointSet<d> surface;
	//KernelSPH<d - 1> sph_kernel;
	real base_mass = 1.;
	real base_rho = 1.;

	void Initialize(real dx, int np_on_h) {
		surface.Initialize(dx, np_on_h, this);
		//sph_kernel.Initialize(surface.t_r);
		nbs_searcher = surface.nbs_searcher;
		Update();
	}

	//set:
	//base_thickness,base_density
	//VL,Gamma,B
	void Set_Values(real _base_rho, real _base_mass) {
		int n = Size();
		base_mass = _base_mass;
		base_rho = _base_rho;
#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			V(i) = VectorD::Zero();
			FauxV(i) = VectorD::Zero();
			LMom(i) = VectorD::Zero();
			LAMom(i) = VectorD::Zero();
			LV_T(i) = VectorD::Zero();
			LV(i) = VectorD::Zero();
			this->Gamma(i) = 1.0;
			Ext(i) = 1;
			Inside_Another(i) = 0;
			BH_Ratio(i) = 0.;
			BH_LV(i) = VectorD::Zero();
		}
	}

	virtual void Update(void) {
		surface.Update();
	}

	VectorT World_Vec_To_Surface(int idx, const VectorD& v_world)const { return PointSet<d>::Rotate_To_TPlane(v_world, this->E(idx)); }//Base::World_Vec_To_Surface(v_world, E(idx)); }
	VectorT Surface_Relative_Vec(int p, int q) const { return World_Vec_To_Surface(q, Base::World_Relative_Vec(p, q)); }
	VectorD Surface_Vector_To_World(const VectorT& v_surface, MatrixD frame, MatrixT metric) {
		//VectorT scaled_v_surface = metric.inverse() * v_surface;
		VectorT scaled_v_surface = v_surface;
		VectorD v_world; surface.Unproject_To_World(scaled_v_surface, frame, v_world);
		return v_world;
	}

	//Math operators
	//template<typename ArrayT> real Surface_Divergence(const int i, const ArrayT& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const;
	template<typename ArrayT> VectorT Surface_Gradient_Symmetric(const int i, const ArrayT& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const;
	VectorT Surface_Gradient_Symmetric(const int i, const std::function<real(const int)> &f, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const { return Surface_Gradient_Symmetric(i, Function_Indexer<real>(f), kernel, radius, kernel_type); }
	template<typename ArrayT> VectorT Surface_Gradient_Difference(const int i, const ArrayT& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const;
	VectorT Surface_Gradient_Difference(const int i, const std::function<real(const int)> &f, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const { return Surface_Gradient_Difference(i, Function_Indexer<real>(f), kernel, radius, kernel_type); }
	template<typename ArrayT> VectorT Surface_Gradient(const int i, const ArrayT& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const;
	Eigen::Matrix<real,d,d-1> Surface_Jacobian_Difference(const int i, const Array<VectorD>& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const;

	template<class T> T Surface_Laplacian(const int i, const T f_i, const Array<int>& nbs, const Array<T>& nbv, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const;
	template<class T> T Surface_Laplacian(const int i, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const;
	template<class T> T Surface_Laplacian(const int i, const std::function<T(const int)>& f, const KernelSPH& kernel, const real radius, const KernelType kernel_type)const;

	template<class T> T XSPH_Smoothing(const real alpha, Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type, std::function<bool(const int)> filter = nullptr)const;

	//Updates
	void Diffuse_Scalar(Array<real>& f, const real& dt, const real& coeff, const KernelSPH& kernel, const real radius, const KernelType &kernel_type);
	void Update_C_Metric_Tensors(const KernelSPH& kernel, const real radius);
	void Update_H_Curvatures(const KernelSPH& kernel, const real radius);
	void Update_C_Curvatures(const KernelSPH& kernel, const real radius);
};


template<int d>
template<typename ArrayT>
inline Vector<real, d - 1> EulerParticles<d>::Surface_Gradient_Symmetric(const int i, const ArrayT& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	const auto& nbs = nbs_searcher->Find_Neighbors(X(i), radius);
	VectorT grad_f = VectorT::Zero();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		if (i == j) continue;
		//VectorD wr_ji = Base::World_Relative_Vec(j, i);
		VectorT lr_ji = Surface_Relative_Vec(j, i);
		VectorT grad_W = kernel.Grad<d - 1>(lr_ji, radius, kernel_type);
		grad_f += SA(j) * ((f_arr[i] + f_arr[j]) / 2) * grad_W;
	}

	return grad_f;
}

template<int d>
template<typename ArrayT>
inline Vector<real, d - 1> EulerParticles<d>::Surface_Gradient_Difference(const int i, const ArrayT& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	const auto& nbs = nbs_searcher->Find_Neighbors(X(i), radius);
	VectorT grad_f = VectorT::Zero();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		real S_j = SA(j);
		VectorT lr_ji = Surface_Relative_Vec(j, i);
		VectorT grad_W = kernel.Grad<d-1>(lr_ji, radius, kernel_type);
		grad_f += S_j * (f_arr[j] - f_arr[i]) * grad_W;
	}

	return grad_f;
}

template<int d>
template<typename ArrayT>
inline Vector<real, d - 1> EulerParticles<d>::Surface_Gradient(const int i, const ArrayT& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	const auto& nbs = nbs_searcher->Find_Neighbors(X(i), radius);
	VectorT grad_f = VectorT::Zero();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		real S_j = SA(j);
		VectorT lr_ji = Surface_Relative_Vec(j, i);
		VectorT grad_W = kernel.Grad<d - 1>(lr_ji, radius, kernel_type);
		grad_f += S_j * (f_arr[j]) * grad_W;
	}
	return grad_f;
}

template<int d>
inline Eigen::Matrix<real, d, d - 1> EulerParticles<d>::Surface_Jacobian_Difference(const int i, const Array<VectorD>& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	const auto& nbs = nbs_searcher->Find_Neighbors(X(i), radius);
	Eigen::Matrix<real,d,d-1> grad_f = Eigen::Matrix<real, d, d - 1>::Zero();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		real S_j = SA(j);
		VectorT lr_ji = Surface_Relative_Vec(j, i);
		VectorT grad_W = kernel.Grad<d - 1>(lr_ji, radius, kernel_type);
		grad_f += S_j * (f_arr[j] - f_arr[i]) * grad_W.transpose();
	}

	return grad_f;
}

template<int d>
template<class T>
inline T EulerParticles<d>::Surface_Laplacian(const int i, const T f_i, const Array<int>& nbs, const Array<T>& nbv, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	T lap = Zero<T>();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		if (i == j) continue;
		real S_j = SA(j);
		VectorT lr_ji = Surface_Relative_Vec(j, i);
		real norm_ji = std::max(lr_ji.norm(), surface.t_dx * 0.01);
		VectorT grad_W = kernel.Grad<d-1>(lr_ji, radius, kernel_type);
		T f_j = nbv[k];
		lap += S_j * (f_j - f_i) * 2 * grad_W.norm() / norm_ji;
	}
	return lap;
}

template<int d>
template<class T>
inline T EulerParticles<d>::Surface_Laplacian(const int i, const Array<T>& f_arr, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	const auto& nbs = nbs_searcher->Find_Neighbors(X(i), radius);
	Array<T> nbv(nbs.size());
	for (int k = 0; k < nbs.size(); k++) { nbv[k] = f_arr[nbs[k]]; }
	return Surface_Laplacian(i, f_arr[i], nbs, nbv, kernel, radius, kernel_type);
}

template<int d>
template<class T>
inline T EulerParticles<d>::Surface_Laplacian(const int i, const std::function<T(const int)>& f, const KernelSPH& kernel, const real radius, const KernelType kernel_type) const
{
	const auto& nbs = nbs_searcher->Find_Neighbors(X(i), radius);
	Array<T> nbv(nbs.size());
	for (int k = 0; k < nbs.size(); k++) { nbv[k] = f(nbs[k]); }
	return Surface_Laplacian(i, f(i), nbs, nbv, kernel, radius, kernel_type);
}

template<int d>
template<class T>
inline T EulerParticles<d>::XSPH_Smoothing(const real alpha, Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type, std::function<bool(const int)> filter)const
{
	Array<T> smoothed_f; smoothed_f.resize(Size());
	AuxFunc::Fill(smoothed_f, Zero<T>());
#pragma omp parallel for
	for (int i = 0; i < Size(); i++) {
		if (filter && !filter(i)) {
			smoothed_f[i] = f_arr[i];
		}
		else {
			Array<int> nbs = nbs_searcher->Find_Neighbors(X(i), radius);
			T sum = Zero<T>();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				VectorD wr_pj = X(j) - X(i);
				VectorT sr_pj = PointSet<d>::Rotate_To_TPlane(wr_pj, this->E(i));
				real w = kernel.Weight<d - 1>(sr_pj, radius, kernel_type);
				sum += SA(j) * (f_arr[j] - f_arr[i]) * w;
			}
			smoothed_f[i] = f_arr[i] + alpha * sum;
		}
	}
#pragma omp parallel for
	for (int i = 0; i < Size(); i++) {
		f_arr[i] = smoothed_f[i];
	}
}

#endif
