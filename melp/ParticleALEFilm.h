//////////////////////////////////////////////////////////////////////////
// ParticleALEFilm
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __ParticleALEFilm_h__
#define __ParticleALEFilm_h__
#include "PointSet.h"
#include "Fluid3DSPH.h"
#include "RandomNumber.h"
#include "ArrayIO.h"
#include "LagrangeParticles.h"
#include "EulerParticles.h"
#include "ALEParams.h"
#include "ArrayFunc.h"
#include "Timer.h"

#include <numeric>      // std::iota
#include <iomanip>
#include <fstream>
#include <omp.h>

template<int d>
struct BlackHoleSeeder {
	Typedef_VectorDii(d);
	
	VectorD loc = VectorD::Zero();
	real r = 0.;
	real period = 100.;
	real phase = 0.;
	real threshold = 0.;
	bool initialized = false;
	BlackHoleSeeder()
	{
	}
	BlackHoleSeeder(VectorD _loc, real _r, real _period, real _phase, real _threshold)
	{
		loc = _loc;
		r = _r;
		period = _period;
		phase = _phase;
		threshold = _threshold;
		initialized = true;
	}
};


template<int d> 
class ParticleALEFilm
{
	Typedef_VectorDii(d);
	Typedef_MatrixD(d);
public:
	using VectorT = Vector<real, d - 1>;
	using MatrixT = Matrix<real, d - 1>;
	Fluid3DSPH<d>* fluid_3d = nullptr;
	LagrangeParticles<d> l_particles;
	EulerParticles<d> e_particles;

	real curr_nden_err = 0.;

	real simulation_scale = 1.; // how big the simulated object is 
	real simulation_fineness = 1.; // how fine the simulated object is 

	real max_vel = std::numeric_limits<real>::max();
	real max_vel_e = std::numeric_limits<real>::max();
	real max_vel_el = std::numeric_limits<real>::max();
	real max_vel_faux_e = std::numeric_limits<real>::max();
	real max_vel_l = std::numeric_limits<real>::max();
	real max_vel_3d = std::numeric_limits<real>::max();

	Array<real> ones_e;
	Array<real> ones_b;
	Array<real> ones_l;
	Array<VectorD> temp_LV;
	Array<real> temp_eh;
	Array<real> eh;
	Array<real> temp_lh;

	Array<VectorD> plateau_normals;

	Array<VectorD> temp_vector;
	Array<real> temp_scalar;
	Array<real> is_rim_confidence;
	Array<real> near_rim_confidence;
	Array<int> to_delete_e;
	Array<int> to_delete_l;

	real avg_gamma = 1.;
	real avg_h = 1.;
	real avg_NDen = 1.;
	real fluid_volume = 1.;
	real fluid_area = 1.;
	real fluid_soap = 1.;
	real avg_l_sa = 1.;
	real avg_e_sa = 1.;

	real e_evenness = std::numeric_limits<real>::max();

	// the two boundaries should represent the same geometric thing
	LagrangeParticles<d> boundary; //boundary is Lagrange particle
	ImplicitShape<d> analytical_boundary;

	VectorD g = 0. * VectorD::Unit(1) * (real)-1.;	////gravity
	VectorD faux_g = 0. * VectorD::Unit(1) * (real)-1.;	////gravity

	NeighborParams<d> neighbor_params;
	ChemicalParams chem_params;
	NumericalParams numeric_params;

	// Center
	VectorD init_COM;
	VectorD curr_COM;
	VectorD curr_COM_velocity = VectorD::Zero();

	//external force func
	std::function<VectorD(const int)> ext_acc_func = nullptr;
	std::function<VectorD(const int)> wind_func = nullptr;
	std::function<real(const int)> temp_func = nullptr;

	//simulation info
	std::string output_dir;
	int curr_frame = 0;
	real cfl;
	real frame_rate;
	real bursting_frame_rate;
	bool reinit_l_geometry = true;
	int verbose = 0;
	bool consider_black_hole = false;
	int l_geometry_mode = 0; // 0 means compute divergence, then times own height. 1 means compute density derivative
	bool l_geometry_sync = true; // if sync L height with E height
	real l_geometry_sync_strength = 0.2;
	bool delete_particles = false;
	bool evaporate = false;
	bool fixed_surface = false;
	bool flat_surface = false;

	bool e_only = false;
	
	//enclosed geometry info
	int enclosed = 0;
	real enclosed_vol = 0.;
	real enclosed_amount = 0.;
	real outer_pressure = 0.; //p_atm

	//advection scheme
	std::string exp_mode = "default";

	//some auxiliary switches
	int max_even_out_e_steps = 3;
	bool use_3d = false;
	Array<VectorD> rim_normals;
	bool first_advance = true;
	bool bursting = false;
	int num_E_preburst = 0;
	real dilution_rate_3d = 1.;
	bool camera_moving_along = false;
	bool euler_feels_gravity = true;
	bool enforce_momentum_conserve = false;
	// tools to seed black holes
	Array<BlackHoleSeeder<d>> bh_seeders;

	void Initialize(const real _cfl, const real _frame_rate, std::string _output_dir)
	{	
		double begin_time, end_time;
		if (verbose) {
			std::cout << "[Initialization]: fluid init begin" << std::endl;
			begin_time = omp_get_wtime();
		}
		cfl = _cfl;
		frame_rate = _frame_rate;
		bursting_frame_rate = 15 * frame_rate;
		output_dir = _output_dir;
		l_particles.Initialize();
		e_particles.Initialize(neighbor_params.e_dx, neighbor_params.sph_radius/neighbor_params.e_dx);
		boundary.Initialize();

		Reinit_Arrays();

		Update_E();
		Project_L(0.); //make sure that L particles are on the surface of E particles
		Update_L(); 

		num_E_preburst = e_particles.Size();

		// This part is done so that the boundary particles will be given the appropriate volume and mass, for L2E
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) {
			// Compute the number density of each E particle
			// Must remember to account for the boundary as well, so that the SA for particles near the boundary is not overly large
			e_particles.NDen(i) = e_particles.Surface_Sum_Value(e_particles.X(i), e_particles.E(i), ones_e, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC)
				+ boundary.Surface_Sum_Value(e_particles.X(i), e_particles.E(i), ones_b, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
			e_particles.SA(i) = 1. / e_particles.NDen(i);
		}
		avg_NDen = AuxFunc::Mean<real>(e_particles.NDenRef());

#pragma omp parallel for
		for (int i = 0; i < boundary.Size(); i++) {
			// this step is necessary for L2E, we need to know the correct volume of boundary particles, so we need to know the SA
			boundary.SA(i) = 1. / avg_NDen;
			boundary.H(i) = fluid_volume / fluid_area;
			boundary.Gamma(i) = fluid_soap / fluid_area;
			boundary.Vol(i) = boundary.SA(i) * boundary.H(i);
			boundary.Soap(i) = boundary.SA(i) * boundary.Gamma(i);
		}

		L2E(false, false); //this step is necessary as mass and volume must be transferred to E particles, in order to compute E particle geometry

		Update_E_Geometry();
		Update_L_Geometry(0.);

		init_COM = Compute_COM();

		if (verbose) {
			end_time = omp_get_wtime();
			std::cout << "[Initialization]: fluid init done. Time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
		}
	}

	void Degenerate(void)
	{
		e_particles.Resize(0);
		l_particles.Resize(0);
		Update_E();
		Update_L();
		Reinit_Arrays();
	}

	//template<int d>
	void Numeric_Check_Vector(Array<VectorD>& vectors_to_check)
	{
		bool success = true;
		int pn = vectors_to_check.size();
		for (int i = 0; i < pn; i++) {
			for (int axis = 0; axis < d; axis++) {
				real val = vectors_to_check[i][axis];
				if (std::isnan(val) || !std::isfinite(val)) {
					std::cerr << "Error: Numerical Check fails for X(" << i << ")[" << axis << "]=" << val << "\n";
					success = false;
					break;
				}
			}
			if (!success) break;
		}
		if (!success) {
			exit(0);
		}
	}

	//template<int d>
	void Numeric_Check_Matrix(Array<Matrix<real, d>>& vectors_to_check)
	{
		bool success = true;
		int pn = vectors_to_check.size();
		for (int i = 0; i < pn; i++) {
			for (int axis = 0; axis < d; axis++) {
				for (int axis2 = 0; axis2 < d; axis2++) {
					real val = vectors_to_check[i](axis, axis2);
					if (std::isnan(val) || !std::isfinite(val)) {
						std::cerr << "Error: Numerical Check fails for X(" << i << ")[" << axis << ", " << axis2 << "]=" << val << "\n";
						success = false;
						break;
					}
				}
			}
			if (!success) break;
		}
		if (!success) {
			exit(0);
		}
	}

	//template<int d>
	void Numeric_Check_Matrix_Inverse(Array<Matrix<real, d>>& vectors_to_check)
	{
		bool success = true;
		int pn = vectors_to_check.size();
		for (int i = 0; i < pn; i++) {
			for (int axis = 0; axis < d; axis++) {
				for (int axis2 = 0; axis2 < d; axis2++) {
					real val = vectors_to_check[i].inverse()(axis, axis2);
					if (std::isnan(val) || !std::isfinite(val)) {
						std::cerr << "Error: Numerical Check fails for X(" << i << ")[" << axis << ", " << axis2 << "]=" << val << "\n";
						success = false;
						break;
					}
				}
			}
			if (!success) break;
		}
		if (!success) {
			exit(0);
		}
	}

	void Numeric_Check_Scalar(Array<real>& scalars_to_check)
	{
		bool success = true;
		int pn = scalars_to_check.size();
		for (int i = 0; i < pn; i++) {
			real val = scalars_to_check[i];
			if (std::isnan(val) || !std::isfinite(val)) {
				std::cerr << "Error: Numerical Check fails for X(" << i << ")=" << val << "\n";
				success = false;
				break;
			}
		}
		if (!success) {
			exit(0);
		}
	}


	// update max velocity
	void Update_Max_Velocity(void) 
	{
		Update_E_Max_Velocity();
		Update_L_Max_Velocity();
		Update_EL_Max_Velocity();

		max_vel = std::max(max_vel_e, max_vel_l);

		//std::cout << "max vel e: " << max_vel_e << std::endl;
		//std::cout << "max vel l: " << max_vel_l << std::endl;
		//std::cout << "max vel el: " << max_vel_el << std::endl;

		if (fluid_3d != nullptr) {
			Update_3D_Max_Velocity();
			max_vel = std::max(max_vel, max_vel_3d);
		}
	}
	// update E max velocity
	void Update_E_Max_Velocity(void)
	{
		max_vel_e = ArrayFunc::Largest_Norm(e_particles.VRef());
	}
	// update Faux E max velocity
	void Update_E_Max_Faux_Velocity(void)
	{
		max_vel_faux_e = ArrayFunc::Largest_Norm(e_particles.FauxVRef());
	}
	// update L max velocity
	void Update_L_Max_Velocity(void)
	{
		max_vel_l = ArrayFunc::Largest_Norm(l_particles.VRef());
	}
	// update L max velocity
	void Update_EL_Max_Velocity(void)
	{
		max_vel_el = ArrayFunc::Largest_Norm(e_particles.LVRef());
	}
	// update 3d max velocity
	void Update_3D_Max_Velocity(void)
	{
		max_vel_3d = ArrayFunc::Largest_Norm(fluid_3d->particles.VRef());
	}

	real Cohesion_Kernel(const real& _r, const real& _h) const {
		real r = _r / _h;
		real h = 1;
		if (0 <= r && 2 * r <= h) {//0~0.5h
			return 2 * pow((h - r) * r, 3) - pow(h, 6) / 64.0;
		}
		else if (2 * r > h && r <= h) {
			return pow((h - r) * r, 3);
		}
		else return 0;
	}

	//void Configure_BH_Seeders(void);
	void Detect_BH(void);
	void Update_E_Positions(const real dt);
	void Update_L_Positions(const real dt);
	void Project_L(const real dt);
	void Even_Out_E(const real _dt, const bool repacking, const int _max_steps);
	void Even_Out_L(const int max_steps);
	void Initial_Packing(const int e_steps, const int l_steps);
	void L2E(bool use_affine = true, bool config_boundary = true);
	void E2L(void);
	void Update_E_Geometry(void);
	void Update_C_Curvature(void);
	void Apply_Rim_Surface_Tension(real dt);
	void Update_L_Geometry(const real dt);
	void Update_Norm_Dynamics(const real dt);
	void Update_Dynamics(const real dt);
	void Update_Burst(const real dt);
	void Advance_E(const real dt, const int iter);
	void Advance_L(const real dt);
	void Update_Tang_Dynamics(const real dt, const bool internal_force_only = false);
	void Update_Dynamics_E_Only(const real dt);
	real Update_Faux_Tang_Dynamics(const real dt, const real elasticity = 0.);
	void Update_Capillary_Forces(void);
	void Compute_Rim_Normals(void);

	void Update_BH_Ratio(void);
	void Apply_BH_Forces(const real dt);
	void Apply_Ext_Forces(const real dt);
	VectorD Compute_COM(void);
	VectorD Compute_COM_Velocity(void);
	void Evaporate(const real dt);

	void Delete_Particles(void);
	void Reinit_Arrays(bool e_changed = true, bool l_changed = true, bool b_changed = true) {
		if (e_changed) {
			to_delete_e.resize(e_particles.Size());
			std::fill(to_delete_e.begin(), to_delete_e.end(), 0);
			ones_e.resize(e_particles.Size());
			std::fill(ones_e.begin(), ones_e.end(), 1.);
			temp_LV.resize(e_particles.Size());
			AuxFunc::Fill(temp_LV, VectorD::Zero());
			temp_eh.resize(e_particles.Size());
			eh.resize(l_particles.Size());
			temp_vector.resize(e_particles.Size());
			AuxFunc::Fill(temp_vector, VectorD::Zero());
			temp_scalar.resize(e_particles.Size());
			AuxFunc::Fill(temp_scalar, 0.);
			is_rim_confidence.resize(e_particles.Size());
			AuxFunc::Fill(is_rim_confidence, 0.);
			near_rim_confidence.resize(e_particles.Size());
			AuxFunc::Fill(near_rim_confidence, 0.);
			rim_normals.resize(e_particles.Size());
			AuxFunc::Fill(rim_normals, VectorD::Zero());
			plateau_normals.resize(e_particles.Size());
			AuxFunc::Fill(plateau_normals, VectorD::Zero());
		}
		if (l_changed) {
			to_delete_l.resize(l_particles.Size());
			std::fill(to_delete_l.begin(), to_delete_l.end(), 0);
			ones_l.resize(l_particles.Size());
			std::fill(ones_l.begin(), ones_l.end(), 1.);
			eh.resize(l_particles.Size());
			AuxFunc::Fill(eh, 0.);
			temp_lh.resize(l_particles.Size());
			AuxFunc::Fill(temp_lh, 0.);
		}
		if (b_changed) {
			ones_b.resize(boundary.Size());
			std::fill(ones_b.begin(), ones_b.end(), 1.);
		}
	}

	real Compute_Enclosed_Volume(void);
	void Update_Enclosed_Forces(void);
	void Configure_Boundary(void);

	bool Held_In_Place(void) { // if the film is being attached to a boundary
		return (analytical_boundary.Available() || boundary.Size() > 0); 
	}

	void Update_E(void);
	void Update_L(void);

	void Update_WM(real radius);

	template<class T> T XSPH_Smoothing(const real alpha, Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type, std::function<bool(const int)> filter = nullptr)const; //modify the original f_arr
	template<class T> T XSPH_Smoothed_Value(const VectorD& pos, const MatrixD& frame, T& value, Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type)const; //return the velocity difference
	template<class T> T Boundary_Adhesion(Array<T>& f_arr) const; //boundary_adhesion
	VectorD XSPH_Smoothing_Tang_Norm(const real alpha_tang, const real alpha_norm, Array<VectorD>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type, std::function<bool(const int)> filter = nullptr)const; //modify the original f_arr
	void Decay_Velocity(Array<VectorD>& vels, const real dt, const real tang_strength, const real norm_strength);

	bool Eulerian_Near_Boundary(const int idx, const real dist) const
	{
		if (analytical_boundary.Available()) {
			return e_particles.phis[idx] > -dist;
		}
		else {
			return false;
		}
	}

	bool Is_Trivial(void) const
	{
		return (e_particles.Size() < 1 && l_particles.Size() < 1);
	}

	template<typename T> void Distribute_L2E(Array<T>& e_arr, const Array<T>& l_arr, real interp_radius, KernelType type = KernelType::QUINTIC);
	template<typename T> void Distribute_B2E(Array<T>& e_arr, const Array<T>& b_arr, real interp_radius, KernelType type = KernelType::QUINTIC); //distribute boundary quantities
	void AffineDistribute_L2E(Array<VectorD>& e_arr, const Array<MatrixT>& l_arr, real interp_radius, KernelType type = KernelType::QUINTIC);

	virtual void Advance(const real dt, const real current_time, const int iter)
	{	
		if (e_only) {
			std::cout << "Advancing E Only" << std::endl;
			Update_E_Geometry();
			Update_Dynamics_E_Only(dt);
			Advance_E(dt, iter);
			Update_Max_Velocity();
		}
		else {
			if (bursting) Compute_Rim_Normals();
#pragma omp parallel for
			for (int i = 0; i < e_particles.Size(); i++) {
				temp_LV[i] = e_particles.LV(i);
			}
			L2E(); // transfer momentum, mass, etc
#pragma omp parallel for
			for (int i = 0; i < e_particles.Size(); i++) {
				temp_LV[i] += e_particles.LV(i);
				temp_LV[i] *= 0.5;
				//if (near_rim_confidence[i] > 0.) temp_LV[i] *= 0.;
			}
			Update_E_Geometry();
			Update_L_Geometry(dt);
			Update_Dynamics(dt);
			E2L(); // transfer velocity

			Advance_E(dt, iter);
			Advance_L(dt);
			Update_Max_Velocity();
		}
	}

	//// test E only
	//virtual void Advance(const real dt, const real current_time, const int iter)
	//{
	//	Update_E_Geometry();

	//	////if we just want to test e redistribution
	//	//Even_Out_E(1., true, 1);

	//	//if we want to advance E without L
	//	Apply_Ext_Forces(dt);
	//	Update_Norm_Dynamics(dt);
	//	//smooth out e_v
	//	XSPH_Smoothing_Tang_Norm(numeric_params.XSPH_V_Tang, numeric_params.XSPH_V_Norm, e_particles.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
	//	for (int i = 0; i < numeric_params.XSPH_V_Passes; i++) {
	//		XSPH_Smoothing_Tang_Norm(0., 0.999, e_particles.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
	//	}


	//	Advance_E(dt, iter);

	//	Update_Max_Velocity();
	//	Decay_Velocity(e_particles.VRef(), dt, numeric_params.Velocity_Decay, numeric_params.Velocity_Decay);
	//}





	//// test E Burst only
	//virtual void Advance(const real dt, const real current_time, const int iter)
	//{
	//	std::cout << "Is bursting? " << bursting << std::endl;
	//	Update_E_Geometry();

	//	////if we just want to test e redistribution
	//	//Even_Out_E(1., true, 1);

	//	//if we want to advance E without L
	//	Apply_Ext_Forces(dt);
	//	if (bursting) {
	//		std::cout << "Updating Burst" << std::endl;
	//		std::cout << "DT is: " << dt << std::endl;
	//		Update_Burst(dt);
	//	}
	//	else {
	//		std::cout << "Updating Norm Dynamics" << std::endl;
	//		std::cout << "DT is: " << dt << std::endl;
	//		Update_Norm_Dynamics(dt);
	//	}
	//	//smooth out e_v
	//	XSPH_Smoothing_Tang_Norm(numeric_params.XSPH_V_Tang, numeric_params.XSPH_V_Norm, e_particles.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
	//	for (int i = 0; i < numeric_params.XSPH_V_Passes; i++) {
	//		XSPH_Smoothing_Tang_Norm(0., 0.999, e_particles.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
	//	}


	//	Advance_E(dt, iter);

	//	Update_Max_Velocity();
	//	Decay_Velocity(e_particles.VRef(), dt, numeric_params.Velocity_Decay, numeric_params.Velocity_Decay);
	//	
	//	if (delete_particles) Delete_Particles();
	//}
};



template<int d>
template<typename T>
inline void ParticleALEFilm<d>::Distribute_L2E(Array<T>& e_arr, const Array<T>& l_arr, real interp_radius, KernelType type)
{
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		Array<int> nbs = l_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), interp_radius);
		for (int k = 0; k < nbs.size(); k++) {
			int l_nb = nbs[k];
			VectorD wr_pj = (l_particles.X(l_nb) - e_particles.X(i));
			e_arr[i] += l_arr[l_nb] * neighbor_params.kernel.template Weight<d>(wr_pj, interp_radius, type) * l_particles.WM(l_nb);
		}
	}
}

template<int d>
template<typename T>
inline void ParticleALEFilm<d>::Distribute_B2E(Array<T>& e_arr, const Array<T>& b_arr, real interp_radius, KernelType type)
{
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		Array<int> nbs = boundary.nbs_searcher->Find_Neighbors(e_particles.X(i), interp_radius);
		T amount_received = Zero<T>();
		for (int k = 0; k < nbs.size(); k++) {
			int b_nb = nbs[k];
			VectorD wr_pj = (boundary.X(b_nb) - e_particles.X(i));
			e_arr[i] += b_arr[b_nb] * neighbor_params.kernel.template Weight<d>(wr_pj, interp_radius, type) * boundary.WM(b_nb);
		}
	}
}

template<int d>
inline void ParticleALEFilm<d>::AffineDistribute_L2E(Array<VectorD>& e_arr, const Array<MatrixT>& l_arr, real interp_radius, KernelType type)
{
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		if (bursting
			&& near_rim_confidence[i] > 0.) {
			//pass
		}
		else {
			Array<int> nbs = l_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), interp_radius);
			for (int k = 0; k < nbs.size(); k++) {
				int l_nb = nbs[k];
				VectorD wr_pj = (e_particles.X(i) - l_particles.X(l_nb));
				VectorT lr_pj = PointSet<d>::Project_To_TPlane(wr_pj, l_particles.E(i));
				VectorT mom_2 = l_arr[l_nb] * (lr_pj);
				VectorD mom_3;
				e_particles.surface.Unproject_To_World(mom_2, l_particles.E(i), mom_3);
				real weight = neighbor_params.kernel.template Weight<d>(wr_pj, interp_radius, type) * l_particles.WM(l_nb);
				e_arr[i] += weight * mom_3;
			}
		}
	}
}

template<int d>
template<class T>
inline T ParticleALEFilm<d>::XSPH_Smoothing(const real alpha, Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type, std::function<bool(const int)> filter)const
{
	Array<T> smoothed_f; smoothed_f.resize(e_particles.Size());
	AuxFunc::Fill(smoothed_f, Zero<T>());
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		if ((filter && !filter(i))|| is_rim_confidence[i] > 0.) {
			smoothed_f[i] = f_arr[i];
		}
		else {
			Array<int> nbs = e_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), radius);
			T sum = Zero<T>();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				VectorD wr_pj = e_particles.X(j) - e_particles.X(i);
				VectorT sr_pj = PointSet<d>::Rotate_To_TPlane(wr_pj, e_particles.E(i));
				real w = kernel.Weight<d - 1>(sr_pj, radius, kernel_type);
				sum += e_particles.SA(j) * (f_arr[j] - f_arr[i]) * w;
			}
			// account for boundary too. Assuming that boundary carries zero value (true for curvature, velocity)
			Array<int> b_nbs = boundary.nbs_searcher->Find_Neighbors(e_particles.X(i), radius);
			for (int k = 0; k < b_nbs.size(); k++) {
				int j = b_nbs[k];
				VectorD wr_pj = boundary.X(j) - e_particles.X(i);
				VectorT sr_pj = PointSet<d>::Rotate_To_TPlane(wr_pj, e_particles.E(i));
				real w = kernel.Weight<d - 1>(sr_pj, radius, kernel_type);
				sum += boundary.SA(j) * (-f_arr[i]) * w;
			}
			smoothed_f[i] = f_arr[i] + alpha * sum;
		}
	}
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		f_arr[i] = smoothed_f[i];
	}
}

template<int d>
Vector<real, d> ParticleALEFilm<d>::XSPH_Smoothing_Tang_Norm(const real alpha_tang, const real alpha_norm, Array<Vector<real,d>>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type, std::function<bool(const int)> filter)const
{
	Array<VectorD> smoothed_f; smoothed_f.resize(e_particles.Size());
	AuxFunc::Fill(smoothed_f, VectorD::Zero());
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		if ((filter && !filter(i)) || is_rim_confidence[i] > 0.) {
			smoothed_f[i] = f_arr[i];
		}
		else {
			Array<int> nbs = e_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), radius);
			VectorD sum = VectorD::Zero();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				VectorD wr_pj = e_particles.X(j) - e_particles.X(i);
				VectorT sr_pj = PointSet<d>::Rotate_To_TPlane(wr_pj, e_particles.E(i));
				real w = kernel.Weight<d - 1>(sr_pj, radius, kernel_type);
				sum += e_particles.SA(j) * (f_arr[j] - f_arr[i]) * w;
			}
			// account for boundary too. Assuming that boundary carries zero value (true for curvature, velocity)
			Array<int> b_nbs = boundary.nbs_searcher->Find_Neighbors(e_particles.X(i), radius);
			for (int k = 0; k < b_nbs.size(); k++) {
				int j = b_nbs[k];
				VectorD wr_pj = boundary.X(j) - e_particles.X(i);
				VectorT sr_pj = PointSet<d>::Rotate_To_TPlane(wr_pj, e_particles.E(i));
				real w = kernel.Weight<d - 1>(sr_pj, radius, kernel_type);
				sum += boundary.SA(j) * (-f_arr[i]) * w;
			}
			VectorD normal = e_particles.Normal(i);
			VectorD sum_norm = AuxFunc::Component_Along(sum, normal);
			VectorD sum_tang = sum - sum_norm;
			smoothed_f[i] = f_arr[i] + alpha_tang * sum_tang + alpha_norm * sum_norm;
		}
	}
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		f_arr[i] = smoothed_f[i];
	}
}

template<int d>
template<class T>
inline T ParticleALEFilm<d>::XSPH_Smoothed_Value(const VectorD& pos, const MatrixD& frame, T& value, Array<T>& f_arr, const KernelSPH& kernel, const real radius, KernelType kernel_type)const
{
	Array<int> nbs = e_particles.nbs_searcher->Find_Neighbors(pos, radius);
	T sum = Zero<T>();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		VectorD wr_pj = e_particles.X(j) - pos;
		VectorT sr_pj = PointSet<d>::Rotate_To_TPlane(wr_pj, frame);
		real w = kernel.Weight<d - 1>(sr_pj, radius, kernel_type);
		sum += e_particles.SA(j) * (f_arr[j] - value) * w;
	}
	// account for boundary too. Assuming that boundary carries zero value (true for curvature, velocity)
	Array<int> b_nbs = boundary.nbs_searcher->Find_Neighbors(pos, radius);
	for (int k = 0; k < b_nbs.size(); k++) {
		int j = b_nbs[k];
		VectorD wr_pj = boundary.X(j) - pos;
		VectorT sr_pj = PointSet<d>::Rotate_To_TPlane(wr_pj, frame);
		real w = kernel.Weight<d - 1>(sr_pj, radius, kernel_type);
		sum += boundary.SA(j) * (-value) * w;
	}
	return sum;
}

template<int d>
template<class T>
inline T ParticleALEFilm<d>::Boundary_Adhesion(Array<T>& f_arr)const
{
	if (analytical_boundary.Available()) {
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
			if (Eulerian_Near_Boundary(i, neighbor_params.sph_radius)) {
				f_arr[i] *= abs(e_particles.phis[i]) / neighbor_params.sph_radius;
			}
		}
	}
}

#endif

