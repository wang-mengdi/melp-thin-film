//////////////////////////////////////////////////////////////////////////
// ParticleALEFilm
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __ParticleALEFoam_h__
#define __ParticleALEFoam_h__
#include "PointSet.h"
#include "Fluid3DSPH.h"
#include "ParticleALEFilm.h"
#include "RandomNumber.h"
#include "ArrayIO.h"
#include "ArrayFunc.h"
#include "Timer.h"
#include "PointSetFunc.h"

#include <numeric>      // std::iota
#include <iomanip>
#include <fstream>

template<int d>
struct SharingInfo {
	Typedef_VectorDii(d);
	Typedef_MatrixD(d);
	int region_idx;
	real signed_distance; //if positive -> inside, negative -> outside
	real sa_ratio; //how does their total surface area compare to my total surface area
	VectorD proj_pos;
	MatrixD proj_frame;
	SharingInfo(int idx, real distance, real ratio, VectorD pos, MatrixD frame)
	{
		region_idx = idx;
		signed_distance = distance;
		sa_ratio = ratio;
		proj_pos = pos;
		proj_frame = frame;
	}
};

template<int d>
class ParticleALEFoam
{
	Typedef_VectorDii(d);
	Typedef_MatrixD(d);
public:
	using VectorT = Vector<real, d - 1>;
	using MatrixT = Matrix<real, d - 1>;

	std::function<real(const VectorD)> temp_func = nullptr;

	//pointer to 3d solver
	Fluid3DSPH<d>* fluid_3d = nullptr;
	Array<std::shared_ptr<ParticleALEFilm<d>>> regions;
	Array<Array<Array<SharingInfo<d>>>> master_list; //the list of inter-belongings, e.g. a point in A is also in B, C, D, etc. 
	Array<real> region_vols;
	Array<real> region_sas;
	Array<VectorD> region_mins;
	Array<VectorD> region_maxs;
	Array<real> region_modulus;
	MatrixXd stresses;

	int num_detonations = 0.;
	int num_detonations_done = 0.;

	bool recompute_init_COM = false;

	//the global boundaries
	// the two boundaries should represent the same geometric thing
	LagrangeParticles<d> boundary; //boundary is Lagrange particle
	ImplicitShape<d> analytical_boundary;
	
	//params for boundary
	real bottom = 0.;
	real top = 0.;
	real left = 0.;
	real right = 0.;
	real front = 0.;
	real back = 0.;

	//simulation info
	std::string output_dir;
	int curr_frame = 0;
	real cfl;
	real frame_rate;
	int verbose = 0;
	real e_dx; //the characteristic particle gap
	real sph_radius; //the characteristic kernel radius
	KernelSPH kernel;
	int simulation_fineness;
	real simulation_scale;
	real outer_pressure = 10.;
	bool e_only = false;
	real decay_strength = 0.;
	bool e_feels_gravity = false;
	real repulsion_multiplier = 1.;

	VectorD g = -1. * VectorD::Unit(1) * simulation_scale;

	VectorD init_COM = VectorD::Zero();
	VectorD curr_COM = VectorD::Zero();
	VectorD curr_COM_velocity = VectorD::Zero();
	Array<int> conserve_momentum_along;

	real max_vel = std::numeric_limits<real>::max();

	std::function<Array<int>(const int)> bubble_defreeze_func = nullptr;
	std::function<void(const int)> bubble_merge_func = nullptr;

	Array<int> freeze;

	void Initialize(const real _cfl, const real _frame_rate, std::string _output_dir)
	{
		double begin_time, end_time;
		if (verbose) {
			std::cout << "[Initialization]: foam init begin" << std::endl;
			begin_time = omp_get_wtime();
		}
		cfl = _cfl;
		frame_rate = _frame_rate;
		output_dir = _output_dir;

		boundary.Initialize();

		Reset_Master_List();
		stresses.resize(Num_Regions(), Num_Regions());
		stresses.setZero();
		Defreeze_Bubbles();
		Update_All_Geometry(0.);

		init_COM = curr_COM;
		std::cout << "[Initialization]:  Init COM: \n" << std::setprecision(5) << init_COM << std::endl;

		region_modulus.resize(Num_Regions());
		for (int i = 0; i < Num_Regions(); i++) {
			region_modulus[i] = 0.; //default to no resisistance. Merge directly
		}

		if (verbose) {
			end_time = omp_get_wtime();
			std::cout << "[Initialization]: foam init done. Time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
		}
	}

	// determine whether it is worth checking the particle sharing between 2 bubbles
	bool Has_Overlap(const int i, const int j)
	{
		//if (fluid_this.Is_Trivial() || fluid_that.Is_Trivial()) return false; // if either one is trivial, then don't bother checking
		if (Is_Region_Trivial(i) || Is_Region_Trivial(j)) return false; // if either one is trivial, then don't bother checking
		return (region_mins[i](0) < region_maxs[j](0) && region_mins[j](0) < region_maxs[i](0)
			&& region_mins[i](1) < region_maxs[j](1) && region_mins[j](1) < region_maxs[i](1)
			&& region_mins[i](2) < region_maxs[j](2) && region_mins[j](2) < region_maxs[i](2));
	}

	// update max velocity
	void Update_Max_Velocity(void)
	{
		// update max velocity of all regions
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			fluid.Update_Max_Velocity();
		}
		// consider all regions
		max_vel = 0.;
		for (int i = 0; i < regions.size(); i++) {
			max_vel = std::max(max_vel, regions[i]->max_vel);
		}
		if (fluid_3d != nullptr) {
			fluid_3d->Update_Max_Velocity();
			max_vel = std::max(max_vel, fluid_3d->max_vel);
		}
	}

	bool Is_Region_Trivial(int i) {
		return (regions[i]->Is_Trivial() || freeze[i]);
	}

	// update max velocity
	bool Is_Trivial(void)
	{
		for (int i = 0; i < regions.size(); i++) {
			if (!Is_Region_Trivial(i)) return false;
		}
		return true;
	}

	// number of regions (bubbles)
	int Num_Regions(void)
	{
		return regions.size();
	}

	// re-catalog the particle sharing
	void Update_Particle_Sharing_Info(void)
	{
		for (int i = 0; i < Num_Regions(); i++) {
			if (Is_Region_Trivial(i)) continue;
			ParticleALEFilm<d>& fluid = *regions[i];
			//AuxFunc::Fill(fluid.temp_vector, VectorD::Zero());
#pragma omp parallel for
			for (int m = 0; m < fluid.e_particles.Size(); m++) {
				master_list[i][m].clear(); //clear the guest list
				real total_sa_mine = fluid.e_particles.World_Sum_Value(fluid.e_particles.X(m), fluid.e_particles.SARef(), kernel, sph_radius, KernelType::QUINTIC);
				for (int j = 0; j < Num_Regions(); j++) {
					if (j == i) continue;
					ParticleALEFilm<d>& other_fluid = *regions[j];
					if (Is_Region_Trivial(j)) continue;
					if (!Has_Overlap(i, j)) continue;
					int other_nearest_nb = other_fluid.e_particles.nbs_searcher->Find_Nearest_Nb(fluid.e_particles.X(m));
					//int num_other_neighbors = other_fluid.e_particles.nbs_searcher->Find_Neighbors(fluid.e_particles.X(m), fluid.neighbor_params.sph_radius).size();
					bool has_nbs = true;
					VectorD other_avg_pos = other_fluid.e_particles.World_Avg_Value_KNN(fluid.e_particles.X(m), other_fluid.e_particles.XRef(), other_fluid.neighbor_params.kernel, other_fluid.neighbor_params.interp_radius, 10, KernelType::QUINTIC, has_nbs);
					if (!has_nbs) continue;
					//VectorD proj_pos = other_fluid.e_particles.X(other_nearest_nb);
					VectorD proj_pos = other_avg_pos;
					MatrixD proj_frame;
					////int not_enough_neighbors = other_fluid.e_particles.surface.Nearest_Geometry(proj_pos, proj_frame, 7, false);
					int not_enough_neighbors = other_fluid.e_particles.surface.Nearest_Geometry_KNN(proj_pos, proj_frame, 20, 7, false);
					//// first make sure that the projection is successful
					//// if can't find enough neighbors in B, then it definitely does not belong to B
					if (not_enough_neighbors > 0) continue;
					VectorD diff = proj_pos - fluid.e_particles.X(m);
					VectorD normal = proj_frame.col(2);
					int sign = diff.normalized().dot(normal) > (real)0 ? 1 : -1; // 1 means inside 2 means outside
					//if (sign > 0 || diff.norm() < sph_radius) { // if is INSIDE or is NEAR
					real signed_distance = sign * diff.norm();
					//if (i == 1 && m == 367) {
					//	std::cout << "FROM INFO!" << std::endl;
					//	std::cout << "now checking bubble: " << j << std::endl;
					//	std::cout << "Projected pos: " << proj_pos << std::endl;
					//	std::cout << "Diff: " << diff << std::endl;
					//	std::cout << "Diff normalized: " << diff.normalized() << std::endl;
					//	std::cout << "normal: " << normal << std::endl;
					//	std::cout << "sign is: " << sign << std::endl;
					//	std::cout << "signed dist is: " << signed_distance << std::endl;
					//	//fluid.temp_vector[m] = diff;
					//}
					if (signed_distance > -sph_radius) { // if is INSIDE or is NEAR (note positive signed distance means inside(BAD), negative means outside) (the more negative means farther away)
					//if (diff.norm() < sph_radius) { // if is INSIDE or is NEAR
						real total_sa_theirs = other_fluid.e_particles.World_Sum_Value(fluid.e_particles.X(m), other_fluid.e_particles.SARef(), kernel, sph_radius, KernelType::QUINTIC);
						real sa_ratio = total_sa_theirs / total_sa_mine;
						master_list[i][m].push_back(SharingInfo<d>(j, signed_distance, sa_ratio, proj_pos, proj_frame));
					}
				}
				if (master_list[i][m].size() >= 1) {
					fluid.temp_scalar[m] = 3.e-7;
				}
				else {
					fluid.temp_scalar[m] = 5.e-7;
				}
			}
		}
	}

	//resize the master list appropriatelyf
	//needs to be called when number of particles change
	void Reset_Master_List(void)
	{
		freeze.resize(Num_Regions());
		AuxFunc::Fill(freeze, 0);
		master_list.resize(Num_Regions());
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			Array<Array<SharingInfo<d>>> region_list;
			region_list.resize(fluid.e_particles.Size());
#pragma omp parallel for
			for (int m = 0; m < fluid.e_particles.Size(); m++) {
				Array<SharingInfo<d>> list;
				region_list[m] = list;
			}
			master_list[i] = region_list;
		}
	}

	// this method computes the curvature near a plateau boarder
	// it applies to those particles that belong to multiple bubbles
	// and yet is not inside any other bubble (if penetration happens, then we use the original curvature instead)
	// or is not an external particle (is face-to-face with other bubbles)
	// when computing the laplacian, we will also take into account all the shared bubbles
	// but we shall also exclude the particles that are not external, or is inside other particles
	VectorD Compute_Plateau_Curvature(const int i, const int m, const VectorD& _external_normal) {
		if (_external_normal.norm() <= 1.e-8) return VectorD::Zero();
		VectorD external_normal = _external_normal / _external_normal.norm();
		real filter_out_distance = 2. * e_dx;
		Array<SharingInfo<d>>& also_in = master_list[i][m];
		//first consider the lap in own bubble
		ParticleALEFilm<d>& fluid = *regions[i];
		real lap = 0.;
		//first find the external-normal vector of the particle
		Array<int> nbs = fluid.e_particles.nbs_searcher->Find_Neighbors(fluid.e_particles.X(m), fluid.neighbor_params.sph_radius);
		nbs = fluid.e_particles.nbs_searcher->Find_Neighbors(fluid.e_particles.X(m), fluid.neighbor_params.sph_radius);
		for (int k = 0; k < nbs.size(); k++) {
			int l = nbs[k]; // l is the neighbor idx
			if (m == l) continue;
			Array<SharingInfo<d>>& also_in_nb = master_list[i][l];
			bool is_external = true;
			for (int q = 0; q < also_in_nb.size(); q++) {
				if (also_in_nb[q].signed_distance > -filter_out_distance) {
					is_external = false;
					break;
				}
			}
			if (!is_external) continue;
			real S_j = fluid.e_particles.SA(l);
			VectorT lr_ji = fluid.e_particles.Surface_Relative_Vec(l, m);
			real norm_ji = std::max(lr_ji.norm(), fluid.e_particles.surface.t_dx * 0.01);
			VectorT grad_W = fluid.neighbor_params.kernel.template Grad<d - 1>(lr_ji, fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
			real f_j = fluid.e_particles.X(l).dot(external_normal);
			real f_i = fluid.e_particles.X(m).dot(external_normal);
			lap += S_j * (f_j - f_i) * 2 * grad_W.norm() / norm_ji;
		}
		//then consider the lap in all shared bubbles
		for (int k = 0; k < also_in.size(); k++) {
			int j = also_in[k].region_idx;
			ParticleALEFilm<d>& other_fluid = *regions[j];
			Array<int> nbs = other_fluid.e_particles.nbs_searcher->Find_Neighbors(fluid.e_particles.X(m), fluid.neighbor_params.sph_radius);
			for (int k = 0; k < nbs.size(); k++) {
				int l = nbs[k]; // l is the neighbor idx
				Array<SharingInfo<d>>& also_in_nb = master_list[j][l];
				bool is_external = true;
				for (int q = 0; q < also_in_nb.size(); q++) {
					if (also_in_nb[q].signed_distance > -filter_out_distance) {
						is_external = false;
						break;
					}
				}
				if (!is_external) continue;
				real S_j = other_fluid.e_particles.SA(l);
				VectorT lr_ji = PointSet<d>::Rotate_To_TPlane(fluid.e_particles.X(m) - other_fluid.e_particles.X(l), fluid.e_particles.E(m));
				//VectorT lr_ji = fluid.e_particles.Surface_Relative_Vec(l, m);
				real norm_ji = std::max(lr_ji.norm(), fluid.e_particles.surface.t_dx * 0.01);
				VectorT grad_W = fluid.neighbor_params.kernel.template Grad<d - 1>(lr_ji, fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
				real f_j = other_fluid.e_particles.X(l).dot(external_normal);
				real f_i = fluid.e_particles.X(m).dot(external_normal);
				lap += S_j * (f_j - f_i) * 2 * grad_W.norm() / norm_ji;
			}
		}
		return lap * _external_normal;
	}

	void Update_Particle_Sharing_Momentum(const real dt)
	{
		// first record the OLD velocities, to ensure that this doesn't depend on order of bubbles
		Array<Array<VectorD>> old_velocities;
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			Array<VectorD> old_velocities_i;
			old_velocities_i.resize(fluid.e_particles.Size());
#pragma omp parallel for
			for (int m = 0; m < fluid.e_particles.Size(); m++) {
				old_velocities_i[m] = fluid.e_particles.V(m);
			}
			old_velocities.push_back(old_velocities_i);
		}
		// done recording OLD velocities

		// define some parameters here
		real close_enough_threshold = 0.1 * simulation_scale;
		real adhesion_param = 0. * 20.;
		real repulsion_param = 3.3e2 * repulsion_multiplier; //3.3e3
		real viscosity = 0.66; //should be between 0 and 1;

		//for each bubble
		for (int i = 0; i < Num_Regions(); i++) {
			if (Is_Region_Trivial(i)) continue;
			ParticleALEFilm<d>& fluid = *regions[i];
			//AuxFunc::Fill(fluid.temp_vector, VectorD::Zero());
			//std::cout << "for i: " << i << " total velocity before Momentum exchange is: " << AuxFunc::Sum(fluid.e_particles.VRef()) << std::endl;
			Array<VectorD> repulsive_adhesive_accs; repulsive_adhesive_accs.resize(fluid.e_particles.Size());
			Array<Array<real>> recorded_stresses; recorded_stresses.resize(fluid.e_particles.Size());
#pragma omp parallel for
			for (int m = 0; m < fluid.e_particles.Size(); m++) {
				Array<real> tmp(Num_Regions(), 0.);
				recorded_stresses[m] = tmp;
			}
#pragma omp parallel for
			for (int m = 0; m < fluid.e_particles.Size(); m++) {
				//VectorD repulsive_adhesive_acc = VectorD::Zero();
				//VectorD total_norm_vel = AuxFunc::Component_Along(fluid.e_particles.V(m), fluid.e_particles.Normal(m));
				//VectorD tang_vel = fluid.e_particles.V(m) - total_norm_vel;
				//Array<SharingInfo<d>>& also_in = master_list[i][m];
				//int num_close_to = 0;
				//for (int k = 0; k < also_in.size(); k++) {
				//	int j = also_in[k].region_idx;
				//	ParticleALEFilm<d>& other_fluid = *regions[j];
				//	VectorD proj_pos = also_in[k].proj_pos;
				//	MatrixD proj_frame = also_in[k].proj_frame;
				//	VectorD proj_normal = proj_frame.col(2);
				//	real signed_dist = also_in[k].signed_distance;
				//	VectorD diff = proj_pos - fluid.e_particles.X(m);
				//	if (diff.norm() > 1. * e_dx) continue;
				//	num_close_to++;
				//	// the repulsion and adhesion are emperical and should be tuned
				//	real avg_sa = region_sas[i] / fluid.e_particles.Size();
				//	if (signed_dist < 0) {//if is outside, then apply adhesion
				//		real adhesion = 20 * diff.norm();
				//		repulsive_adhesive_acc += adhesion * diff.normalized();
				//		//recorded_stresses[m][j] += avg_sa + adhesion;
				//	}
				//	else if (signed_dist > 0) {// if inside, then apply repulsion
				//		real repulsion = 1.e5 * MathFunc::Quick_Pow(diff.norm(), 2);
				//		repulsive_adhesive_acc += repulsion * diff.normalized();
				//		recorded_stresses[m][j] += avg_sa + repulsion;
				//	}
				//	VectorD collected_vel = other_fluid.e_particles.Surface_Convolution_Value(proj_pos, proj_frame, old_velocities[j], other_fluid.neighbor_params.kernel, other_fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
				//	total_norm_vel += AuxFunc::Component_Along(collected_vel, fluid.e_particles.Normal(m));
				//}
				//total_norm_vel /= (num_close_to + 1);
				//total_norm_vel = 0.6 * total_norm_vel + 0.4 * AuxFunc::Component_Along(fluid.e_particles.V(m), fluid.e_particles.Normal(m));
				//fluid.e_particles.V(m) = total_norm_vel + tang_vel;
				//fluid.e_particles.V(m) += dt * repulsive_adhesive_acc;

				VectorD repulsive_adhesive_acc = VectorD::Zero();
				VectorD total_norm_vel = AuxFunc::Component_Along(fluid.e_particles.V(m), fluid.e_particles.Normal(m));
				VectorD total_tang_vel = fluid.e_particles.V(m) - total_norm_vel;
				VectorD tang_vel = fluid.e_particles.V(m) - total_norm_vel;
				Array<SharingInfo<d>>& also_in = master_list[i][m];
				int num_close_to = 0;
				//compute total sa_ratio
				real total_sa_ratio = 1.;
				for (int k = 0; k < also_in.size(); k++) {
					int j = also_in[k].region_idx;
					ParticleALEFilm<d>& other_fluid = *regions[j];
					VectorD proj_pos = also_in[k].proj_pos;
					MatrixD proj_frame = also_in[k].proj_frame;
					VectorD proj_normal = proj_frame.col(2);
					real signed_dist = also_in[k].signed_distance;
					VectorD diff = proj_pos - fluid.e_particles.X(m);
					if (diff.norm() > close_enough_threshold) continue;
					num_close_to++;
					// the repulsion and adhesion are emperical and should be tuned
					real avg_sa = region_sas[i] / fluid.e_particles.Size();
					if (signed_dist < 0) {//if is outside, then apply adhesion
						if ((stresses(i, j) + stresses(j, i)) / 2. > (region_modulus[i] + region_modulus[j]) / 2.) {
						//if (true) {
							real adhesion = adhesion_param * diff.norm();
							repulsive_adhesive_acc += adhesion * diff.normalized();
						}
					}
					else if (signed_dist > 0) {// if inside, then apply repulsion
						real repulsion = repulsion_param * diff.norm();
						repulsive_adhesive_acc += repulsion * diff.normalized();
						//recorded_stresses[m][j] += avg_sa  repulsion;
						recorded_stresses[m][j] += avg_sa * repulsion;
					}

					VectorD other_smoothed_V = other_fluid.XSPH_Smoothed_Value(proj_pos, proj_frame, old_velocities[i][m], old_velocities[j], other_fluid.neighbor_params.kernel, other_fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
					VectorD collected_vel = fluid.e_particles.V(m) + viscosity * other_fluid.XSPH_Smoothed_Value(proj_pos, proj_frame, old_velocities[i][m], old_velocities[j], other_fluid.neighbor_params.kernel, other_fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
					VectorD collected_vel_norm = AuxFunc::Component_Along(collected_vel, fluid.e_particles.Normal(m));
					VectorD collected_vel_tang = collected_vel - collected_vel_norm;
					/*if (fluid.e_particles.Normal(m).dot(collected_vel_norm) < 0.) {*/
					if (true) {
						total_sa_ratio += also_in[k].sa_ratio;
						total_norm_vel += also_in[k].sa_ratio * collected_vel_norm;
						total_tang_vel += also_in[k].sa_ratio * collected_vel_tang;
					}
				}

				if (repulsive_adhesive_acc.dot(fluid.e_particles.Normal(m)) < 0.) {
					repulsive_adhesive_accs[m] = repulsive_adhesive_acc;
				}
				else {
					repulsive_adhesive_accs[m] = VectorD::Zero();
				}

				total_norm_vel /= total_sa_ratio;
				total_tang_vel /= total_sa_ratio;

				fluid.e_particles.V(m) = total_norm_vel + total_tang_vel;
				fluid.e_particles.V(m) += dt * repulsive_adhesive_accs[m];
			}
			for (int m = 0; m < fluid.e_particles.Size(); m++) {
				for (int j = 0; j < Num_Regions(); j++) {
					stresses(i, j) += recorded_stresses[m][j];
				}
			}
		}

//		// define some parameters here
//		real close_enough_threshold = 0.1 * simulation_scale;
//		real adhesion_param = 3. * 20.;
//		real repulsion_param = 3.3e2 * repulsion_multiplier; //3.3e3
//		real viscosity = 0.999; //should be between 0 and 1;
//
//		//for each bubble
//		for (int i = 0; i < Num_Regions(); i++) {
//			if (Is_Region_Trivial(i)) continue;
//			ParticleALEFilm<d>& fluid = *regions[i];
//			//AuxFunc::Fill(fluid.temp_vector, VectorD::Zero());
//			//std::cout << "for i: " << i << " total velocity before Momentum exchange is: " << AuxFunc::Sum(fluid.e_particles.VRef()) << std::endl;
//			Array<VectorD> repulsive_adhesive_accs; repulsive_adhesive_accs.resize(fluid.e_particles.Size());
//			Array<Array<real>> recorded_stresses; recorded_stresses.resize(fluid.e_particles.Size());
//#pragma omp parallel for
//			for (int m = 0; m < fluid.e_particles.Size(); m++) {
//				Array<real> tmp(Num_Regions(), 0.);
//				recorded_stresses[m] = tmp;
//			}
//#pragma omp parallel for
//			for (int m = 0; m < fluid.e_particles.Size(); m++) {
//				VectorD repulsive_adhesive_acc = VectorD::Zero();
//				VectorD total_norm_vel = AuxFunc::Component_Along(fluid.e_particles.V(m), fluid.e_particles.Normal(m));
//				VectorD total_tang_vel = fluid.e_particles.V(m) - total_norm_vel;
//				VectorD tang_vel = fluid.e_particles.V(m) - total_norm_vel;
//				Array<SharingInfo<d>>& also_in = master_list[i][m];
//				int num_close_to = 0;
//				//compute total sa_ratio
//				real total_sa_ratio = 1.;
//				for (int k = 0; k < also_in.size(); k++) {
//					int j = also_in[k].region_idx;
//					ParticleALEFilm<d>& other_fluid = *regions[j];
//					VectorD proj_pos = also_in[k].proj_pos;
//					MatrixD proj_frame = also_in[k].proj_frame;
//					VectorD proj_normal = proj_frame.col(2);
//					real signed_dist = also_in[k].signed_distance;
//					VectorD diff = proj_pos - fluid.e_particles.X(m);
//					if (diff.norm() > close_enough_threshold) continue;
//					num_close_to++;
//					// the repulsion and adhesion are emperical and should be tuned
//					real avg_sa = region_sas[i] / fluid.e_particles.Size();
//					if (signed_dist < 0) {//if is outside, then apply adhesion
//						//if ((stresses(i, j) + stresses(j, i)) / 2. > (region_modulus[i] + region_modulus[j]) / 2.) {
//						if (true) {
//							real adhesion = adhesion_param * diff.norm();
//							repulsive_adhesive_acc += adhesion * diff.normalized();
//						}
//						//real adhesion = adhesion_param * diff.norm();
//						//repulsive_adhesive_acc += adhesion * diff.normalized();
//						//recorded_stresses[m][j] += avg_sa + adhesion;
//					}
//					else if (signed_dist > 0) {// if inside, then apply repulsion
//						real repulsion = repulsion_param * diff.norm();
//						repulsive_adhesive_acc += repulsion * diff.normalized();
//						//recorded_stresses[m][j] += avg_sa  repulsion;
//						recorded_stresses[m][j] += avg_sa * repulsion;
//					}
//					//if (i == 1 && m == 274) {
//					//	std::cout << "checking i == 1, m == 274" << std::endl;
//					//	std::cout << "j is: " << j << std::endl;
//					//	std::cout << "signed_dist: " << signed_dist << std::endl;
//					//}
//					VectorD other_smoothed_V = other_fluid.XSPH_Smoothed_Value(proj_pos, proj_frame, old_velocities[i][m], old_velocities[j], other_fluid.neighbor_params.kernel, other_fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
//					VectorD collected_vel = fluid.e_particles.V(m) + viscosity * other_fluid.XSPH_Smoothed_Value(proj_pos, proj_frame, old_velocities[i][m], old_velocities[j], other_fluid.neighbor_params.kernel, other_fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
//					VectorD collected_vel_norm = AuxFunc::Component_Along(collected_vel, fluid.e_particles.Normal(m));
//					VectorD collected_vel_tang = collected_vel - collected_vel_norm;
//					/*if (fluid.e_particles.Normal(m).dot(collected_vel_norm) < 0.) {*/
//					if (true){
//						total_sa_ratio += also_in[k].sa_ratio;
//						total_norm_vel += also_in[k].sa_ratio * collected_vel_norm;
//						total_tang_vel += also_in[k].sa_ratio * collected_vel_tang;
//					}
//				}
//
//				if (repulsive_adhesive_acc.dot(fluid.e_particles.Normal(m)) < 0.) {
//					repulsive_adhesive_accs[m] = repulsive_adhesive_acc;
//				}
//				else {
//					repulsive_adhesive_accs[m] = VectorD::Zero();
//				}
//
//				total_norm_vel /= total_sa_ratio;
//				total_tang_vel /= total_sa_ratio;
//				//total_norm_vel = viscosity * total_norm_vel + (1. - viscosity) * AuxFunc::Component_Along(fluid.e_particles.V(m), fluid.e_particles.Normal(m));
//				//fluid.e_particles.V(m) = total_norm_vel + total_tang_vel;
//				
//				//fluid.e_particles.V(m) = total_norm_vel + (0.3 * tang_vel + 0.7 * total_tang_vel);
//				
//				fluid.e_particles.V(m) = total_norm_vel + tang_vel;
//				
//				//fluid.temp_vector[m] = fluid.e_particles.V(m) - fluid.temp_vector[m];
//				//fluid.e_particles.V(m) += dt * repulsive_adhesive_acc;
//				//fluid.temp_vector[m] = repulsive_adhesive_acc * dt;
//				//fluid.temp_vector[m] = fluid.e_particles.V(m) - fluid.temp_vector[m];
//				//fluid.temp_vector[m] = repulsive_adhesive_accs[m] * dt;
//			}
//			fluid.XSPH_Smoothing(0.999, repulsive_adhesive_accs, fluid.neighbor_params.kernel, fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
//
//#pragma omp parallel for
//			for (int m = 0; m < fluid.e_particles.Size(); m++) {
//				fluid.e_particles.V(m) += dt * AuxFunc::Component_Along(repulsive_adhesive_accs[m], fluid.e_particles.Normal(m));
//			}
//			for (int m = 0; m < fluid.e_particles.Size(); m++) {
//				for (int j = 0; j < Num_Regions(); j++) {
//					stresses(i, j) += recorded_stresses[m][j];
//				}
//			}
//			std::function<bool(const int)> near_junction = [&](const int idx)->real {return (master_list[i][idx].size() > 0); };
//			fluid.XSPH_Smoothing_Tang_Norm(0., 0.999, fluid.e_particles.VRef(), fluid.neighbor_params.kernel, fluid.neighbor_params.sph_radius, KernelType::QUINTIC, near_junction);
//			//std::cout << "for i: " << i << " total velocity after Momentum exchange is: " << AuxFunc::Sum(fluid.e_particles.VRef()) << std::endl;
//		}
		//std::cout << "Stress record: \n" << stresses << std::endl;
		for (int i = 0; i < Num_Regions(); i++) {
			if (Is_Region_Trivial(i)) continue;
			ParticleALEFilm<d>& fluid = *regions[i];
			for (int m = 0; m < fluid.e_particles.Size(); m++) {
				Array<SharingInfo<d>> new_also_in;
				Array<SharingInfo<d>>& also_in = master_list[i][m];
				for (int k = 0; k < also_in.size(); k++) {
					if ((stresses(i, also_in[k].region_idx) + stresses(also_in[k].region_idx, i))/2. > (region_modulus[i] + region_modulus[also_in[k].region_idx])/2.) {
						new_also_in.push_back(also_in[k]);
					}
				}
				also_in = new_also_in;
			}
		}
		//std::cout << "done updating momentum" << std::endl;
	}


//void Update_Particle_Sharing_Momentum(const real dt)
//{
//	// first record the OLD velocities, to ensure that this doesn't depend on order of bubbles
//	Array<Array<VectorD>> old_velocities;
//	for (int i = 0; i < Num_Regions(); i++) {
//		ParticleALEFilm<d>& fluid = *regions[i];
//		Array<VectorD> old_velocities_i;
//		old_velocities_i.resize(fluid.e_particles.Size());
//#pragma omp parallel for
//		for (int m = 0; m < fluid.e_particles.Size(); m++) {
//			old_velocities_i[m] = fluid.e_particles.V(m);
//		}
//		old_velocities.push_back(old_velocities_i);
//	}
//	// done recording OLD velocities
//
//	//for each bubble
//	for (int i = 0; i < Num_Regions(); i++) {
//		ParticleALEFilm<d>& fluid = *regions[i];
//
//		Array<Array<real>> recorded_stresses; recorded_stresses.resize(fluid.e_particles.Size());
//#pragma omp parallel for
//		for (int m = 0; m < fluid.e_particles.Size(); m++) {
//			Array<real> tmp(Num_Regions(), 0.);
//			recorded_stresses[m] = tmp;
//		}
//#pragma omp parallel for
//		for (int m = 0; m < fluid.e_particles.Size(); m++) {
//			VectorD repulsive_adhesive_acc = VectorD::Zero();
//			VectorD total_norm_vel = AuxFunc::Component_Along(fluid.e_particles.V(m), fluid.e_particles.Normal(m));
//			VectorD tang_vel = fluid.e_particles.V(m) - total_norm_vel;
//			Array<SharingInfo<d>>& also_in = master_list[i][m];
//			int num_close_to = 0;
//			for (int k = 0; k < also_in.size(); k++) {
//				int j = also_in[k].region_idx;
//				ParticleALEFilm<d>& other_fluid = *regions[j];
//				VectorD proj_pos = also_in[k].proj_pos;
//				MatrixD proj_frame = also_in[k].proj_frame;
//				VectorD proj_normal = proj_frame.col(2);
//				real signed_dist = also_in[k].signed_distance;
//				VectorD diff = proj_pos - fluid.e_particles.X(m);
//				if (diff.norm() > 1. * e_dx) continue;
//				num_close_to++;
//				// the repulsion and adhesion are emperical and should be tuned
//				real avg_sa = region_sas[i] / fluid.e_particles.Size();
//				if (signed_dist < 0) {//if is outside, then apply adhesion
//					real adhesion = 20 * diff.norm();
//					repulsive_adhesive_acc += adhesion * diff.normalized();
//					//recorded_stresses[m][j] += avg_sa + adhesion;
//				}
//				else if (signed_dist > 0) {// if inside, then apply repulsion
//					real repulsion = 1.e5 * MathFunc::Quick_Pow(diff.norm(), 2);
//					repulsive_adhesive_acc += repulsion * diff.normalized();
//					recorded_stresses[m][j] += avg_sa + repulsion;
//				}
//				VectorD collected_vel = other_fluid.e_particles.Surface_Convolution_Value(proj_pos, proj_frame, old_velocities[j], other_fluid.neighbor_params.kernel, other_fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
//				total_norm_vel += AuxFunc::Component_Along(collected_vel, fluid.e_particles.Normal(m));
//			}
//			total_norm_vel /= (num_close_to + 1);
//			total_norm_vel = 0.6 * total_norm_vel + 0.4 * AuxFunc::Component_Along(fluid.e_particles.V(m), fluid.e_particles.Normal(m));
//			fluid.e_particles.V(m) = total_norm_vel + tang_vel;
//			fluid.e_particles.V(m) += dt * repulsive_adhesive_acc;
//		}
//		for (int m = 0; m < fluid.e_particles.Size(); m++) {
//			for (int j = 0; j < Num_Regions(); j++) {
//				stresses(i, j) += recorded_stresses[m][j];
//			}
//		}
//	}
//	std::cout << "Stress record: \n" << stresses << std::endl;
//	for (int i = 0; i < Num_Regions(); i++) {
//		ParticleALEFilm<d>& fluid = *regions[i];
//		for (int m = 0; m < fluid.e_particles.Size(); m++) {
//			Array<SharingInfo<d>> new_also_in;
//			Array<SharingInfo<d>>& also_in = master_list[i][m];
//			for (int k = 0; k < also_in.size(); k++) {
//				if ((stresses(i, also_in[k].region_idx) + stresses(also_in[k].region_idx, i)) / 2. > (region_modulus[i] + region_modulus[also_in[k].region_idx]) / 2.) {
//					new_also_in.push_back(also_in[k]);
//				}
//			}
//			also_in = new_also_in;
//		}
//	}
//	std::cout << "done updating momentum" << std::endl;
//}


	void Interchange_L_Particles(const real dt)
	{
		// next owner initialize
		Array<Array<int>> next_owner;
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			Array<int> owner;
			owner.resize(fluid.l_particles.Size());
			AuxFunc::Fill(owner, i);
			next_owner.push_back(owner);
		}
		// next owner initialized

		//for each bubble, for each l particles, decide next owner based on the Gamma --- surfactant concentration
		//if current owner has much higher surfactant concentration than the next owner
		//there is a larger chance that it will be tranferred over to that new owner
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			Array<real> probs; probs.resize(fluid.e_particles.Size()); AuxFunc::Fill(probs, 0.);
			Array<int> nexts; nexts.resize(fluid.e_particles.Size()); AuxFunc::Fill(nexts, -1);
#pragma omp parallel for
			for (int m = 0; m < fluid.e_particles.Size(); m++) {
				Array<SharingInfo<d>>& also_in = master_list[i][m];
				int most_likely_region = i;
				real most_likely_prob = 0.;
				for (int k = 0; k < also_in.size(); k++) { //loop over all other particles
					real percentage = 0.;
					// first compute percentage due to surfactant concentration
					int j = also_in[k].region_idx;
					ParticleALEFilm<d>& other_fluid = *regions[j];
					int their_nearest_nb = other_fluid.e_particles.nbs_searcher->Find_Nearest_Nb(fluid.e_particles.X(m));
					real my_gamma = fluid.e_particles.H(m);
					real their_gamma = other_fluid.e_particles.H(their_nearest_nb);
					real ratio = 1.;
					if (my_gamma > 0.) { ratio = their_gamma / my_gamma; }
					ratio = std::min<real>(ratio, 1.);
					real strength = 0.25;
					real soap_percentage = 0.5 * strength * (1. - ratio); //if strength = 1, and if I have 1 unit, the other bubble have 0 unit, I will give half of mine to it.
					soap_percentage *= also_in[k].sa_ratio;
					// second compute percentage due to velocity
					//real velocity_percentage = 0.;

					//// Compute distance before and after
					////VectorD proj_pos = fluid.e_particles.X(m);
					////MatrixD proj_frame;
					////other_fluid.e_particles.surface.Nearest_Geometry_KNN(proj_pos, proj_frame, 10, 7, false);
					////real distance_before = (proj_pos - fluid.e_particles.X(m)).norm();
					////VectorD proj_pos2 = fluid.e_particles.X(m) + 5 * dt * fluid.e_particles.LV(m);
					////MatrixD proj_frame2;
					////other_fluid.e_particles.surface.Nearest_Geometry_KNN(proj_pos2, proj_frame2, 10, 7, false);
					////real distance_after = (proj_pos2 - fluid.e_particles.X(m)).norm();
					//// 
					//real distance_before = (other_fluid.e_particles.X(their_nearest_nb) - fluid.e_particles.X(m)).norm();
					//real distance_after = (other_fluid.e_particles.X(their_nearest_nb) - (fluid.e_particles.X(m) + dt * fluid.e_particles.LV(m))).norm();
					////
					//real distance_ratio = 0.;
					//if (distance_before > 0.01 * e_dx) {
					//	distance_ratio = 1. - std::min<real>(distance_after / distance_before, 1.);
					//} // a score between 0-1, 0 means stay in the same bubble, 1 means go to the other bubble
					//// add a criterion for alignednes
					////VectorD my_normal = fluid.e_particles.Normal(m);
					////VectorD their_normal = other_fluid.e_particles.Normal(their_nearest_nb);	
					//VectorD my_normal = fluid.plateau_normals[m];
					//VectorD their_normal = other_fluid.plateau_normals[their_nearest_nb];
					//real gauge = (my_normal + their_normal).norm() / 2;
					////if ((my_normal + their_normal).norm() < 0.1) { velocity_percentage = 0.; }
					////else {
					////	real alignedness = my_normal.dot(their_normal); // the bigger (closer to 1) the better (easier to transfer)
					////	if (alignedness < 0.) distance_ratio *= 0.; //1. + alignedness; // when compleltely opposite, then don't transfer
					////	velocity_percentage = distance_ratio;
					////}
					//distance_ratio *= gauge;
					//real max_transport = dt * 3;
					////std::cout << "Max Transport? " << max_transport << std::endl;
					//distance_ratio = std::min<real>(distance_ratio, max_transport);
					//velocity_percentage = distance_ratio;
					//// done
					//percentage = soap_percentage + velocity_percentage;
					percentage = soap_percentage;
					// multiply by sa_ratio
					//real weighted_percentage = std::min<real>(1., percentage * also_in[k].sa_ratio); //need to multiply by the confidence that this e particle is shared between the two regions
					if (percentage > most_likely_prob) {
						most_likely_prob = percentage;
						most_likely_region = j;
					}
				}
				nexts[m] = most_likely_region;
				probs[m] = most_likely_prob;
				//fluid.temp_vector[m] = probs[m] * fluid.e_particles.Normal(m);
			}
#pragma omp parallel for
			for (int m = 0; m < fluid.l_particles.Size(); m++) {
				int nearest_nb = fluid.e_particles.nbs_searcher->Find_Nearest_Nb(fluid.l_particles.X(m));
				int most_likely_region = nexts[nearest_nb];
				real most_likely_prob = probs[nearest_nb];
				//Array<SharingInfo<d>>& also_in = master_list[i][nearest_nb];
				//int most_likely_region = i;
				//real most_likely_prob = 0.;
				//for (int k = 0; k < also_in.size(); k++) { //loop over all other particles
				//	real percentage = 0.;
				//	// first compute percentage due to surfactant concentration
				//	int j = also_in[k].region_idx;
				//	ParticleALEFilm<d>& other_fluid = *regions[j];
				//	int their_nearest_nb = other_fluid.e_particles.nbs_searcher->Find_Nearest_Nb(fluid.l_particles.X(m));
				//	real my_gamma = fluid.e_particles.Gamma(nearest_nb);
				//	real their_gamma = other_fluid.e_particles.Gamma(their_nearest_nb);
				//	real ratio = 1.;
				//	if (my_gamma > 0.) { ratio = their_gamma / my_gamma; }
				//	ratio = std::min<real>(ratio, 1.);
				//	real strength = .0;
				//	real soap_percentage = 0.5 * strength * (1. - ratio); //if strength = 1, and if I have 1 unit, the other bubble have 0 unit, I will give half of mine to it.
				//	soap_percentage *= also_in[k].sa_ratio;
				//	// second compute percentage due to velocity
				//	real velocity_percentage = 0.;

				//	// Compute distance before and afert
				//	//VectorD proj_pos; 
				//	//other_fluid.e_particles.surface.Project_To_Surface(pos);
				//	//VectorD v = proj_pos - pos;
				//	//real distance_before v.norm();
				//	real distance_before = (other_fluid.e_particles.X(their_nearest_nb) - fluid.l_particles.X(m)).norm();
				//	real distance_after = (other_fluid.e_particles.X(their_nearest_nb) - (fluid.l_particles.X(m) + dt * fluid.l_particles.V(m))).norm();
				//	//if (distance_after < distance_before) {
				//	//VectorD there = other_fluid.e_particles.X(their_nearest_nb) - fluid.l_particles.X(m);
				//	//if (there.dot(fluid.l_particles.V(m)) > 0.) {
				//	//	VectorD my_normal = fluid.plateau_normals[nearest_nb];
				//	//	VectorD their_normal = other_fluid.plateau_normals[their_nearest_nb];
				//	//	//VectorD my_normal = fluid.e_particles.Normal(nearest_nb);
				//	//	//VectorD their_normal = other_fluid.e_particles.Normal(their_nearest_nb);			
				//	//	real alignedness = my_normal.dot(their_normal); // the bigger (closer to 1) the better (easier to transfer)
				//	//	//alignedness = std::max<real>(alignedness, 0.); // a score between 0-1. 0 means opposite, 1 means aligned
				//	//	alignedness = 0.5 * (alignedness+1.); // a score between 0-1. 0 means opposite, 1 means aligned
				//	//	velocity_percentage = alignedness * also_in[k].sa_ratio;
				//	//	//velocity_percentage = also_in[k].sa_ratio;
				//	//	//velocity_percentage *= 0.5;
				//	//}
				//	real distance_ratio = 0.;
				//	if (distance_before > 0.01 * e_dx) {
				//		distance_ratio = 1.-std::min<real>(distance_after / distance_before, 1.);
				//	} // a score between 0-1, 0 means stay in the same bubble, 1 means go to the other bubble
				//	// add a criterion for alignednes
				//	//VectorD my_normal = fluid.e_particles.Normal(nearest_nb);
				//	//VectorD their_normal = other_fluid.e_particles.Normal(their_nearest_nb);	
				//	VectorD my_normal = fluid.plateau_normals[nearest_nb];
				//	VectorD their_normal = other_fluid.plateau_normals[their_nearest_nb];
				//	real alignedness = my_normal.dot(their_normal); // the bigger (closer to 1) the better (easier to transfer)
				//	if (alignedness < 0.) distance_ratio *= 1. + alignedness; // when compleltely opposite, then don't transfer
				//	// done
				//	//velocity_percentage = distance_ratio;
				//	//VectorD my_normal = fluid.plateau_normals[nearest_nb];
				//	//VectorD their_normal = other_fluid.plateau_normals[their_nearest_nb];
				//	//real alignedness = my_normal.dot(their_normal); // the bigger (closer to 1) the better (easier to transfer)
				//	//alignedness = 0.5 * (alignedness + 1.); // a score between 0-1. 0 means opposite, 1 means aligned
				//	//real velocity_percentage = distance_ratio * alignedness;
				//	//VectorD unit_velocity = fluid.l_particles.V(i);
				//	//VectorD my_normal = fluid.e_particles.Normal(nearest_nb);
				//	//VectorD their_normal = fluid.e_particles.Normal(their_nearest_nb);
				//	////VectorD my_cross_product = unit_velocity.cross(my_normal);
				//	////VectorD their_cross_product = unit_velocity.cross(their_normal);
				//	////real score = my_cross_product.dot(their_cross_product.normalized()); //a score between -1 and 1
				//	//real score = my_normal.dot(their_normal);
				//	////real velocity_percentage = 0.5 * (score + 1.);
				//	//real velocity_percentage = std::max<real>(score, 0.);
				//	//velocity_percentage *= also_in[k].sa_ratio;
				//	//percentage = std::max<real>(soap_percentage, velocity_percentage);
				//	//percentage = soap_percentage;
				//	//percentage = velocity_percentage;
				//	percentage = soap_percentage + velocity_percentage;
				//	percentage = std::min<real>(percentage, 1.0);
				//	// multiply by sa_ratio
				//	//real weighted_percentage = std::min<real>(1., percentage * also_in[k].sa_ratio); //need to multiply by the confidence that this e particle is shared between the two regions
				//	if (percentage > most_likely_prob) {
				//		most_likely_prob = percentage;
				//		most_likely_region = j;
				//	}
				//}
				if (RandomFunc::Random_Real(0., 1.) < most_likely_prob) {
					next_owner[i][m] = most_likely_region;
				}
			}
		}
		Array<bool> has_modified; has_modified.resize(Num_Regions()); AuxFunc::Fill(has_modified, false);
		// after all particles to transfer are marked, add them to new owner
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			//Array<int> to_delete_l; to_delete_l.resize(fluid.l_particles.Size()); AuxFunc::Fill(to_delete_l, 0);
			// note that here I am using next_owner[i].size instead of l_particles.Size()
			// that is because, say, bubble B has 100 particles, and I added 10 from A to B
			// then I don't want to check the newly added particles, whether or not they need to be deleted or transferred
			// so I want to loop the first 100 rather than the first 110.
			// but this also assumes that the particles are added to the tail
			int number_transferred = 0;
			AuxFunc::Fill(fluid.l_particles.StatusRef(), 0);
			for (int m = 0; m < next_owner[i].size(); m++) {
				int next = next_owner[i][m];
				if (next != i) {
					// note that both "i" and "next" are modified
					has_modified[i] = true;
					has_modified[next] = true;
					number_transferred++;
					fluid.l_particles.Status(m) = 1;
					fluid.l_particles.WM(m) *= -1.; // since this L particle will be deleted, we need to compute the LOSS
					ParticleALEFilm<d>& next_fluid = *regions[next];
					int idx = next_fluid.l_particles.Add_Element();
					next_fluid.l_particles.Copy_Element_From(idx, fluid.l_particles, m);
					//better project the l particle here
					bool vb = false;
					bool has_nbs = true;
					int flg = next_fluid.e_particles.surface.Nearest_Geometry_KNN(next_fluid.l_particles.X(idx), next_fluid.l_particles.E(idx), 20, 7, vb);
					next_fluid.l_particles.V(idx) = AuxFunc::Eliminate_Unit_Component(next_fluid.l_particles.V(idx), next_fluid.l_particles.Normal(idx));
					VectorD nearby_e_V = next_fluid.e_particles.World_Avg_Value_KNN(next_fluid.l_particles.X(idx), next_fluid.e_particles.VRef(), next_fluid.neighbor_params.kernel, next_fluid.neighbor_params.interp_radius, 3, KernelType::QUINTIC, has_nbs);
					next_fluid.l_particles.V(idx) += AuxFunc::Component_Along(nearby_e_V, next_fluid.l_particles.Normal(idx));
					//next_fluid.l_particles.V(idx) = next_fluid.l_particles.World_Avg_Value_KNN(next_fluid.l_particles.X(idx), next_fluid.l_particles.VRef(), next_fluid.neighbor_params.kernel, next_fluid.neighbor_params.interp_radius, 20, KernelType::QUINTIC, has_nbs);
					//next_fluid.l_particles.H(idx) = next_fluid.l_particles.H(nearest_l_nb);
					next_fluid.l_particles.H(idx) = next_fluid.l_particles.World_Avg_Value_KNN(next_fluid.l_particles.X(idx), next_fluid.l_particles.HRef(), next_fluid.neighbor_params.kernel, next_fluid.neighbor_params.interp_radius, 1, KernelType::QUINTIC, has_nbs);
					next_fluid.l_particles.Status(idx) = 0; //make sure that the newly-received ones don't get deleted instantly.
					// update WM
					Array<int> nbs = next_fluid.e_particles.nbs_searcher->Find_Neighbors(next_fluid.l_particles.X(idx), next_fluid.neighbor_params.interp_radius);
					real total_w = 0.;
					for (int k = 0; k < nbs.size(); k++) {
						int e_nb = nbs[k];
						total_w += next_fluid.neighbor_params.kernel.template Weight<d>(next_fluid.l_particles.X(idx) - next_fluid.e_particles.X(e_nb), next_fluid.neighbor_params.interp_radius, KernelType::QUINTIC);
					}
					nbs = next_fluid.boundary.nbs_searcher->Find_Neighbors(next_fluid.l_particles.X(idx), next_fluid.neighbor_params.interp_radius);
					for (int k = 0; k < nbs.size(); k++) {
						int b_nb = nbs[k];
						total_w += next_fluid.neighbor_params.kernel.template Weight<d>(next_fluid.l_particles.X(idx) - next_fluid.boundary.X(b_nb), next_fluid.neighbor_params.interp_radius, KernelType::QUINTIC);
					}
					if (total_w > 0.) {
						next_fluid.l_particles.WM(idx) = 1. / total_w;
					}
					else {
						next_fluid.l_particles.WM(idx) = 0.;
					}
					// update WM done
				}
				else {
					fluid.l_particles.WM(m) = 0.;//this is to prepare for the accumulation of volumn gained/lost to the e_particles during this update.
				}
			}
			//std::cout << "Done" << std::endl;
			//don't delete l particles here just yet. Because the next bubble might come back looking for the particle index in your particles
			//you should keep it and delete everything later.
			//std::cout << "for region: " << i << " number transferred is: " << number_transferred << std::endl;
		}
		// For optimized performance, we may want to update only when particles are added/deleted to this bubble
		for (int i = 0; i < Num_Regions(); i++) {
			if (!has_modified[i]) continue;
			ParticleALEFilm<d>& fluid = *regions[i];
			int temp_num = 0;
			int temp_num_1 = 0;
			int temp_num_2 = 0;
			for (int m = 0; m < fluid.l_particles.Size(); m++) {
				if (fluid.l_particles.WM(m) < 0.) temp_num++;
				if (fluid.l_particles.WM(m) == 0.) temp_num_1++;
				if (fluid.l_particles.WM(m) > 0.) temp_num_2++;
			}
			//std::cout << "for fluid: " << i << " num with negative WM: " << temp_num << std::endl;
			//std::cout << "for fluid: " << i << " num with 0 WM: " << temp_num_1 << std::endl;
			//std::cout << "for fluid: " << i << " num with positive WM: " << temp_num_2 << std::endl;

			// before we delete we need to fix the L particle heights. With the INWARD or OUTWARD flux, the E particles will naturally rise or decrease in height
			// but the L particles don't know that. So we need to let them know
			fluid.Update_L();
			Array<real> E_changed_vols; E_changed_vols.resize(fluid.e_particles.Size()); AuxFunc::Fill(E_changed_vols, 0.);
			fluid.Distribute_L2E(E_changed_vols, fluid.l_particles.VolRef(), fluid.neighbor_params.interp_radius, KernelType::QUINTIC);
			//std::cout << "Total next influx for fluid: " << i << " is: " << AuxFunc::Sum(E_changed_vols) << std::endl;
#pragma omp parallel for
			for (int m = 0; m < fluid.e_particles.Size(); m++) {
				E_changed_vols[m] /= fluid.e_particles.SA(m);
				//if (E_changed_vols[m] < 0) std::cout << "yessir" << std::endl;
				fluid.e_particles.H(m) += E_changed_vols[m];
				//fluid.temp_vector[m] = 1000000. * E_changed_vols[m] * fluid.e_particles.Normal(m);
			}
			fluid.XSPH_Smoothing(0.99, E_changed_vols, fluid.neighbor_params.kernel, fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
#pragma omp parallel for
			for (int m = 0; m < fluid.l_particles.Size(); m++) {
				/*bool has_nbs;
				fluid.l_particles.H(m) += fluid.e_particles.World_Avg_Value(fluid.l_particles.X(m), E_changed_vols, fluid.neighbor_params.kernel, fluid.neighbor_params.sph_radius, KernelType::QUINTIC, has_nbs);*/
				fluid.l_particles.H(m) += fluid.e_particles.Surface_Convolution_Value(fluid.l_particles.X(m), fluid.l_particles.E(m), E_changed_vols, fluid.neighbor_params.kernel, fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
			}
#pragma omp parallel for
			for (int m = 0; m < fluid.l_particles.Size(); m++) {
				real e_h = fluid.e_particles.Surface_Convolution_Value(fluid.l_particles.X(m), fluid.l_particles.E(m), fluid.e_particles.HRef(), fluid.neighbor_params.kernel, fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
				real clipped_h = std::min<real>(fluid.l_particles.H(m), 2. * e_h);
				clipped_h = std::max<real>(clipped_h, .5 * e_h);
				real near_boundary_extent = 0.;
				int nearest_e_nb = fluid.e_particles.nbs_searcher->Find_Nearest_Nb(fluid.l_particles.X(m));
				for (int u = 0; u < master_list[i][nearest_e_nb].size(); u++) {
					near_boundary_extent = std::max<real>(near_boundary_extent, master_list[i][nearest_e_nb][u].sa_ratio);
				}
				near_boundary_extent = std::min<real>(near_boundary_extent, 1.);
				//fluid.l_particles.H(m) = (1. - near_boundary_extent) * fluid.l_particles.H(m) + near_boundary_extent * clipped_h;
				fluid.l_particles.H(m) = (1. - near_boundary_extent) * fluid.l_particles.H(m) + near_boundary_extent * e_h;
			}
		}
		for (int i = 0; i < Num_Regions(); i++) {
			if (!has_modified[i]) continue;
			ParticleALEFilm<d>& fluid = *regions[i];
			fluid.l_particles.Delete_Elements_Safe(fluid.l_particles.StatusRef());
			fluid.Update_L();
			fluid.Reinit_Arrays(false, true, false);
		}


		int total_l_particles = 0;
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			total_l_particles += fluid.l_particles.Size();
			std::cout << "for region: " << i << " number of L particles is: " << fluid.l_particles.Size() << std::endl;
		}
		std::cout << "total number of L particles is: " << total_l_particles << std::endl;
	}

	// update max velocity
	void Decay_Velocity(const real dt)
	{
		//for each bubble
		for (int i = 0; i < Num_Regions(); i++) {
			if (Is_Region_Trivial(i)) {
				continue;
			}
			ParticleALEFilm<d>& fluid = *regions[i];
			real tmp_decay_strength = decay_strength;
			//temporary!!
			//if (i == 2) {
			//	tmp_decay_strength = 2.;
			//}
			//temporary!!
#pragma omp parallel for
			for (int m = 0; m < fluid.e_particles.Size(); m++) {
				fluid.e_particles.V(m) *= exp(-tmp_decay_strength * dt);
			}
		}
	}

	// update max velocity
	void Update_All_BBox(void)
	{
		region_mins.clear();
		region_maxs.clear();
		// update max velocity of all regions
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			if (Is_Region_Trivial(i)) {
				region_mins.push_back(VectorD::Zero());
				region_maxs.push_back(VectorD::Zero());
				continue;
			}
			VectorD bbox_min, bbox_max;
			VectorD margin; margin << sph_radius, sph_radius, sph_radius;
			AuxFunc::Min_And_Max(fluid.e_particles.XRef(), bbox_min, bbox_max);
			region_mins.push_back(bbox_min-margin);
			region_maxs.push_back(bbox_max+margin);
		}
	}

	// update max velocity
	void Update_All_Geometry(const real dt)
	{
		region_vols.clear();
		region_sas.clear();
		// update max velocity of all regions
		real total_enclosed_vol = 0.;
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			if (Is_Region_Trivial(i)) {
				region_vols.push_back(0.);
				region_sas.push_back(0.);
				continue;
			}
			fluid.Update_E_Geometry();
			if (!e_only) fluid.Update_L_Geometry(dt);
			//compute enclosed volume
			real vol = fluid.Compute_Enclosed_Volume();
			real sa = AuxFunc::Sum(fluid.e_particles.SARef());
			region_vols.push_back(vol);
			total_enclosed_vol += vol;
			region_sas.push_back(sa);
		}
		Update_All_BBox();
		curr_COM = Compute_All_COM();
		curr_COM_velocity = Compute_All_COM_Velocity();
		//::cout << "Curr total enclosed_vol: \n" << std::setprecision(5) << total_enclosed_vol << std::endl;
		//std::cout << "Curr COM: \n" << std::setprecision(5) << curr_COM << std::endl;
		//std::cout << "Curr COM_velocity: \n" << std::setprecision(5) << curr_COM_velocity << std::endl;
	}

	// update max velocity
	void Update_All_Dynamics(const real dt) {
		Update_All_Norm_Dynamics(dt);
		//Apply_Boundary_Adhesion(dt); //adhere to the boundary particles
		if (!e_only) Update_All_Tang_Dynamics(dt);
		// update velocity of all regions
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			if (Is_Region_Trivial(i)) continue;
			fluid.XSPH_Smoothing_Tang_Norm(fluid.numeric_params.XSPH_V_Tang, fluid.numeric_params.XSPH_V_Norm, fluid.e_particles.VRef(), fluid.neighbor_params.kernel, fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
			for (int i = 0; i < fluid.numeric_params.XSPH_V_Passes; i++) {
				fluid.XSPH_Smoothing_Tang_Norm(0., fluid.numeric_params.XSPH_V_Passes_Strength, fluid.e_particles.VRef(), fluid.neighbor_params.kernel, fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
			}
			//Combine normal and tang velocities to form LV
#pragma omp parallel for
			for (int i = 0; i < fluid.e_particles.Size(); i++) { // for each E particle
				fluid.e_particles.LV(i) = AuxFunc::Component_Along(fluid.e_particles.V(i), fluid.e_particles.Normal(i)) + fluid.e_particles.LV_T(i);
			}
//			//Add artificial velocity to test only
//#pragma omp parallel for
//			for (int i = 0; i < fluid.e_particles.Size(); i++) { // for each E particle
//				VectorD vec = VectorD::Unit(1);
//				if (fluid.e_particles.X(i)[0] < 0.) { vec *= -1; }
//				vec *= 0.1;
//				fluid.e_particles.LV(i) += dt * AuxFunc::Eliminate_Unit_Component(vec, fluid.e_particles.Normal(i));
//			}
		}
	}

	void Update_All_Tang_Dynamics(const real dt)
	{
		// update max velocity of all regions
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			if (Is_Region_Trivial(i)) continue;
			fluid.Update_Tang_Dynamics(dt);
		}
	}

	// update norm dynamics (accounting for multi-bubble interactions)
	void Update_All_Norm_Dynamics(const real dt)
	{
		// update dynamics for all regions
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			if (Is_Region_Trivial(i)) continue;
			//std::cout << "For region: " << i << " enclosed amount is: " << fluid.enclosed_amount << std::endl;
			//clear force
			AuxFunc::Fill(fluid.e_particles.FRef(), VectorD::Zero());
			//compute enclosed volume
			real avg_sa = region_sas[i] / fluid.e_particles.Size();
			//std::cout << "AVG SA IS: " << avg_sa << std::endl;
			real vol = region_vols[i];
			real pressure = (vol > 1e-8) ? (fluid.enclosed_amount) / vol : outer_pressure;
			//First add the surface tension (inside layer) and enclosed pressure that belong to themselves (not considering guests)
#pragma omp parallel for
			for (int m = 0; m < fluid.e_particles.Size(); m++) {
				//surface tension (one half of it), since we are only considering the inner interface
				real pressure_difference = 100 * .5 * fluid.e_particles.ST(m) * fluid.e_particles.KS(m);//capillary pressure difference
				VectorD force = pressure_difference * avg_sa * fluid.e_particles.Normal(m);
				fluid.e_particles.F(m) += force;
				pressure_difference = pressure; // note that we are not subtracting the outer pressure, since we are considering only the inner interface
				force = pressure_difference * avg_sa * fluid.e_particles.Normal(m);
				fluid.e_particles.F(m) += force;
			}
			//also we consider the other interfaces
#pragma omp parallel for
			for (int m = 0; m < fluid.e_particles.Size(); m++) {
				Array<SharingInfo<d>>& also_in = master_list[i][m];
				for (int k = 0; k < also_in.size(); k++) { //loop over all other particles
					int j = also_in[k].region_idx;
					ParticleALEFilm<d>& other_fluid = *regions[j];
					VectorD proj_pos = also_in[k].proj_pos;
					MatrixD proj_frame = also_in[k].proj_frame;
					VectorD proj_normal = proj_frame.col(2);
					int nearest_nb = other_fluid.e_particles.nbs_searcher->Find_Nearest_Nb(fluid.e_particles.X(m));
					real KS = other_fluid.e_particles.KS(nearest_nb);
					real pressure_difference = 100 * 0.5 * fluid.e_particles.ST(m) * KS;//capillary pressure difference
					VectorD force = pressure_difference * avg_sa * proj_normal;
					force = AuxFunc::Component_Along(force, fluid.e_particles.Normal(m));
					fluid.e_particles.F(m) += force;
					//Then add the pressure force
					real vol = region_vols[j]; // the enclosed volume of the other fluid
					real pressure = (vol > 1e-8) ? (other_fluid.enclosed_amount) / vol : outer_pressure;
					force = pressure * avg_sa * proj_normal;
					force = AuxFunc::Component_Along(force, fluid.e_particles.Normal(m));
					fluid.e_particles.F(m) += force * also_in[k].sa_ratio;
				}
			}
//			//interface with air
#pragma omp parallel for
			for (int m = 0; m < fluid.e_particles.Size(); m++) {
				Array<SharingInfo<d>>& also_in = master_list[i][m];
				VectorD ext_normal = fluid.e_particles.Normal(m);
				bool is_inside_any = false;
				for (int q = 0; q < also_in.size(); q++) {
					if (also_in[q].signed_distance > 0. * e_dx) {
						is_inside_any = true;
						break;
					}
				}

				fluid.plateau_normals[m] = fluid.e_particles.Normal(m);
				if (also_in.size() < 1 || is_inside_any) {
					fluid.e_particles.F(m) += 100 * 0.5 * fluid.e_particles.ST(m) * avg_sa * fluid.e_particles.KS(m) * ext_normal;
					fluid.e_particles.F(m) += outer_pressure * avg_sa * -ext_normal;
					fluid.plateau_normals[m] = VectorD::Zero(); // no plateau normal if is not near plateau boarder or is inside any other bubble
				}
				else {
					for (int k = 0; k < also_in.size(); k++) { //loop over all other particles
						int j = also_in[k].region_idx;
						ext_normal += also_in[k].proj_frame.col(2) * also_in[k].sa_ratio;
					}
					fluid.plateau_normals[m] = ext_normal;
					//real norm = ext_normal.norm();
					//if (norm > 1.e-8)
					//	fluid.plateau_normals[m] = ext_normal / norm;
					//else
					//	fluid.plateau_normals[m] = VectorD::Zero();
					VectorD KS_Vector = Compute_Plateau_Curvature(i, m, ext_normal);
					VectorD ext_ST = 100 * 0.5 * fluid.e_particles.ST(m) * avg_sa * KS_Vector;
					VectorD totality = ext_ST + (outer_pressure * avg_sa * -ext_normal);
					fluid.e_particles.F(m) += AuxFunc::Component_Along(totality, fluid.e_particles.Normal(m));
				}
			}
			VectorD avg_F = AuxFunc::Mean(fluid.e_particles.FRef());
			//std::cout << "for i: " << i << " avg_F is: " << avg_F << std::endl;
			//std::cout << "for i: " << i << " total velocity before is: " << AuxFunc::Sum(fluid.e_particles.VRef()) << std::endl;
#pragma omp parallel for
			for (int m = 0; m < fluid.e_particles.Size(); m++) { // for each E particle
				fluid.e_particles.V(m) += (fluid.e_particles.F(m) - avg_F) / fluid.e_particles.M(m) * dt;
			}
			if (e_feels_gravity) {
				VectorD temp_g = g;
#pragma omp parallel for
				for (int m = 0; m < fluid.e_particles.Size(); m++) { // for each E particle
					fluid.e_particles.V(m) += temp_g * dt; //add in gravity
				}
			}
		}
	}
	
	real Adhesion_Kernel(const real& r, const real& h) const {
		//real alpha = 0;
		//if (r < h && 2 * r > h) {//0~0.5h
		//	alpha = 0.007 / pow(h, 3.25) * pow(-4 * pow(r, 2) / h + 6 * r - 2 * h, 1. / 4.);
		//}
		//return alpha;
		real alpha = 0;
		real tmp = pow(1. - r / h, 2.);
		return tmp;
	}


	// update norm dynamics (accounting for multi-bubble interactions)
	void Apply_Boundary_Adhesion(const real dt)
	{
//		// update dynamics for all regions
//		for (int i = 0; i < Num_Regions(); i++) {
//			ParticleALEFilm<d>& fluid = *regions[i];
//			if (Is_Region_Trivial(i)) continue;
//			real adhesion_radius = 3. * e_dx;
//			real adhesion_param = .33;
//			Array<VectorD> adhesions; 
//			adhesions.resize(fluid.e_particles.Size());
//			if (analytical_boundary.Available()) { fluid.e_particles.Update_Phis(analytical_boundary); }
//#pragma omp parallel for
//			for (int m = 0; m < fluid.e_particles.Size(); m++) {
//				VectorD aggregate_adhesion = VectorD::Zero();
//				Array<int> boundary_nbs = boundary.nbs_searcher->Find_Neighbors(fluid.e_particles.X(m), adhesion_radius);
//				for (int k = 0; k < boundary_nbs.size(); k++) {
//					int j = boundary_nbs[k];
//					VectorD r_ij = fluid.e_particles.X(m) - boundary.X(j);
//					aggregate_adhesion += -adhesion_param * Adhesion_Kernel(r_ij.norm(), adhesion_radius) * r_ij.normalized();
//				}
//
//				if (analytical_boundary.Available()) {
//					real friction_radius = fluid.neighbor_params.e_dx;
//					real closest_dist = abs(fluid.e_particles.phis[m]);
//					if (closest_dist < friction_radius) {
//						real ratio = closest_dist / friction_radius;
//						real multiplier = (1. - ratio); //1 when 0, 0 when friction_radius
//						real decay_strength = 2.5;
//						fluid.e_particles.V(m) *= exp(-decay_strength * multiplier * dt);
//					}
//				}
//
//				adhesions[m] = AuxFunc::Component_Along(aggregate_adhesion, fluid.e_particles.Normal(m));
//			}
//			fluid.XSPH_Smoothing(0.999, adhesions, fluid.neighbor_params.kernel, fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
////#pragma omp parallel for
////			for (int m = 0; m < fluid.e_particles.Size(); m++) {
////				fluid.temp_vector[m] = adhesions[m];
////				fluid.e_particles.V(m) += adhesions[m] * dt;
////			}
//		}
	}

	void Shift_Positions(Array<VectorD>& positions)
	{
		real strength = 0.1;
		VectorD offset = VectorD::Zero();
		for (int axis = 0; axis < d; axis++) {
			if (std::count(conserve_momentum_along.begin(), conserve_momentum_along.end(), axis)) {
				offset[axis] = init_COM[axis] - curr_COM[axis];
			}
		}
		offset *= strength;
#pragma omp parallel for
		for (int m = 0; m < positions.size(); m++) {
			positions[m] += offset;
		}
	}

	void Shift_Velocities(Array<VectorD>& velocities, real dt = 0., int i = -1)
	{
		real strength = 1.;
		VectorD offset = VectorD::Zero();
		VectorD damping_velocity = curr_COM_velocity;
		for (int axis = 0; axis < d; axis++) {
			if (std::count(conserve_momentum_along.begin(), conserve_momentum_along.end(), axis)) {
				offset[axis] = init_COM[axis] - curr_COM[axis];
			}
			else {
				damping_velocity[axis] = 0.;
			}
		}
		offset *= strength;
#pragma omp parallel for
		for (int m = 0; m < velocities.size(); m++) {
			velocities[m] += offset;
			velocities[m] -= 0.33 * damping_velocity;
		}
		if (i >= 0) {
//			ParticleALEFilm<d>& fluid = *regions[i];
//			real decay_strength = 0.1;
////			if (boundary.Size()) {
////				real bubble_bottom = region_mins[i](1) + (sph_radius + e_dx);
////				real bubble_left = region_mins[i](0) + (sph_radius + e_dx);
////				real bubble_right = region_maxs[i](0) - (sph_radius + e_dx);
////				real bubble_back = region_mins[i](2) + (sph_radius + e_dx);
////				real bubble_front = region_maxs[i](2) - (sph_radius + e_dx);
////				//check if a bubble is below top
////				if ((bubble_back > back && bubble_front < front) &&
////					(bubble_left > left && bubble_right < right)) {
////					decay_strength = 0.12;
////				}	
////			}
//#pragma omp parallel for
//			for (int m = 0; m < velocities.size(); m++) { 
//				velocities[m] -= decay_strength * dt * fluid.curr_COM_velocity;
//			}
		}
	}

	// update max velocity
	void Advance_All(const real dt, const int iter)
	{
		// update max velocity of all regions
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			if (Is_Region_Trivial(i)) continue;
			/*fluid.Advance_E(dt, iter);*/

			// Advance E part
			// Update position of each E particle
			Shift_Velocities(fluid.e_particles.VRef(), dt, i);
//			//temporary!!!!!!!!
//#pragma omp parallel for
//			for (int m = 0; m < fluid.e_particles.Size(); m++) {
//				// real highest = 20. * fluid.neighbor_params.e_dx;
//				//if (fluid.e_particles.X(m)[1] > 0 && fluid.e_particles.X(m)[1] < highest) {
//				//	real ratio = (fluid.e_particles.X(m)[1])/highest;
//				//	VectorD norm_part = AuxFunc::Component_Along(fluid.e_particles.V(m), fluid.e_particles.Normal(m));
//				//	VectorD tang_part = fluid.e_particles.V(m) - norm_part;
//				//	fluid.e_particles.V(m) = norm_part + tang_part * pow(ratio, 0.5);
//				//}
//				real highest = 10. * fluid.neighbor_params.e_dx;
//				if (fluid.e_particles.X(m)[1] > 0 && fluid.e_particles.X(m)[1] < highest) {
//					real ratio = (fluid.e_particles.X(m)[1]) / highest;
//					real multiplier = 1. - pow(1. - ratio, 3.); //0 when 0, 1 when highest
//					fluid.e_particles.V(m) *= 0.8 + 0.2 * multiplier;
//				}
//			}
//			//temporary!!!!!!!!
#pragma omp parallel for
			for (int m = 0; m < fluid.e_particles.Size(); m++) {
				fluid.e_particles.X(m) += fluid.e_particles.V(m) * dt;
			}
			//Shift_Positions(fluid.e_particles.XRef());
			if (analytical_boundary.Available()) fluid.e_particles.Correct_Position_With_Analytical_Boundary(dt, simulation_scale, analytical_boundary);
			//Manual_Enforce_Wall_BC();
			// Update E NB searcher
			fluid.Update_E();
			// Even_out
			if (!fluid.bursting) { 
				//std::cout << "This is called!" << std::endl;
				fluid.Even_Out_E(dt, false, fluid.max_even_out_e_steps); 
			} //if bursting or fixed_surface then we don't even out E
			else{
				//std::cout << "This is not called!" << std::endl;
			}

			if (!e_only) {
				// Advance L part
				//fluid.Advance_L(dt);
				// Update position of each L particle
				Shift_Velocities(fluid.l_particles.VRef(), dt, i);
#pragma omp parallel for
				for (int m = 0; m < fluid.l_particles.Size(); m++) {
					fluid.l_particles.X(m) += fluid.l_particles.V(m) * dt;
				}
				//Shift_Positions(fluid.l_particles.XRef());
				fluid.Project_L(dt); // project L, unless the surface is perfectly 2d
				if (analytical_boundary.Available()) fluid.l_particles.Correct_Position_With_Analytical_Boundary(dt, simulation_scale, analytical_boundary);
				// Update L NB searcher
				fluid.Update_L();
				// Update Temperature
				if (temp_func != nullptr) {
#pragma omp parallel for
					for (int m = 0; m < fluid.l_particles.Size(); m++) {
						fluid.l_particles.Temp(m) += 1. / fluid.chem_params.heat_capacity_param * dt * (temp_func(fluid.l_particles.X(m)) - fluid.l_particles.Temp(m));
						//fluid.l_particles.Temp(m) = temp_func(fluid.l_particles.X(m));
					}
				}
			}
		}
	}

	// L2E for all
	void L2E_All()
	{
		// update max velocity of all regions
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			if (Is_Region_Trivial(i)) continue;
//#pragma omp parallel for
//			for (int i = 0; i < fluid.e_particles.Size(); i++) {
//				fluid.temp_LV[i] = fluid.e_particles.LV(i);
//			}
			fluid.L2E(); // transfer momentum, mass, etc
#pragma omp parallel for
			for (int i = 0; i < fluid.e_particles.Size(); i++) {
				fluid.temp_LV[i] = fluid.e_particles.LV(i);
				//fluid.temp_LV[i] *= 0.5;
			}
		}
	}

	// L2E for all
	void E2L_All()
	{
		// update max velocity of all regions
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			if (Is_Region_Trivial(i)) continue;
			fluid.E2L();
		}
	}

	void Defreeze_Bubbles(void)
	{
		if (!bubble_defreeze_func) return;
		freeze = bubble_defreeze_func(curr_frame);
	}

	void Merge_Bubbles(void)
	{
		if (!bubble_merge_func) return;
		bubble_merge_func(curr_frame);
	}

	void Merge_Two_Bubbles(int i, int j)
	{
		ParticleALEFilm<d>& fluid_this = *regions[i];
		ParticleALEFilm<d>& fluid_that = *regions[j];
		if (fluid_that.Is_Trivial()) return;
		std::shared_ptr<ParticleALEFilm<d>> new_fluid_ptr = std::make_shared<ParticleALEFilm<d>>();
		//transfer points
		for (int m = 0; m < fluid_this.e_particles.Size(); m++) {
			int idx = new_fluid_ptr->e_particles.Add_Element();
			new_fluid_ptr->e_particles.Copy_Element_From(idx, fluid_this.e_particles, m);
		}
		for (int m = 0; m < fluid_that.e_particles.Size(); m++) {
			int idx = new_fluid_ptr->e_particles.Add_Element();
			new_fluid_ptr->e_particles.Copy_Element_From(idx, fluid_that.e_particles, m);
		}
		if (!e_only) {
			for (int m = 0; m < fluid_this.l_particles.Size(); m++) {
				int idx = new_fluid_ptr->l_particles.Add_Element();
				new_fluid_ptr->l_particles.Copy_Element_From(idx, fluid_this.l_particles, m);
			}
			for (int m = 0; m < fluid_that.l_particles.Size(); m++) {
				int idx = new_fluid_ptr->l_particles.Add_Element();
				new_fluid_ptr->l_particles.Copy_Element_From(idx, fluid_that.l_particles, m);
			}
		}
		//transfer points
		//set up fluid
		new_fluid_ptr->simulation_scale = fluid_this.simulation_scale;
		new_fluid_ptr->simulation_fineness = fluid_this.simulation_fineness;
		new_fluid_ptr->verbose = fluid_this.verbose;
		real e_dx = fluid_this.neighbor_params.e_dx;
		real e_np_on_h = 4;
		real interp_np_on_h = 4;
		new_fluid_ptr->neighbor_params = NeighborParams<d>(e_dx, e_np_on_h, interp_np_on_h);
		//---configure chem and numeric params---//
		new_fluid_ptr->chem_params = ChemicalParams();
		new_fluid_ptr->numeric_params = NumericalParams();
		new_fluid_ptr->fluid_volume = fluid_this.fluid_volume + fluid_that.fluid_volume;
		new_fluid_ptr->fluid_area = fluid_this.fluid_area + fluid_that.fluid_area;
		new_fluid_ptr->fluid_soap = fluid_this.fluid_soap + fluid_that.fluid_soap;
		new_fluid_ptr->l_geometry_sync = fluid_this.l_geometry_sync;
		new_fluid_ptr->g = g;
		new_fluid_ptr->numeric_params.XSPH_V_Passes = fluid_this.numeric_params.XSPH_V_Passes;
		new_fluid_ptr->numeric_params.XSPH_V_Passes_Strength = fluid_this.numeric_params.XSPH_V_Passes_Strength;
		new_fluid_ptr->numeric_params.XSPH_KS_Passes = fluid_this.numeric_params.XSPH_KS_Passes;
		new_fluid_ptr->e_particles.base_mass = fluid_this.e_particles.base_mass;
		new_fluid_ptr->e_particles.base_rho = fluid_this.e_particles.base_rho;
		//init
		new_fluid_ptr->Initialize(cfl, frame_rate, output_dir);
		//done init
		new_fluid_ptr->enclosed = 1;
		new_fluid_ptr->enclosed_vol = fluid_this.enclosed_vol + fluid_that.enclosed_vol;
		new_fluid_ptr->enclosed_amount = fluid_this.enclosed_amount + fluid_that.enclosed_amount;
		//degenerate old ones
		fluid_this.Degenerate();
		fluid_that.Degenerate();
		regions[i] = new_fluid_ptr;
		Reset_Master_List();
	}

	VectorD Compute_All_COM(void)
	{
		VectorD avg_scaled_position = VectorD::Zero();
		real total_SA = 0.;
		// update max velocity of all regions
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			if (Is_Region_Trivial(i)) continue;
			//std::cout << "i? " << i << std::endl;
			//std::cout << "COM? \n" << fluid.curr_COM << std::endl;
			//std::cout << "fluid area? \n" << fluid.fluid_area << std::endl;
			avg_scaled_position += fluid.curr_COM * fluid.fluid_area;
			total_SA += fluid.fluid_area;
		}
		if (total_SA >= 1.e-8) {
			avg_scaled_position /= total_SA;
		}
		else {
			avg_scaled_position = VectorD::Zero();
		}
		return avg_scaled_position;
	}

	VectorD Compute_All_COM_Velocity(void)
	{
		VectorD avg_scaled_velocity = VectorD::Zero();
		real total_SA = 0.;
		// update max velocity of all regions
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			if (Is_Region_Trivial(i)) continue;
			avg_scaled_velocity += fluid.curr_COM_velocity * fluid.fluid_area;
			total_SA += fluid.fluid_area;
		}
		if (total_SA >= 1.e-8) {
			avg_scaled_velocity /= total_SA;
		}
		else {
			avg_scaled_velocity = VectorD::Zero();
		}
		return avg_scaled_velocity;
	}

	void Prune_Self_Intersections(void)
	{
		for (int i = 0; i < Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *regions[i];
			if (Is_Region_Trivial(i)) continue;
			// now delete some particles facing each other
			AuxFunc::Fill(fluid.to_delete_e, 0);
			real threshold = 1.1e5;
			Array<real> values; values.resize(fluid.e_particles.Size());
#pragma omp parallel for
			for (int i = 0; i < fluid.e_particles.Size(); i++) {
				std::function<VectorD(const int)> normal_func = [&](const int idx)->VectorD {return fluid.e_particles.Normal(idx); };
				VectorD normal_laplacian = fluid.e_particles.Surface_Laplacian(i, normal_func, fluid.neighbor_params.kernel, fluid.neighbor_params.sph_radius, KernelType::QUINTIC);
				values[i] = normal_laplacian.norm();
				if (values[i] > threshold) fluid.to_delete_e[i] = 1;
			}
			int num_deletion = AuxFunc::Sum(fluid.to_delete_e);
			if (num_deletion) {
				int remaining_e = fluid.e_particles.Delete_Elements_Safe(fluid.to_delete_e);
				fluid.Update_E();
				fluid.Reinit_Arrays(true, false, false);
			}
			//std::cout << "max norm: " << AuxFunc::Max(values) << std::endl;
			//std::cout << "mean norm: " << AuxFunc::Mean(values) << std::endl;
		}
	}

	void Burst_Bubbles(void)
	{
		//bool detonate = false;
		//if (num_detonations_done < num_detonations) {
		//	detonate = true;
		//	num_detonations_done++;
		//}
		//for (int i = 0; i < Num_Regions(); i++) {
		//	ParticleALEFilm<d>& fluid = *regions[i];
		//	if (Is_Region_Trivial(i)) continue;
		//	bool to_burst = false;
		//	if (fluid.curr_nden_err > 4.) {
		//		to_burst = true;
		//	}
		//	if (boundary.Size()) {
		//		real bubble_bottom = region_mins[i](1) + (sph_radius + e_dx);
		//		real bubble_left = region_mins[i](0) + (sph_radius + e_dx);
		//		real bubble_right = region_maxs[i](0) - (sph_radius + e_dx);
		//		real bubble_back =  region_mins[i](2) + (sph_radius + e_dx);
		//		real bubble_front = region_maxs[i](2) - (sph_radius + e_dx);
		//		//check if a bubble is below top
		//		if (bubble_bottom < top) {
		//			if ((bubble_left < left && bubble_right > left && bubble_back > back && bubble_front < front) ||
		//				(bubble_left < right && bubble_right > right && bubble_back > back && bubble_front < front) ||
		//				(bubble_back < back && bubble_front > back && bubble_left > left && bubble_right < right) ||
		//				(bubble_back < front && bubble_front > front && bubble_left > left && bubble_right < right)) {
		//				to_burst = true;
		//			}
		//		}
		//		//burst near bottom
		//		if (detonate) {
		//			if ((bubble_back > back && bubble_front < front) &&
		//				(bubble_left > left && bubble_right < right)) {
		//				real distance_to_bottom = (fluid.curr_COM - (-0. * VectorD::Unit(1) + 0.5 * VectorD::Unit(2))).norm();
		//				real multiplier = std::max<real>(1. - distance_to_bottom, 0.) / 1.;
		//				if (RandomFunc::Random_Real(0., 1.) < .5 * multiplier) {
		//					to_burst = true;
		//				}
		//			}
		//			//if ((bubble_back > back && bubble_front < front) &&
		//			//	(bubble_left > left && bubble_right < right)) {
		//			//	real distance_to_center = sqrt(fluid.curr_COM[0] * fluid.curr_COM[0] + fluid.curr_COM[2] * fluid.curr_COM[2]);
		//			//	real max_distance_to_center = sqrt(bubble_right * bubble_right + bubble_back * bubble_back);
		//			//	real ratio = 1.-(distance_to_center / max_distance_to_center);
		//			//	//ratio = ratio * ratio;
		//			//	real bottom = -0.3;
		//			//	real multiplier = ratio * (1.-(fluid.curr_COM[1] - bottom)/(1.-bottom));
		//			//	if (RandomFunc::Random_Real(0., 1.) < 1.4 * multiplier) {
		//			//		to_burst = true;
		//			//	}
		//			//}
		//		}
		//	}
		//	if (to_burst) {
		//		fluid.Degenerate();
		//	}
		//}
	}

	void Manual_Enforce_Wall_BC(void)
	{
//		if (!boundary.Size()) return;
//		for (int i = 0; i < Num_Regions(); i++) {
//			ParticleALEFilm<d>& fluid = *regions[i];
//			if (Is_Region_Trivial(i)) continue;
//#pragma omp parallel for
//			for (int m = 0; m < fluid.e_particles.Size(); m++) {
//				if (fluid.e_particles.X(m)[1] >= top) {
//					//do nothing
//				}
//				else {
//					// 0->left, 1->right, 2->back, 3->front
//					Array<real> distances; 
//					distances.push_back(left-fluid.e_particles.X(m)[0]);//left distance
//					distances.push_back(fluid.e_particles.X(m)[0]-right);//right distance
//					distances.push_back(back - fluid.e_particles.X(m)[2]);//back distance
//					distances.push_back(fluid.e_particles.X(m)[2] - front);//front distance
//
//					int closest_boundary = -1;
//					real closest_dist = 0.;
//					for (int k = 0; k < distances.size(); k++) {
//						if (closest_boundary < 0 || distances[k] > closest_dist) {
//							closest_dist = distances[k];
//							closest_boundary = k;
//						}
//					}
//					real cushion = 0.5 * e_dx;
//					if (closest_dist > -cushion && closest_dist < cushion) {
//						if (closest_boundary == 0) {
//							VectorD normal = VectorD::Unit(0);
//							VectorD normal_velocity = fluid.e_particles.V(m).dot(normal) * normal;
//							VectorD tangential_velocity = fluid.e_particles.V(m) - normal_velocity;
//							tangential_velocity += std::max<real>(fluid.e_particles.V(m).dot(normal), 0.) * normal;
//							fluid.e_particles.V(m) = tangential_velocity;
//							if (closest_dist > 0) {
//								fluid.e_particles.X(m)[0] = left;
//							}
//						}
//						else if (closest_boundary == 1) {
//							VectorD normal = -VectorD::Unit(0);
//							VectorD normal_velocity = fluid.e_particles.V(m).dot(normal) * normal;
//							VectorD tangential_velocity = fluid.e_particles.V(m) - normal_velocity;
//							tangential_velocity += std::max<real>(fluid.e_particles.V(m).dot(normal), 0.) * normal;
//							fluid.e_particles.V(m) = tangential_velocity;
//							if (closest_dist > 0) {
//								fluid.e_particles.X(m)[0] = right;
//							}
//						}
//						else if (closest_boundary == 2) {
//							VectorD normal = VectorD::Unit(2);
//							VectorD normal_velocity = fluid.e_particles.V(m).dot(normal) * normal;
//							VectorD tangential_velocity = fluid.e_particles.V(m) - normal_velocity;
//							tangential_velocity += std::max<real>(fluid.e_particles.V(m).dot(normal), 0.) * normal;
//							fluid.e_particles.V(m) = tangential_velocity;
//							if (closest_dist > 0) {
//								fluid.e_particles.X(m)[2] = back;
//							}
//						}
//						else if (closest_boundary == 3) {
//							VectorD normal = -VectorD::Unit(2);
//							VectorD normal_velocity = fluid.e_particles.V(m).dot(normal) * normal;
//							VectorD tangential_velocity = fluid.e_particles.V(m) - normal_velocity;
//							tangential_velocity += std::max<real>(fluid.e_particles.V(m).dot(normal), 0.) * normal;
//							fluid.e_particles.V(m) = tangential_velocity;
//							if (closest_dist > 0) {
//								fluid.e_particles.X(m)[2] = front;
//							}
//						}
//					}
//				}
//			}
//		}
	}


	virtual void Advance(const real dt, const real current_time, const int iter)
	{
		Defreeze_Bubbles();
		Decay_Velocity(dt);
		Update_Particle_Sharing_Info();
		Update_Particle_Sharing_Momentum(dt);

		if (!e_only) L2E_All();
		Update_All_Geometry(dt);
		Update_All_Dynamics(dt);
		if (!e_only) E2L_All();

		Advance_All(dt, iter);
		Update_Max_Velocity();
		if (!e_only) Interchange_L_Particles(dt);
	
		//this is how to get rid of a bubble
		//Burst_Bubbles();
		//ParticleALEFilm<d>& fluid = *regions[0];
		//fluid.Degenerate();
		//
		Merge_Bubbles();
		Prune_Self_Intersections();
	}
};

#endif
