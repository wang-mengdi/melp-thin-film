#include "ParticleALEFilm.h"
#include "RandomNumber.h"

// apply black hole forces
template<int d>
void ParticleALEFilm<d>::Evaporate(const real dt)
{
#pragma omp parallel for
	for (int i = 0; i < l_particles.Size(); i++) {
		real rate = .1;
		real factor = 1 - dt * .1;
		if (temp_func && chem_params.T_0 > 0.) factor *= l_particles.Temp(i) / chem_params.T_0;
		l_particles.Vol(i) *= factor;
		l_particles.H(i) *= factor;
	}
}

// apply black hole forces
template<int d>
void ParticleALEFilm<d>::Apply_Ext_Forces(const real dt)
{
	if (!ext_acc_func) return;
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		VectorD ext_acc = ext_acc_func(i);
		ext_acc *= dt;
		e_particles.V(i) += ext_acc;
		e_particles.LV(i) += ext_acc;
		e_particles.LV_T(i) += AuxFunc::Eliminate_Unit_Component(ext_acc, e_particles.Normal(i));
	}
}

template<int d>
Vector<real, d> ParticleALEFilm<d>::Compute_COM(void)
{
	Array<VectorD> scaled_positions; scaled_positions.resize(e_particles.Size());
	AuxFunc::Fill(scaled_positions, VectorD::Zero());
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		scaled_positions[i] = e_particles.SA(i) * e_particles.X(i);
	}
	VectorD avg_scaled_position = AuxFunc::Sum(scaled_positions);
	avg_scaled_position /= fluid_area;
	return avg_scaled_position;
}

template<int d>
Vector<real, d> ParticleALEFilm<d>::Compute_COM_Velocity(void)
{
	Array<VectorD> scaled_velocities; scaled_velocities.resize(e_particles.Size());
	AuxFunc::Fill(scaled_velocities, VectorD::Zero());
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		scaled_velocities[i] = e_particles.SA(i) * e_particles.V(i);
	}
	VectorD avg_scaled_velocity = AuxFunc::Sum(scaled_velocities);
	avg_scaled_velocity /= fluid_area;
	return avg_scaled_velocity;
}

// update max velocity
template<int d>
void ParticleALEFilm<d>::Detect_BH(void)
{
	for (int i = 0; i < bh_seeders.size(); i++) {
		if (!bh_seeders[i].initialized) continue;
		// this is to determine when the seeder is active
		if (sin(2 * pi * (curr_frame/ bh_seeders[i].period + bh_seeders[i].phase)) < bh_seeders[i].threshold) continue;
		VectorD loc = bh_seeders[i].loc;
		real r = bh_seeders[i].r;
		Array<int> nbs = l_particles.nbs_searcher->Find_Neighbors(loc, r);
		for (int k = 0; k < nbs.size(); k++) {
			int j = nbs[k];
			if (l_particles.Is_BH(j) < 1.) {
				l_particles.Is_BH(j) = 1.; 
			}
		}
	}
#pragma omp parallel for
	for (int i = 0; i < l_particles.Size(); i++) {
		if (curr_frame > 0 && l_particles.H(i) <= 80 * 1.E-9) { //if is low in height, directly convert to black hole particle
			l_particles.Is_BH(i) = 1.;
		}
		else {
			real BH_ratio = 0.;
			Array<int> nbs = l_particles.nbs_searcher->Find_Neighbors(l_particles.X(i), 1. * neighbor_params.e_dx);
			if (nbs.size() >= 1) {
				for (int k = 0; k < nbs.size(); k++) {
					int j = nbs[k];
					if (l_particles.Is_BH(j) > 0.) {
						BH_ratio++;
					}
				}
				BH_ratio /= nbs.size();
			}
			if (BH_ratio > 0.66) l_particles.Is_BH(i) = 1.;
		}
	}
}


// update black hole ratio
template<int d>
void ParticleALEFilm<d>::Update_BH_Ratio(void)
{
	Array<real> all_particles; all_particles.resize(e_particles.Size()); AuxFunc::Fill(all_particles, 0.);
	Array<real> bh_particles; bh_particles.resize(e_particles.Size()); AuxFunc::Fill(bh_particles, 0.);
	Distribute_L2E(all_particles, ones_l, neighbor_params.interp_radius, KernelType::QUINTIC);
	Distribute_L2E(bh_particles, l_particles.Is_BHRef(), neighbor_params.interp_radius, KernelType::QUINTIC);
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		if (all_particles[i] <= 1.e-8) {
			e_particles.BH_Ratio(i) = 0.; // if there is no L particles around, then we say its not a black hole --- or maybe it is?
		}
		else {
			e_particles.BH_Ratio(i) = bh_particles[i] / all_particles[i];
		}
	}
	Detect_BH();
}

//// apply black hole forces
//template<int d>
//void ParticleALEFilm<d>::Apply_BH_Forces(const real dt)
//{
//	Array<VectorD> SNs; SNs.resize(e_particles.Size()); AuxFunc::Fill(SNs, VectorD::Zero());
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) {
//		if (e_particles.BH_Ratio(i) > 0.) {
//			VectorT SN_2d = e_particles.Surface_Gradient_Difference(i, e_particles.BH_RatioRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//			VectorD SN_3d; e_particles.surface.Unproject_To_World(SN_2d, e_particles.E(i), SN_3d);
//			// this is to make sure it is fineness invariant. Say a black hole has radius 0.1 in a rim with radius 1
//			// if we scale up to be twice the fineness, the black hole will have a radius 0.05, in a rim with radiu 1
//			// since everything is compressed, the gradient will be twice
//			// so we need to divide it by 2, and multiplying the sph radius will do that
//			SNs[i] = neighbor_params.sph_radius * SN_3d;
//		}
//	}
//	real cohesion_parameter = 300.;
//	real curvature_parameter = 22.;
//	real boundary_repulsion_parameter = 30.;
//	real white_multiplier = 0.2;
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) {
//		e_particles.BH_LV(i) = VectorD::Zero();
//		if (e_particles.BH_Ratio(i) > 0.) {
//			Array<int> nbs = e_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.interp_radius);
//			VectorD total_cohesion = VectorD::Zero();
//			for (int k = 0; k < nbs.size(); k++) {
//				int j = nbs[k];
//				VectorD r_ij = e_particles.X(i) - e_particles.X(j);
//				real blackness = 0.5 * (e_particles.BH_Ratio(i) + e_particles.BH_Ratio(j));
//				VectorD diff = e_particles.X(j) - e_particles.X(i);
//				if (diff.dot(SNs[i]) > 0) {
//					total_cohesion += blackness * (-2. * e_particles.Size()) / (e_particles.NDen(i) + e_particles.NDen(j)) * r_ij.normalized() * Cohesion_Kernel(r_ij.norm(), neighbor_params.interp_radius);
//				}
//				else if (diff.dot(SNs[i]) < 0) {
//					total_cohesion += (1. - e_particles.BH_Ratio(i)) * -1. * blackness * (-2. * e_particles.Size()) / (e_particles.NDen(i) + e_particles.NDen(j)) * r_ij.normalized() * Cohesion_Kernel(r_ij.norm(), neighbor_params.interp_radius);
//				}
//			}
//			e_particles.BH_LV(i) += dt * cohesion_parameter * total_cohesion; 
//
//			VectorD total_curvature = VectorD::Zero();
//			for (int k = 0; k < nbs.size(); k++) {
//				int j = nbs[k];
//				if (e_particles.BH_Ratio(j) > 0.) {
//					real blackness = 0.5 * (e_particles.BH_Ratio(i) + e_particles.BH_Ratio(j));
//					VectorD diff = e_particles.X(j) - e_particles.X(i);
//					if (diff.dot(SNs[i])>0) {
//						VectorD normal_diff = SNs[i] - SNs[j];
//						if (normal_diff.dot(SNs[i]) > 0)
//						total_curvature += blackness * (-2. * e_particles.Size()) / (e_particles.NDen(i) + e_particles.NDen(j)) * (SNs[i] - SNs[j]);
//					}
//				}
//			}
//			e_particles.BH_LV(i) += curvature_parameter * dt * total_curvature * (1.- e_particles.BH_Ratio(i));
//
//			VectorD total_cohesion2 = VectorD::Zero();
//			for (int k = 0; k < nbs.size(); k++) {
//				int j = nbs[k];
//				VectorD r_ij = e_particles.X(i) - e_particles.X(j);
//				real whiteness = 0.5 * ((1. - e_particles.BH_Ratio(i)) + (1. - e_particles.BH_Ratio(j)));
//				VectorD diff = e_particles.X(j) - e_particles.X(i);
//				if (diff.dot(-SNs[i]) > 0) {
//					total_cohesion2 += whiteness * (-2. * e_particles.Size()) / (e_particles.NDen(i) + e_particles.NDen(j)) * r_ij.normalized() * Cohesion_Kernel(r_ij.norm(), neighbor_params.interp_radius);
//				}
//				else if (diff.dot(SNs[i]) < 0) {
//					total_cohesion2 += (e_particles.BH_Ratio(i)) * -1. * whiteness * (-2. * e_particles.Size()) / (e_particles.NDen(i) + e_particles.NDen(j)) * r_ij.normalized() * Cohesion_Kernel(r_ij.norm(), neighbor_params.interp_radius);
//				}
//			}
//			e_particles.BH_LV(i) += white_multiplier * dt * total_cohesion2 * cohesion_parameter;
//			// surface tension for the non-BH phase
//			VectorD total_curvature2 = VectorD::Zero();
//			for (int k = 0; k < nbs.size(); k++) {
//				int j = nbs[k];
//				VectorD r_ij = e_particles.X(i) - e_particles.X(j);
//				real whiteness = 0.5 * ((1.-e_particles.BH_Ratio(i)) + (1.-e_particles.BH_Ratio(j)));
//				VectorD diff = e_particles.X(j) - e_particles.X(i);
//				if (diff.dot(-SNs[i]) > 0) {
//					total_curvature2 += whiteness * (-2. * e_particles.Size()) / (e_particles.NDen(i) + e_particles.NDen(j)) * ((-SNs[i]) - (-SNs[j]));
//				}
//			}
//			e_particles.BH_LV(i) += white_multiplier * dt * total_curvature2 * curvature_parameter * e_particles.BH_Ratio(i);
//
//			//add boundary repulsion
//			Array<int> b_nbs = boundary.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.interp_radius);
//			VectorD total_boundary_repulsion = VectorD::Zero();
//			for (int k = 0; k < b_nbs.size(); k++) {
//				int j = b_nbs[k];
//				VectorD r_ij = e_particles.X(i) - boundary.X(j);
//				real norm = r_ij.norm();
//				if (norm < 1.e-8) continue;
//				real blackness = e_particles.BH_Ratio(i);
//				total_boundary_repulsion += blackness * (1. * e_particles.Size()) / (e_particles.NDen(i)) * r_ij.normalized() * cos(norm / neighbor_params.interp_radius * pi/2.);
//			}
//			e_particles.BH_LV(i) += dt * boundary_repulsion_parameter * total_boundary_repulsion;
//			
//			e_particles.BH_LV(i) += dt * -1. * AuxFunc::Eliminate_Unit_Component(g, e_particles.Normal(i));
//		}
//	}
//	XSPH_Smoothing<VectorD>(0.99, e_particles.BH_LVRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
//		//first decay velocity near BH
//		//then add the BH forces
//		VectorD norm_component = e_particles.LV(i) - e_particles.LV_T(i);
//		e_particles.LV_T(i) += e_particles.BH_LV(i);
//		e_particles.LV_T(i) *= exp(-3 * (e_particles.BH_Ratio(i)) * dt);
//		e_particles.LV(i) = norm_component + e_particles.LV_T(i);
//	}
//}


// apply black hole forces
// NEW
template<int d>
void ParticleALEFilm<d>::Apply_BH_Forces(const real dt)
{
	Array<VectorD> SNs; SNs.resize(e_particles.Size()); AuxFunc::Fill(SNs, VectorD::Zero());
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		if (e_particles.BH_Ratio(i) > 0.) {
			VectorT SN_2d = e_particles.Surface_Gradient_Difference(i, e_particles.BH_RatioRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
			VectorD SN_3d; e_particles.surface.Unproject_To_World(SN_2d, e_particles.E(i), SN_3d);
			// this is to make sure it is fineness invariant. Say a black hole has radius 0.1 in a rim with radius 1
			// if we scale up to be twice the fineness, the black hole will have a radius 0.05, in a rim with radiu 1
			// since everything is compressed, the gradient will be twice
			// so we need to divide it by 2, and multiplying the sph radius will do that
			SNs[i] = neighbor_params.sph_radius * SN_3d;
		}
	}
	real cohesion_parameter = 300.;
	real curvature_parameter = 22.;
	real boundary_repulsion_parameter = 30.;
	real white_multiplier = 0.2;
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		e_particles.BH_LV(i) = VectorD::Zero();
		if (e_particles.BH_Ratio(i) > 0.) {
			Array<int> nbs = e_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.interp_radius);
			VectorD total_cohesion = VectorD::Zero();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				VectorD r_ij = e_particles.X(i) - e_particles.X(j);
				real blackness = 0.5 * (e_particles.BH_Ratio(i) + e_particles.BH_Ratio(j));
				VectorD diff = e_particles.X(j) - e_particles.X(i);
				if (diff.dot(SNs[i]) > 0) {
					total_cohesion += blackness * (-2. * e_particles.Size()) / (e_particles.NDen(i) + e_particles.NDen(j)) * r_ij.normalized() * Cohesion_Kernel(r_ij.norm(), neighbor_params.interp_radius);
				}
				else if (diff.dot(SNs[i]) < 0) {
					total_cohesion += (1. - e_particles.BH_Ratio(i)) * -1. * blackness * (-2. * e_particles.Size()) / (e_particles.NDen(i) + e_particles.NDen(j)) * r_ij.normalized() * Cohesion_Kernel(r_ij.norm(), neighbor_params.interp_radius);
				}
			}
			e_particles.BH_LV(i) += dt * cohesion_parameter * total_cohesion; 

			VectorD total_curvature = VectorD::Zero();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				if (e_particles.BH_Ratio(j) > 0.) {
					real blackness = 0.5 * (e_particles.BH_Ratio(i) + e_particles.BH_Ratio(j));
					VectorD diff = e_particles.X(j) - e_particles.X(i);
					if (diff.dot(SNs[i])>0) {
						VectorD normal_diff = SNs[i] - SNs[j];
						if (normal_diff.dot(SNs[i]) > 0)
						total_curvature += blackness * (-2. * e_particles.Size()) / (e_particles.NDen(i) + e_particles.NDen(j)) * (SNs[i] - SNs[j]);
					}
				}
			}
			e_particles.BH_LV(i) += curvature_parameter * dt * total_curvature * (1.- e_particles.BH_Ratio(i));

			VectorD total_cohesion2 = VectorD::Zero();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				VectorD r_ij = e_particles.X(i) - e_particles.X(j);
				real whiteness = 0.5 * ((1. - e_particles.BH_Ratio(i)) + (1. - e_particles.BH_Ratio(j)));
				VectorD diff = e_particles.X(j) - e_particles.X(i);
				if (diff.dot(-SNs[i]) > 0) {
					total_cohesion2 += whiteness * (-2. * e_particles.Size()) / (e_particles.NDen(i) + e_particles.NDen(j)) * r_ij.normalized() * Cohesion_Kernel(r_ij.norm(), neighbor_params.interp_radius);
				}
				else if (diff.dot(SNs[i]) < 0) {
					total_cohesion2 += (e_particles.BH_Ratio(i)) * -1. * whiteness * (-2. * e_particles.Size()) / (e_particles.NDen(i) + e_particles.NDen(j)) * r_ij.normalized() * Cohesion_Kernel(r_ij.norm(), neighbor_params.interp_radius);
				}
			}
			e_particles.BH_LV(i) += white_multiplier * dt * total_cohesion2 * cohesion_parameter;
			// surface tension for the non-BH phase
			VectorD total_curvature2 = VectorD::Zero();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				VectorD r_ij = e_particles.X(i) - e_particles.X(j);
				real whiteness = 0.5 * ((1.-e_particles.BH_Ratio(i)) + (1.-e_particles.BH_Ratio(j)));
				VectorD diff = e_particles.X(j) - e_particles.X(i);
				if (diff.dot(-SNs[i]) > 0) {
					total_curvature2 += whiteness * (-2. * e_particles.Size()) / (e_particles.NDen(i) + e_particles.NDen(j)) * ((-SNs[i]) - (-SNs[j]));
				}
			}
			e_particles.BH_LV(i) += white_multiplier * dt * total_curvature2 * curvature_parameter * e_particles.BH_Ratio(i);

			//add boundary repulsion
			Array<int> b_nbs = boundary.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.interp_radius);
			VectorD total_boundary_repulsion = VectorD::Zero();
			for (int k = 0; k < b_nbs.size(); k++) {
				int j = b_nbs[k];
				VectorD r_ij = e_particles.X(i) - boundary.X(j);
				real norm = r_ij.norm();
				if (norm < 1.e-8) continue;
				real blackness = e_particles.BH_Ratio(i);
				total_boundary_repulsion += blackness * (1. * e_particles.Size()) / (e_particles.NDen(i)) * r_ij.normalized() * cos(norm / neighbor_params.interp_radius * pi/2.);
			}
			e_particles.BH_LV(i) += dt * boundary_repulsion_parameter * total_boundary_repulsion;
			
			e_particles.BH_LV(i) += dt * -1. * AuxFunc::Eliminate_Unit_Component(g, e_particles.Normal(i));
		}
	}
	XSPH_Smoothing<VectorD>(0.99, e_particles.BH_LVRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
		//first decay velocity near BH
		//then add the BH forces
		VectorD norm_component = e_particles.LV(i) - e_particles.LV_T(i);
		e_particles.LV_T(i) += e_particles.BH_LV(i);
		e_particles.LV_T(i) *= exp(-3 * (e_particles.BH_Ratio(i)) * dt);
		e_particles.LV(i) = norm_component + e_particles.LV_T(i);
	}
}



// update max velocity
template<int d>
void ParticleALEFilm<d>::Decay_Velocity(Array<VectorD>& vels, const real dt, const real tang_strength, const real norm_strength)
{
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		VectorD norm_V = AuxFunc::Component_Along(vels[i], e_particles.Normal(i));
		VectorD tang_V = vels[i] - norm_V;
		norm_V *= exp(-norm_strength * dt);
		tang_V *= exp(-tang_strength * dt);
		vels[i] = norm_V + tang_V;
	}
}

template<int d>
void ParticleALEFilm<d>::Update_WM(real radius)
{
	//for each L particle,
	//loop through each E particle that "claims possession of it" (will find this L using the E's kernel)
	//and gauge the relationship (weight)
	//finally decide oa multiplier of how much to give to it
	//for example, a L particle has 1 unit of mass,
	//there are 3 E particles: E1, E2, E3 that will find this L
	//they each evaluate to a relationship (weight): 200, 300, 500
	//then the value stored will be 1/(200 + 300 + 500) = 0.001
	//so that in the time of actual distribution, E1 will get 0.001 = 0.2, E2 will get 0.3, E3 will get 0.5, and the total sum is 1
#pragma omp parallel for
	for (int i = 0; i < l_particles.Size(); i++) {
		Array<int> nbs = e_particles.nbs_searcher->Find_Neighbors(l_particles.X(i), radius);
		real total_w = 0.;
		for (int k = 0; k < nbs.size(); k++) {
			int e_nb = nbs[k];
			total_w += neighbor_params.kernel.template Weight<d>(l_particles.X(i) - e_particles.X(e_nb), radius, KernelType::QUINTIC);
		}
		nbs = boundary.nbs_searcher->Find_Neighbors(l_particles.X(i), radius);
		for (int k = 0; k < nbs.size(); k++) {
			int b_nb = nbs[k];
			total_w += neighbor_params.kernel.template Weight<d>(l_particles.X(i) - boundary.X(b_nb), radius, KernelType::QUINTIC);
		}
		if (total_w > 0.) {
			l_particles.WM(i) = 1. / total_w;
		}
		else {
			l_particles.WM(i) = 0.;
		}
	}
	if (l_particles.Size() > 0) {
		real avg_wm = AuxFunc::Mean<real>(l_particles.WMRef());
#pragma omp parallel for
		for (int i = 0; i < boundary.Size(); i++) {
			boundary.WM(i) = avg_wm;
		}
	}
	else {
#pragma omp parallel for
		for (int i = 0; i < boundary.Size(); i++) {
			boundary.WM(i) = 0.;
		}
	}
}

// this is to make sure that the boundary particles "mirrors" the fluid inside
template<int d>
void ParticleALEFilm<d>::Configure_Boundary(void)
{
	avg_gamma = AuxFunc::Mean<real>(e_particles.GammaRef());
	avg_h = AuxFunc::Mean<real>(e_particles.HRef());
#pragma omp parallel for
	for (int i = 0; i < boundary.Size(); i++) {
		bool has_nbs = true;
		Array<int> e_nbs = e_particles.nbs_searcher->Find_Neighbors(boundary.X(i), neighbor_params.sph_radius);
		real min_h = avg_h;
		real min_gamma = avg_gamma;
		for (int k = 0; k < e_nbs.size(); k++) {
			int e_nb = e_nbs[k];
			if (k == 0 || e_particles.H(e_nb) < min_h) min_h = e_particles.H(e_nb);
			if (k == 0 || e_particles.Gamma(e_nb) < min_gamma) min_gamma = e_particles.Gamma(e_nb);
		}
		boundary.H(i) = min_h;
		boundary.Gamma(i) = min_gamma;
		boundary.Vol(i) = boundary.H(i) * boundary.SA(i);
		boundary.M(i) = boundary.Vol(i) * chem_params.rho_0;
		boundary.Soap(i) = boundary.Gamma(i) * boundary.SA(i);
	}
}

template<int d>
void ParticleALEFilm<d>::Update_E(void)
{
	e_particles.Update();
	if (analytical_boundary.Available()) e_particles.Update_Phis(analytical_boundary);
}

template<int d>
void ParticleALEFilm<d>::Update_L(void)
{
	l_particles.Update();
	if (analytical_boundary.Available()) l_particles.Update_Phis(analytical_boundary);
}

template<int d>
void ParticleALEFilm<d>::Update_E_Positions(const real dt)
{
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		e_particles.X(i) += e_particles.V(i) * dt;
	}
	if (camera_moving_along && !ext_acc_func) {
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) {
			e_particles.X(i) += init_COM - curr_COM;
		}
	}
	// enforce boundary condition
	if (analytical_boundary.Available()) e_particles.Correct_Position_With_Analytical_Boundary(dt, simulation_scale, analytical_boundary);
}

template<int d>
void ParticleALEFilm<d>::Update_L_Positions(const real dt)
{
#pragma omp parallel for
	for (int i = 0; i < l_particles.Size(); i++) {
		l_particles.X(i) += l_particles.V(i) * dt;
	}
	if (camera_moving_along && !ext_acc_func) {
#pragma omp parallel for
		for (int i = 0; i < l_particles.Size(); i++) {
			l_particles.X(i) += init_COM - curr_COM;
		}
	}

	// Project to E
	if (!flat_surface) Project_L(dt); // project L, unless the surface is perfectly 2d

	if (analytical_boundary.Available()) l_particles.Correct_Position_With_Analytical_Boundary(dt, simulation_scale, analytical_boundary);
}

template<int d>
void ParticleALEFilm<d>::Project_L(const real dt)
{
	double begin_time, end_time;
	if (verbose) {
		std::cout << "[L Project] begin" << std::endl;
		begin_time = omp_get_wtime();
	}
#pragma omp parallel for
	for (int i = 0; i < l_particles.Size(); i++) {
		//the number 7 here is carefully chosen. MLS<2,2> has 6 parameters, to avoid Runge phenomenon we force at least 7 points
		bool vb = false; //verbose?
		//int flg = e_particles.surface.Nearest_Geometry(l_particles.X(i), l_particles.E(i), 7, vb);
		int flg = e_particles.surface.Nearest_Geometry_KNN(l_particles.X(i), l_particles.E(i), 10, 7, vb);
		if (flg) to_delete_l[i] = 1;
		if (!bursting) {
			VectorD norm_component = AuxFunc::Component_Along(l_particles.V(i), l_particles.Normal(i));
			VectorD tang_component = l_particles.V(i) - norm_component;
			//revise velocity after they obtained new reference frame according to E's
			VectorD new_norm_component = norm_component.norm() * AuxFunc::Component_Along(norm_component, l_particles.Normal(i)).normalized();
			VectorD new_tang_component = tang_component.norm() * AuxFunc::Eliminate_Unit_Component(tang_component, l_particles.Normal(i)).normalized();
			l_particles.V(i) = new_norm_component + new_tang_component;
			l_particles.Mom(i) = l_particles.V(i) * l_particles.Vol(i) * chem_params.rho_0;
		}
	}
	if (verbose) {
		end_time = omp_get_wtime();
		std::cout << "[L Project] done. Time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
	}
}

template<int d>
void ParticleALEFilm<d>::L2E(bool use_affine, bool config_boundary) {
	double begin_time, end_time;
	if (verbose) {
		std::cout << "[L2E] begin" << std::endl;
		begin_time = omp_get_wtime();
	}
	if (config_boundary) Configure_Boundary();
	//prepare multiplier so that the sum of weight is 1 (necessary for MPM style transfer)
	//this is done for both the L and Boundary particles
	Update_WM(neighbor_params.interp_radius);

	//Zero out the arrays to be tranferred, mass, volume, affine momentum
	AuxFunc::Fill(e_particles.LVolRef(), 0.);
	AuxFunc::Fill(e_particles.SoapRef(), 0.);
	AuxFunc::Fill(e_particles.LMomRef(), VectorD::Zero());
	AuxFunc::Fill(e_particles.LAMomRef(), VectorD::Zero());

	// Compute the mass and volume from each L particle
	Distribute_L2E(e_particles.LVolRef(), l_particles.VolRef(), neighbor_params.interp_radius, KernelType::QUINTIC);
	Distribute_L2E(e_particles.SoapRef(), l_particles.SoapRef(), neighbor_params.interp_radius, KernelType::QUINTIC);

	// Compute the mass and volume of each B particle (Don't forget about the boundary particles)
	Distribute_B2E(e_particles.LVolRef(), boundary.VolRef(), neighbor_params.interp_radius, KernelType::QUINTIC);
	Distribute_B2E(e_particles.SoapRef(), boundary.SoapRef(), neighbor_params.interp_radius, KernelType::QUINTIC);


	// Compute the momentum from each L particle
	Distribute_L2E(e_particles.LMomRef(), l_particles.MomRef(), neighbor_params.interp_radius, KernelType::QUINTIC);
	// Compute the AFFINE momentum from each L particle
	AffineDistribute_L2E(e_particles.LAMomRef(), l_particles.AMomRef(), neighbor_params.interp_radius, KernelType::QUINTIC);

	if (consider_black_hole) Update_BH_Ratio();

	// Process what we have collected
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		e_particles.LM(i) = e_particles.LVol(i) * chem_params.rho_0;
		// First we deal with the fluid momentum (L Mom)

		// The total (fluid) momentum is the sum of momentum and affine momentum (same as APIC)
		if (use_affine) e_particles.LMom(i) = e_particles.LMom(i) + e_particles.LAMom(i);
		VectorD LMom_N = AuxFunc::Component_Along(e_particles.LMom(i), e_particles.Normal(i));
		VectorD LMom_T = e_particles.LMom(i) - LMom_N;
		// Then we deal with the sack momentum (E Mom), which is the sum of both momentums from E and L
		//e_particles.V(i) = (LMom_N + e_particles.V(i) * e_particles.base_mass) / (e_particles.LM(i) + e_particles.base_mass);
		// if there is no L particles nearby, then the collected tangential velocity is zero (beware of devision by 0) 
		if (e_particles.LM(i) > 0.)
			e_particles.LV_T(i) = LMom_T / e_particles.LM(i);
		else
			e_particles.LV_T(i) = VectorD::Zero();
		e_particles.LV(i) = AuxFunc::Component_Along(e_particles.V(i), e_particles.Normal(i)) + e_particles.LV_T(i); // the full velocity is the sum of tangential and normal components.
	
		// collect temperature
		bool has_nbs = true;
		e_particles.Temp(i) = l_particles.World_Avg_Value(e_particles.X(i), l_particles.TempRef(), neighbor_params.kernel, neighbor_params.interp_radius, KernelType::QUINTIC, has_nbs);
		if (!has_nbs) {
			e_particles.Temp(i) = chem_params.T_0;
		}
	}

	if (verbose) {
		end_time = omp_get_wtime();
		std::cout << "[L2E] done. Time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
	}
}

template<int d>
void ParticleALEFilm<d>::Update_C_Curvature(void) {
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		std::function<real(const int)> z_func_i = [&](const int idx)->real {return e_particles.X(idx).dot(e_particles.Normal(i)); };
		real KS = e_particles.Surface_Laplacian(i, z_func_i, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
		std::function<real(const int)> z_func_i_b = [&](const int idx)->real {return boundary.X(idx).dot(e_particles.Normal(i)); };
		real boundary_KS = boundary.Surface_Laplacian(e_particles.X(i), e_particles.E(i), e_particles.X(i).dot(e_particles.Normal(i)), z_func_i_b, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
		real total_KS = KS + boundary_KS;
		if (Eulerian_Near_Boundary(i, neighbor_params.sph_radius)) {
			real max_abs_KS = 2. / simulation_scale;
			if (abs(total_KS) > max_abs_KS) { total_KS *= max_abs_KS / abs(total_KS); }
		}
		e_particles.KS(i) = total_KS;
	}
	////smooth out KS
	for (int i = 0; i < numeric_params.XSPH_KS_Passes; i++) {
		XSPH_Smoothing<real>(0.999, e_particles.KSRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
	}

}

template<int d>
void ParticleALEFilm<d>::Update_E_Geometry(void) {
	double begin_time, end_time;
	if (verbose) {
		std::cout << "[E Geometry] begin" << std::endl;
		begin_time = omp_get_wtime();
	}
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		// Compute the number density of each E particle
		// Must remember to account for the boundary as well, so that the SA for particles near the boundary is not overly large
		e_particles.NDen(i) = e_particles.Surface_Sum_Value(e_particles.X(i), e_particles.E(i), ones_e, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC)
			+ boundary.Surface_Sum_Value(e_particles.X(i), e_particles.E(i), ones_b, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
		// SA = 1/NDen
		e_particles.SA(i) = 1. / e_particles.NDen(i);
		// H = LVol/SA
		e_particles.H(i) = e_particles.LVol(i) / e_particles.SA(i);
		// Clip to avoid singularity for things like 1/h
		real e_min_h = 0.1 * fluid_volume / fluid_area;
		e_particles.H(i) = std::max<real>(e_particles.H(i), e_min_h);
		// Den = LM/SA
		e_particles.Den(i) = e_particles.LM(i) / e_particles.SA(i);
		// Gamma = Soap/SA
		e_particles.Gamma(i) = e_particles.Soap(i) / e_particles.SA(i);
		// ST = st_water - RT * gamma
		e_particles.ST(i) = 7.275 * 0.01 - 8.314 * 298.15 * e_particles.Gamma(i);
		// The mass of a e_particle is the summation of its own mass and the mass of the fluid contained in it
		e_particles.M(i) = e_particles.LM(i) + e_particles.base_mass;
		// Same with the the volume
		e_particles.Vol(i) = e_particles.LVol(i) + e_particles.base_mass / e_particles.base_rho;
	}

	// if bursting, then near the rim, the estimation will be inaccurate
	// hence we set the Gamma to be the average value.
	if (bursting) {
		avg_gamma = AuxFunc::Mean<real>(e_particles.GammaRef());
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) {
			if (near_rim_confidence[i] > 0.)
				e_particles.Gamma(i) = avg_gamma;
		}
		std::function<bool(const int)> near_rim = [&](const int idx)->bool {return (near_rim_confidence[idx]>0.);};
		for (int pass = 0; pass < 3; pass++) {
			XSPH_Smoothing(0.999, e_particles.GammaRef(), neighbor_params.kernel, 1. * neighbor_params.sph_radius, KernelType::QUINTIC, near_rim);
		}
	}

	Update_C_Curvature();

	fluid_area = AuxFunc::Sum(e_particles.SARef());
	avg_l_sa = fluid_area / l_particles.Size();
	avg_e_sa = fluid_area / e_particles.Size();

	curr_COM = Compute_COM();
	curr_COM_velocity = Compute_COM_Velocity();
	//std::cout << "Current COM: \n" << curr_COM << std::endl;

	if (verbose) {
		end_time = omp_get_wtime();
		std::cout << "[E Geometry] done. Time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
	}
}


template<int d>
void ParticleALEFilm<d>::Update_L_Geometry(const real dt) {
	double begin_time, end_time;
	if (verbose) {
		std::cout << "[L Geometry] begin" << std::endl;
		begin_time = omp_get_wtime();
	}

//#pragma omp parallel for
//	for (int i = 0; i < l_particles.Size(); i++) {
//		if (temp_func) { l_particles.Temp(i) = temp_func(i); }
//		l_particles.H(i) = 1.e-9 * l_particles.Temp(i);
//		l_particles.SA(i) = avg_l_sa;
//	}
//	return;

	if (curr_frame == 0 && reinit_l_geometry) { //initialize in the beginning
#pragma omp parallel for
		for (int i = 0; i < l_particles.Size(); i++) {
			if (l_particles.Is_BH(i)>0.) {
				l_particles.H(i) = 0.;
				l_particles.SA(i) = avg_l_sa;
			}
			else {
				const VectorD& pos = l_particles.X(i);
				bool has_nbs = true;
				l_particles.H(i) = e_particles.World_Avg_Value(pos, e_particles.HRef(), neighbor_params.kernel, .5 * neighbor_params.interp_radius, KernelType::QUINTIC, has_nbs);
				l_particles.SA(i) = avg_l_sa;
			}
		}
	}
	else { //evolve afterwards
		// prepare for the decay
		real decay_factor = std::max(0., 1. - dt * l_geometry_sync_strength);
		if (l_geometry_sync) {
#pragma omp parallel for
			for (int i = 0; i < e_particles.Size(); i++) {
				// the interpolated value of l.h on E
				bool has_nbs = true;
				temp_eh[i] = l_particles.World_Avg_Value(e_particles.X(i), l_particles.HRef(), neighbor_params.kernel, neighbor_params.interp_radius, KernelType::QUINTIC, has_nbs);
			}
		}
		// done prepare for the decay
#pragma omp parallel for
		for (int i = 0; i < l_particles.Size(); i++) {
			// if bursting and a L particle is near the rim, then we skip
			bool skipping = false;
			if (bursting) {
				real tmp = e_particles.World_Sum_Value(l_particles.X(i), near_rim_confidence, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
				if (tmp > 0.) skipping = true;
			}
			// if is black hole particle
			if (l_particles.Is_BH(i)>0.) {
				l_particles.H(i) *= 0.;
				l_particles.SA(i) = avg_l_sa;
			}
			else if (skipping) {
				//pass
			}
			else {
				const VectorD& pos = l_particles.X(i);
				real dhdt = 0;

				// WAY 0
				if (l_geometry_mode == 0) {
					real dhdt2 = e_particles.Surface_Divergence(pos, l_particles.E(i), l_particles.V(i), temp_LV, neighbor_params.kernel, 2. * neighbor_params.sph_radius, KernelType::QUINTIC) +
						boundary.Surface_Divergence(pos, l_particles.E(i), l_particles.V(i), boundary.VRef(), neighbor_params.kernel, 2. * neighbor_params.sph_radius, KernelType::QUINTIC);
					dhdt2 *= -l_particles.H(i);
					dhdt = dhdt2;
				}
				// WAY 1
				else if (l_geometry_mode == 1) {
					if (consider_black_hole) {
						dhdt = e_particles.Surface_Density_Derivative(pos, l_particles.E(i), l_particles.V(i), e_particles.LVRef(), neighbor_params.kernel, 2. * neighbor_params.sph_radius, KernelType::QUINTIC) +
							boundary.Surface_Density_Derivative(pos, l_particles.E(i), l_particles.V(i), boundary.VRef(), neighbor_params.kernel, 2. * neighbor_params.sph_radius, KernelType::QUINTIC);
					}
					else {
						dhdt = e_particles.Surface_Density_Derivative(pos, l_particles.E(i), l_particles.V(i), e_particles.LVRef(), neighbor_params.kernel, 2. * neighbor_params.sph_radius, KernelType::QUINTIC) +
							boundary.Surface_Density_Derivative(pos, l_particles.E(i), l_particles.V(i), boundary.VRef(), neighbor_params.kernel, 2. * neighbor_params.sph_radius, KernelType::QUINTIC);
					}
				}


				l_particles.H(i) += dt * dhdt;
				// h has a minimal value
				real min_l_h = 20 * 1.E-9;
				real max_l_h = 8000 * 1.E-9;
				l_particles.H(i) = std::max<real>(min_l_h, l_particles.H(i));
				l_particles.H(i) = std::min<real>(max_l_h, l_particles.H(i));
				l_particles.SA(i) = avg_l_sa;

				if (l_geometry_sync) {
					//std::cout << "This part is getting callled!" << std::endl;
					bool has_nbs = true;
					temp_lh[i] = e_particles.World_Avg_Value(l_particles.X(i), temp_eh, neighbor_params.kernel, neighbor_params.interp_radius, KernelType::QUINTIC, has_nbs);
					eh[i] = e_particles.World_Avg_Value(l_particles.X(i), e_particles.HRef(), neighbor_params.kernel, neighbor_params.interp_radius, KernelType::QUINTIC, has_nbs);
					real err = eh[i] - temp_lh[i];
					real err_decayed = err * decay_factor;
					// decay the error
					l_particles.H(i) = (eh[i] - err_decayed) + (l_particles.H(i) - temp_lh[i]);
					l_particles.SA(i) = avg_l_sa;
				}
			}
		}
	}
	if (verbose) {
		end_time = omp_get_wtime();
		std::cout << "[L Geometry] done. Time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
	}
}

template<int d>
void ParticleALEFilm<d>::Compute_Rim_Normals(void)
{
	AuxFunc::Fill(is_rim_confidence, 0.);
	AuxFunc::Fill(near_rim_confidence, 0.);
	Array<real> temp_array; temp_array.resize(e_particles.Size());
	real threshold = 60;
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		VectorT norm_2d = e_particles.Surface_Gradient(i, ones_e, neighbor_params.kernel, 1.2 * neighbor_params.sph_radius, KernelType::QUINTIC);
		VectorD norm_3d; e_particles.surface.Unproject_To_World(norm_2d, e_particles.E(i), norm_3d);
		rim_normals[i] = norm_3d;
	}
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		if (rim_normals[i].norm() > threshold) {
			is_rim_confidence[i] = 1.;
			temp_array[i] = 1.;
		}
		else {
			is_rim_confidence[i] = 0.;
			temp_array[i] = 0.;
		}
	}
	for (int pass = 0; pass < (int)simulation_fineness; pass++) {
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) {
			near_rim_confidence[i] = e_particles.World_Sum_Value(e_particles.X(i), temp_array, neighbor_params.kernel, neighbor_params.interp_radius, KernelType::QUINTIC);
		}
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) {
			temp_array[i] = near_rim_confidence[i];
		}
	}
}

template<int d>
void ParticleALEFilm<d>::Apply_Rim_Surface_Tension(real dt)
{
	real threshold = 60;
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		if (is_rim_confidence[i] > 0.) {
			e_particles.F(i) = rim_normals[i] * 0.001;
			VectorD burst_acc = dt * rim_normals[i] * e_particles.ST(i) * 2. * (e_particles.Temp(i) / chem_params.T_0);
			e_particles.V(i) += burst_acc; // APPLY Burst acceleration to E.V only
		}
	}
}

template<int d>
void ParticleALEFilm<d>::Delete_Particles(void)
{
	real avg_NDen = AuxFunc::Mean<real>(e_particles.NDenRef());
	real delete_threshold = 150;
	// Determine which E particles should be deleted
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		if (e_particles.Size() < 0.02 * num_E_preburst) {
			to_delete_e[i] = 1;
		}
		real scores_loc = abs(e_particles.KS(i));
		if (scores_loc > delete_threshold
			|| e_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.sph_radius).size() < 7
			|| e_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), 3. * neighbor_params.sph_radius).size() < 200)// || e_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.sph_radius).size() < 10)
			to_delete_e[i] = 1;
	}

	// Add deleted E particles to 3d solver
	for (int i = 0; i < e_particles.Size(); i++) {
		if (to_delete_e[i]) {
			if (RandomFunc::Random_Real(0., 1.) < (1./(real)dilution_rate_3d))
			fluid_3d->Add_Particle(e_particles.X(i), .6 * e_particles.V(i) + -4. * rim_normals[i].normalized() * simulation_scale, 2.5684e-08);
		}
	}

	if (!fluid_3d->Is_Trivial()) {
		fluid_3d->Update_Max_Velocity(); //need to update max velocity once the new particles are transferred over
	}
	int original_e = e_particles.Size();
	int original_l = l_particles.Size();
	int remaining_e = e_particles.Delete_Elements_Safe(to_delete_e);
	if (remaining_e < original_e) {
		Update_E();
		if (!bursting) {
			bursting = true;
		}
	}

	// Determine which L particles should be deleted
#pragma omp parallel for
	for (int i = 0; i < l_particles.Size(); i++) {
		Array<int> nbs = e_particles.nbs_searcher->Find_Neighbors(l_particles.X(i), 0.5 * neighbor_params.sph_radius);
		if (nbs.size() < 1) to_delete_l[i] = 1;
	}
	int remaining_l = l_particles.Delete_Elements_Safe(to_delete_l);
	if (remaining_l < original_l) {
		Update_L();
	}

	if (verbose) {
		std::cout << "[Delete Particles] Deleted: " << original_e - remaining_e << " E particles" << std::endl;
		std::cout << "[Delete Particles] Deleted: " << original_l - remaining_l << " L particles" << std::endl;
	}

	if (remaining_l < original_l || remaining_e < original_e) {
		Reinit_Arrays();
	}
}

template<int d>
void ParticleALEFilm<d>::Update_Dynamics(const real dt) {
	// this is to compute raw velocity (before solving)
	Apply_Ext_Forces(dt);
	// raw velocity computed
	 
//	//Update tang velocity (LV_T)
	Update_Tang_Dynamics(dt);
//	//Update norm velocity
	Update_Norm_Dynamics(dt);
	//smooth out e_v
	XSPH_Smoothing_Tang_Norm(numeric_params.XSPH_V_Tang, numeric_params.XSPH_V_Norm, e_particles.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
	//additional smoothing of E.V
	for (int i = 0; i < numeric_params.XSPH_V_Passes; i++) {
		XSPH_Smoothing_Tang_Norm(0., numeric_params.XSPH_V_Passes_Strength, e_particles.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
	}
//	//Combine normal and tang velocities to form LV
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
		e_particles.LV(i) = AuxFunc::Component_Along(e_particles.V(i), e_particles.Normal(i)) + e_particles.LV_T(i);
	}
}


template<int d>
void ParticleALEFilm<d>::Update_Dynamics_E_Only(const real dt) {
	// this is to compute raw velocity (before solving)
	Apply_Ext_Forces(dt);

	Update_Norm_Dynamics(dt);
	//smooth out e_v
	XSPH_Smoothing_Tang_Norm(numeric_params.XSPH_V_Tang, numeric_params.XSPH_V_Norm, e_particles.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
	//additional smoothing of E.V
	for (int i = 0; i < numeric_params.XSPH_V_Passes; i++) {
		XSPH_Smoothing_Tang_Norm(0., numeric_params.XSPH_V_Passes_Strength, e_particles.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
	}
}

template<int d>
void ParticleALEFilm<d>::Update_Burst(const real dt) {
	Apply_Ext_Forces(dt);
	Update_Tang_Dynamics(dt);
	Update_Norm_Dynamics(dt);
	Apply_Rim_Surface_Tension(dt);
	XSPH_Smoothing_Tang_Norm(0.999, 0., e_particles.VRef(), neighbor_params.kernel, 1. * neighbor_params.sph_radius, KernelType::QUINTIC);
	//Combine normal and tang velocities to form LV
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
		e_particles.LV(i) = AuxFunc::Component_Along(e_particles.V(i), e_particles.Normal(i)) +e_particles.LV_T(i);
	}
}

template<int d>
void ParticleALEFilm<d>::Update_Norm_Dynamics(const real dt) {
	AuxFunc::Fill(e_particles.FRef(), VectorD::Zero());

	Update_Capillary_Forces();
	if (enclosed > 0 && !bursting) Update_Enclosed_Forces();

	if (!(enclosed > 0) && !camera_moving_along && euler_feels_gravity) {
		faux_g = g; // record this so that Even_Out_E knows
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
			// this is the physical force
			VectorD gravity = e_particles.M(i) * g;
			e_particles.F(i) += gravity;
		}
	}
	else {
		faux_g = VectorD::Zero();
	}

	// account for wind
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
		VectorD vel_air = (wind_func != nullptr) ? wind_func(i) : e_particles.V(i);
		VectorD vel_diff = vel_air - e_particles.V(i);
		//impact area
		//note that we assume that wind is applied ONLY IF the normal is opposing the wind direction (relative velocity)
		VectorD direction = vel_diff.normalized();
		real impact_area = avg_e_sa * std::max<real>(0., (-e_particles.Normal(i)).dot(direction));
		e_particles.F(i) += chem_params.drag_C * impact_area * vel_diff;
	}

#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
		e_particles.V(i) += e_particles.F(i)/e_particles.M(i) * dt;
	}

	if (camera_moving_along) {
		// if the camera is moving along with the center of the bubble
		// then the E dynamics is ACTUALLY falling, but is not reflected in simulation
		// hence the L velocity (LV_T) which is in the same reference frame as E, need to be compensated.
		// if is FALLING SLOWER, then it is acutually RISING in the falling reference frame.
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
			e_particles.LV_T(i) += dt * -AuxFunc::Eliminate_Unit_Component(g, e_particles.Normal(i));
		}
	}
}

template<int d>
void ParticleALEFilm<d>::E2L(void) {
	double begin_time, end_time;
	if (verbose) {
		std::cout << "[E2L] begin" << std::endl;
		begin_time = omp_get_wtime();
	}
	// Interp velocity for each L particle in the APIC way
#pragma omp parallel for
	for (int i = 0; i < l_particles.Size(); i++) {
		Array<int> nbs = e_particles.nbs_searcher->Find_Neighbors(l_particles.X(i), neighbor_params.interp_radius);
		l_particles.DM(i) *= 0.;
		l_particles.V(i) *= 0.;
		l_particles.BM(i) *= 0;
		for (int k = 0; k < nbs.size(); k++) {
			int e_nb = nbs[k];
			VectorD wr_pj = e_particles.X(e_nb) - l_particles.X(i);
			VectorT lr_pj = PointSet<d>::Rotate_To_TPlane(wr_pj, l_particles.E(i));
			real weight = l_particles.WM(i) * neighbor_params.kernel.template Weight<d>(wr_pj, neighbor_params.interp_radius, KernelType::QUINTIC);
			VectorT tang_LV = PointSet<d>::Project_To_TPlane(e_particles.LV(e_nb), l_particles.E(i));
			l_particles.V(i) += weight * e_particles.LV(e_nb);
			l_particles.BM(i) += weight * tang_LV * (lr_pj).transpose();
			l_particles.DM(i) += weight * (lr_pj) * (lr_pj).transpose();
		}
		// need to consider the boundary particles for the DM (inertia like tensor)
		nbs = boundary.nbs_searcher->Find_Neighbors(l_particles.X(i), neighbor_params.interp_radius);
		for (int k = 0; k < nbs.size(); k++) {
			int b_nb = nbs[k];
			VectorD wr_pj = boundary.X(b_nb) - l_particles.X(i);
			VectorT lr_pj = PointSet<d>::Rotate_To_TPlane(wr_pj, l_particles.E(i));
			real weight = l_particles.WM(i) * neighbor_params.kernel.template Weight<d>(wr_pj, neighbor_params.interp_radius, KernelType::QUINTIC);
			l_particles.DM(i) += weight * (lr_pj) * (lr_pj).transpose();
		}
	}

	// put together Momentum and Affine Momentum
#pragma omp parallel for
	for (int i = 0; i < l_particles.Size(); i++) {
		l_particles.Mom(i) = l_particles.V(i) * l_particles.Vol(i) * chem_params.rho_0;
		Matrix<real, d - 1> DM_inverse = l_particles.DM(i).inverse();
		bool legit = true;
		for (int axis = 0; axis < d-1; axis++) {
			for (int axis2 = 0; axis2 < d-1; axis2++) {
				real val = DM_inverse(axis, axis2);
				if (std::isnan(val) || !std::isfinite(val)) {
					legit = false;
					break;
				}
			}
			if (!legit) break;
		}
		if (!legit) l_particles.AMom(i) = 0. * l_particles.BM(i);
		else l_particles.AMom(i) = l_particles.BM(i) * DM_inverse * l_particles.Vol(i) * chem_params.rho_0;
	}

	if (verbose) {
		end_time = omp_get_wtime();
		std::cout << "[E2L] done, time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
	}
}




template<int d>
void ParticleALEFilm<d>::Advance_E(const real dt, const int iter) {
	// Update position of each E particle
	Update_E_Positions(dt);
	// Update E NB searcher
	Update_E();
	// Similar to pressure projection, but done after the new surface is formed
	if (!bursting) Even_Out_E(dt, false, max_even_out_e_steps); //if bursting or fixed_surface then we don't even out E
}

template<int d>
void ParticleALEFilm<d>::Advance_L(const real dt) {
	// Update position of each L particle
	Update_L_Positions(dt);
	// Update L NB searcher
	Update_L();
	// Update Temperature
	if (temp_func != nullptr) {
#pragma omp parallel for
		for (int i = 0; i < l_particles.Size(); i++) {
			l_particles.Temp(i) += 1./chem_params.heat_capacity_param * dt * (temp_func(i) - l_particles.Temp(i));
		}
	}
}


template<int d>
void ParticleALEFilm<d>::Update_Capillary_Forces(void)
{
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		real pressure_difference = e_particles.ST(i) * e_particles.KS(i);//capillary pressure difference
		VectorD force = pressure_difference * avg_e_sa * e_particles.Normal(i);
		e_particles.F(i) += 100 * force;
	}
}

// Make L (Material Mass/Volume) Evenly Distributed
template<int d>
void ParticleALEFilm<d>::Even_Out_L(const int max_steps)
{
	double begin_time, end_time;
	if (verbose) {
		std::cout << "[L Redistribute] begin" << std::endl;
		begin_time = omp_get_wtime();
	}
	real avg_err = -1.;
	int num_iter = 0;
	for (int k = 0; k < max_steps; k++) {
		//double begin_time = omp_get_wtime();
		real dt = (cfl * neighbor_params.e_dx) / std::max<real>(cfl * neighbor_params.e_dx * frame_rate, max_vel_l);
		L2E(); // transfer momentum, mass, etc
		Update_E_Geometry();
		//Update tang velocity
		Update_Tang_Dynamics(dt, true);
		//Combine normal and tang velocities to form LV
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
			e_particles.LV(i) = e_particles.LV_T(i);
		}
		Decay_Velocity(e_particles.LVRef(), dt, 10., 10.);
		E2L(); // transfer velocity
		Advance_L(dt);
		real new_avg_err = AuxFunc::Max<real>(e_particles.GammaRef()) - AuxFunc::Min<real>(e_particles.GammaRef());
		num_iter++;

		Update_L_Max_Velocity();

		if (verbose)
			std::cout << "[L Redistribute] At iter: " << num_iter << ", height error: " << std::setprecision(2) << new_avg_err << std::endl;
		if (avg_err >= 0.) {
			real err_diff = abs((new_avg_err - avg_err)) / avg_err;
			if (err_diff < 1.e-8) {
				if (verbose)
					std::cout << "[L Redistribute] Converged." << std::endl;
				break;
			}
		}
		avg_err = new_avg_err;
	}
	Update_L_Geometry(0.);
	if (verbose) {
		end_time = omp_get_wtime();
		std::cout << "[L Redistribute] done. Num iters: " << num_iter << ", time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
	}
}

// Make L (Material Mass/Volume) Evenly Distributed
// if "repacking" is set to false, then we don't apply the e_particles.V(i), since we assume that it's in the middle of simulation
// and the velocity is applied in advance_E
// if repacking is true, then we apply the tangential component of e_particles.V(i), so that the velocity is not forgotten each step
template<int d>
void ParticleALEFilm<d>::Even_Out_E(const real _dt, const bool repacking, const int _max_steps)
{
	double begin_time, end_time;
	if (verbose) {
		std::cout << "[E redistribute] begin" << std::endl;
		begin_time = omp_get_wtime();
	}
	real avg_err = -1.; // initial error is 0

	int max_steps = _max_steps;
	if (bursting) max_steps = 1;

	int num_iter = 0;
	for (int k = 0; k < max_steps; k++) {
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) {
			e_particles.FauxV(i) = VectorD::Zero();
		}
		real dt = _dt;
		if (repacking) { // if repacking, then the "dt" is solely dependent upon the max e velocity
			dt = (cfl * neighbor_params.e_dx) / std::max<real>(cfl * neighbor_params.e_dx * frame_rate, max_vel_e);
			// this below basically does what Advance_E does before Even_out_E
			// but with velocity only acting along the normal direction
#pragma omp parallel for
			for (int i = 0; i < e_particles.Size(); i++) {
				//we will perform an extra projection of velocity, since we are not technically inside simulation yet...!
				e_particles.V(i) = AuxFunc::Eliminate_Unit_Component(e_particles.V(i), e_particles.Normal(i));
				e_particles.X(i) += e_particles.V(i) * dt;
			}
			// enforce boundary condition
			if (analytical_boundary.Available()) e_particles.Correct_Position_With_Analytical_Boundary(dt, simulation_scale, analytical_boundary);
			Update_E();
		}

		real elas_coeff = 20. / e_evenness; //may want to modify 0.8 here
		std::cout << "elas_coeff without clip: " << elas_coeff << std::endl;
		real max_elas_coeff = 100./AuxFunc::Mean(e_particles.NDenRef());
		elas_coeff = std::min<real>(elas_coeff, max_elas_coeff); // may want to modify 100 here
		std::cout << "max allowed elas_coeff: " << max_elas_coeff << std::endl;
		std::cout << "elas_coeff: " << elas_coeff << std::endl;
		//elas_coeff *= 1. / 4.;
		real solver_err = Update_Faux_Tang_Dynamics(dt, elas_coeff);
		//real solver_err = 0;

		Update_E_Max_Velocity();
		if (repacking) Decay_Velocity(e_particles.VRef(), dt, 10., 10.);

		// <-----Recompute NDen----->
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) {
			// Compute the number density of each E particle
			// Must remember to account for the boundary as well, so that the SA for particles near the boundary is not overly large
			e_particles.NDen(i) = e_particles.Surface_Sum_Value(e_particles.X(i), e_particles.E(i), ones_e, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC)
				+ boundary.Surface_Sum_Value(e_particles.X(i), e_particles.E(i), ones_b, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
			// SA = 1/NDen
			e_particles.SA(i) = 1. / e_particles.NDen(i);
		}
		real new_avg_err = 1./e_particles.Size() * (AuxFunc::Max<real>(e_particles.NDenRef()) - AuxFunc::Min<real>(e_particles.NDenRef()));

		num_iter++;
		curr_nden_err = new_avg_err;
		if (verbose || 1==1) 
			std::cout << "[E Redistribute] At iter: " << num_iter << ", nden error: " << std::setprecision(2) << curr_nden_err << std::endl;
		if (avg_err >= 0.) {
			real err_diff = abs((new_avg_err - avg_err)) / avg_err;
			if (err_diff < 1.e-3) {
				if (verbose)
					std::cout << "[E Redistribute] Converged." << std::endl;
				//break;
			}
		}
		avg_err = new_avg_err;
		std::cout << "[E redistribute] finished iter: " << k << " error is: " << avg_err << std::endl;
		if (!repacking && avg_err <= 3.) break;
	}
	if (verbose) {
		end_time = omp_get_wtime();
		std::cout << "[E Redistribute] done. Num iters: " << num_iter << ", time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
	}
}

template<int d>
void ParticleALEFilm<d>::Initial_Packing(const int e_steps, const int l_steps)
{
	//Even out E and L
	Even_Out_E(1., true, e_steps);
	Even_Out_L(l_steps);
	//Clear velocities used in the process
	AuxFunc::Fill(e_particles.VRef(), VectorD::Zero());
	AuxFunc::Fill(l_particles.VRef(), VectorD::Zero());
	AuxFunc::Fill(l_particles.MomRef(), VectorD::Zero());
	AuxFunc::Fill(l_particles.AMomRef(), MatrixT::Zero());
	AuxFunc::Fill(e_particles.LVRef(), VectorD::Zero());
}

template<int d>
real ParticleALEFilm<d>::Compute_Enclosed_Volume(void) {
	VectorD origin = VectorD::Zero();
	real vol = 0.;
	for (int i = 0; i < e_particles.Size(); i++) {
		real height = fabs(e_particles.Normal(i).dot(e_particles.X(i) - origin)); // volume of the skewed cone with the origin
		real cone_vol = real(1. / d) * e_particles.SA(i) * height;
		if ((e_particles.X(i)-origin).dot(e_particles.Normal(i))>0.) {
			vol += cone_vol;
		}
		else {
			vol -= cone_vol;
		}
	}
	//std::cout << "Current enclosed volume: \n" << vol << std::endl;
	return vol;
}

template<int d>
void ParticleALEFilm<d>::Update_Enclosed_Forces(void)
{
	//compute enclosed volume
	real vol = Compute_Enclosed_Volume();
	if (verbose) std::cout << "[Enclosed] curr enclosed volume: " << vol << std::endl;
	real pressure = (vol > 1e-8) ? (enclosed_amount)/vol : outer_pressure;
	if (verbose) std::cout << "[Enclosed] curr enclosed pressure: " << pressure << std::endl;
	VectorD center = AuxFunc::Mean<VectorD>(e_particles.XRef());
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) {
		real sgn = 1;
		if ((e_particles.X(i)-center).dot(e_particles.Normal(i)) < 0) sgn = -1;
		real pressure_difference = sgn * (pressure - outer_pressure);
		VectorD force = pressure_difference * avg_e_sa * e_particles.Normal(i);
		e_particles.F(i) += force;
	}
}

template<int d>
real ParticleALEFilm<d>::Update_Faux_Tang_Dynamics(const real dt, const real elasticity) {
	if (dt < 1.e-8) return -1.; // if dt is 0. then return. Otherwise would be numerically troublesome

	Array<real> P; P.resize(e_particles.Size()); //temporary P
	real omega = 0.5; //Jacobi Relaxation
	int max_iter = 30;
	int max_max_iter = 100;
	Array<VectorT> temp; temp.resize(e_particles.Size()); //temp = [sum over j] {SA_j * 1/2 * grad_ij} (for computing diagonal terms)
	Array<real> a_ii; a_ii.resize(e_particles.Size()); //diagonal elements of matrix A
	Array<real> s_i; s_i.resize(e_particles.Size()); //source term
	Array<VectorD> grad_p; grad_p.resize(e_particles.Size()); //grad_p = Surface_Gradient_Symmetric(P)
	Array<real> LHS; LHS.resize(e_particles.Size()); //Left hand side
	Array<real> err; err.resize(e_particles.Size()); //error
	real avg_err = 0.;

	real Ma = elasticity;

#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
		P[i] = e_particles.NDen(i); // Here we are just using P as a placeholder for the revised/updated NDen
	}

#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
		//Compute temp
		//temp[i] = e_particles.Surface_Sum_Grad(e_particles.X(i), e_particles.E(i), e_particles.SARef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
		//temp[i] += boundary.Surface_Sum_Grad(e_particles.X(i), e_particles.E(i), boundary.SARef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
		//temp[i] *= -1;
		temp[i] = VectorT::Zero();
	}

#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
		Array<int> e_nbs = e_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.sph_radius);
		Array<int> b_nbs = boundary.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.sph_radius);
		a_ii[i] = 0.;
		//Compute a_ii
		real a_ii_0 = 0.; // zero order term
		real a_ii_1 = 0.; // first order term
		real a_ii_2 = 0.; // second order term
		// zero order
		a_ii_0 = -1. / (dt * e_particles.NDen(i));

		for (int k = 0; k < e_nbs.size(); k++) {
			int j = e_nbs[k];
			VectorD wr_ji = e_particles.X(i) - e_particles.X(j);
			VectorT sr_ji = PointSet<d>::Rotate_To_TPlane(wr_ji, e_particles.E(i));
			VectorT grad_W = neighbor_params.kernel.template Grad<d-1>(sr_ji, neighbor_params.sph_radius, KernelType::QUINTIC);
			a_ii_2 += e_particles.SA(j) * (temp[i] + e_particles.SA(i) * grad_W).dot(grad_W);
		}
		for (int k = 0; k < b_nbs.size(); k++) {
			int j = b_nbs[k];
			VectorD wr_ji = e_particles.X(i) - boundary.X(j);
			VectorT sr_ji = PointSet<d>::Rotate_To_TPlane(wr_ji, e_particles.E(i));
			VectorT grad_W = neighbor_params.kernel.template Grad<d - 1>(sr_ji, neighbor_params.sph_radius, KernelType::QUINTIC);
			a_ii_2 += boundary.SA(j) * (temp[i] + e_particles.SA(i) * grad_W).dot(grad_W);
		}

		//Compute Diagonal Elements a_ii
		a_ii_2 *= -1.;
		a_ii_2 *= dt * Ma;
		a_ii[i] = a_ii_0 + a_ii_1 + a_ii_2;
		//Compute source term s_i
		real velocity_divergence = e_particles.Surface_Divergence(e_particles.X(i), e_particles.E(i), e_particles.V(i), e_particles.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC) +
			boundary.Surface_Divergence(e_particles.X(i), e_particles.E(i), e_particles.V(i), boundary.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
		real gravity_divergence = boundary.Surface_Divergence(e_particles.X(i), e_particles.E(i), faux_g, boundary.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
		s_i[i] = velocity_divergence - 1. / dt + dt * gravity_divergence;
	}

	// Iterate
	int num_iters = 0;
	if (bursting) max_max_iter = 0;
	for (int l = 0; l < max_max_iter; l++) {
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
			VectorT grad2d = VectorT::Zero();
			//Compute grad p
			grad2d = e_particles.Surface_Gradient_Difference(i, P, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
			VectorD grad3d;
			e_particles.surface.Unproject_To_World(grad2d, e_particles.E(i), grad3d);
			grad_p[i] = grad3d;
		}
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
			//Compute LHS
			// zero order
			LHS[i] = -1. / (dt * e_particles.NDen(i)) * P[i];
			// second order
			real laplace = e_particles.Surface_Divergence(e_particles.X(i), e_particles.E(i), grad_p[i], grad_p, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC) +
				boundary.Surface_Divergence(e_particles.X(i), e_particles.E(i), grad_p[i], boundary.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
			LHS[i] += dt * Ma * laplace;
			//Update pressure
			//e_particles.P(i) = std::max<real>(0, e_particles.P(i) + omega / a_ii[i] * (s_i[i] - LHS[i]));
			P[i] = std::max<real>(0., P[i] + omega / a_ii[i] * (s_i[i] - LHS[i]));
			err[i] = abs(LHS[i] - s_i[i]);
			//temp_scalar[i] = 1.e-8 * P[i];
		}

		num_iters++;
		real new_avg_err = AuxFunc::Mean<real>(err);
		real err_change = avg_err > 0. ? (avg_err - new_avg_err) / avg_err : 0.;
		if (verbose)
			std::cout << "[E Redistribute Solver] at iter: " << num_iters << " error: " << new_avg_err << " improvement: " << err_change * 100 << "%" << std::endl;
		if (avg_err > 0. && (abs(err_change) < 1.e-4 || err_change < 0.)) {
			if (verbose)
				std::cout << "[E Redistribute Solver] Converged." << std::endl;
			break;
		}
		if (num_iters > max_iter && new_avg_err < 0.1) {
			break;
		}
		avg_err = new_avg_err;
	}

	e_evenness = ArrayFunc::Largest_Norm(grad_p);

	//Advance the velocity
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
		//Compute pressure acceleration ap (ap3d)
		VectorT ap_2d = VectorT::Zero();
		VectorD ap_3d;
		ap_2d = -Ma * e_particles.Surface_Gradient_Difference(i, P, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
		e_particles.surface.Unproject_To_World(ap_2d, e_particles.E(i), ap_3d);
		//temp_vector[i] = 0.01 * ap_3d;
		VectorD vp_3d = dt * ap_3d;
		real max_displacement = 0.2 * neighbor_params.e_dx;
		if (vp_3d.norm() > max_displacement / dt) vp_3d = vp_3d.normalized() * max_displacement / dt;
		if (is_rim_confidence[i] > 0.) { vp_3d *= 0.; } // if is boundary, then don't bother evening it out!!!
		e_particles.FauxV(i) = vp_3d;
	}

	XSPH_Smoothing<VectorD>(0.99, e_particles.FauxVRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);

	if (enforce_momentum_conserve) {
		VectorD mean_FauxV = AuxFunc::Mean(e_particles.FauxVRef());
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
			e_particles.FauxV(i) -= 0.5 * mean_FauxV;
		}
	}

#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
		e_particles.X(i) += e_particles.FauxV(i) * dt;
		//e_particles.V(i) += e_particles.FauxV(i);
	}

	if (analytical_boundary.Available()) e_particles.Correct_Position_With_Analytical_Boundary(dt, simulation_scale, analytical_boundary);
	Update_E();

	return avg_err;
}

// MAIN TANGENTIAL ACCELERATION COMPUTATION, difference form
template<int d>
void ParticleALEFilm<d>::Update_Tang_Dynamics(const real dt, const bool internal_force_only) {// Solve for Gamma is a fashion similar to Huang 2020
	double begin_time, end_time;
	if (verbose) {
		std::cout << "[Tang Dynamics] begin" << std::endl;
		begin_time = omp_get_wtime();
	}
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
		// Initialize P to Gamma
		e_particles.P(i) = e_particles.Gamma(i); // Here we are just using P as a placeholder for the revised/updated Gamma
	}
	if (dt < 1.e-8) return; // if dt is 0. then return.

	real Ma = 1000 * 2. * 0.5 * chem_params.R * chem_params.T_0 / chem_params.rho_0; // Marangoni param
	real Dr = chem_params.drag_C / chem_params.rho_0; // Drag param

	//Begin Solve
	//real omega = 0.2; //Jacobi Relaxation
	real omega = 0.2; //try a more conservative omega??
	//int max_iter = 2. * simulation_scale / neighbor_params.e_dx; //max_iter
	int max_iter = 30.;
	int max_max_iter = 100.;
	Array<VectorT> temp; temp.resize(e_particles.Size()); //temp = [sum over j] {SA_j * 1/2 * grad_ij} (for computing diagonal terms)
	Array<real> Bs; Bs.resize(e_particles.Size()); //the bouyancy variable B
	Array<VectorT> grad_inv_h; grad_inv_h.resize(e_particles.Size()); //grad of 1/h
	Array<VectorD> vel_diff; vel_diff.resize(e_particles.Size()); //vel_diff
	Array<real> a_ii; a_ii.resize(e_particles.Size()); //diagonal elements of matrix A
	Array<real> s_i; s_i.resize(e_particles.Size()); //source term
	Array<VectorD> grad_p; grad_p.resize(e_particles.Size()); //grad_p = Surface_Gradient_Symmetric(P)
	Array<real> LHS; LHS.resize(e_particles.Size()); //Left hand side
	Array<real> err; err.resize(e_particles.Size()); //error

	real avg_err = 0.;

#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
		//Compute temp
		//temp[i] = e_particles.Surface_Sum_Grad(e_particles.X(i), e_particles.E(i), e_particles.SARef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
		//temp[i] += boundary.Surface_Sum_Grad(e_particles.X(i), e_particles.E(i), boundary.SARef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
		//temp[i] *= -1;
		temp[i] = VectorT::Zero();
		//Compute B 
		//Bs[i] = 1.;
		Bs[i] = 1 - chem_params.expansion * (e_particles.Temp(i) - chem_params.T_0);
		if (i == 0) std::cout << "l particle temp? " << l_particles.Temp(0) << std::endl;
		if (i == 0) std::cout << "BS? " << Bs[i] << std::endl;
		if (i == 0) std::cout << "e_particles.Temp(i)? " << e_particles.Temp(i) << std::endl;
		if (i == 0) std::cout << "chem_params.T_0? " << chem_params.T_0 << std::endl;
		if (i == 0) std::cout << "chem_params.expansion? " << chem_params.expansion << std::endl;
		//Compute Velocity difference
		VectorD vel_air = (wind_func != nullptr) ? wind_func(i) : e_particles.LV(i);
		VectorD diff = vel_air - e_particles.LV(i);
		vel_diff[i] = AuxFunc::Eliminate_Unit_Component(diff, e_particles.Normal(i));
	}

#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
		Array<int> e_nbs = e_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.sph_radius);
		Array<int> b_nbs = boundary.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.sph_radius);
		a_ii[i] = 0.;
		//Compute a_ii
		real a_ii_0 = 0.; // zero order term
		real a_ii_1 = 0.; // first order term
		real a_ii_2 = 0.; // second order term
		// zero order
		a_ii_0 = -1. / (dt * e_particles.Gamma(i));
		// first order // boundary is already taken care of in "temp"
		std::function<real(const int)> inv_h = [&](const int idx)->real {return 1. / e_particles.H(idx); };
		grad_inv_h[i] = e_particles.Surface_Gradient_Difference(i, inv_h, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC); //No need to compensate for the difference form
		VectorT grad_B = e_particles.Surface_Gradient_Difference(i, Bs, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
		a_ii_1 = temp[i].dot(dt * Ma * grad_inv_h[i]);

		for (int k = 0; k < e_nbs.size(); k++) {
			int j = e_nbs[k];
			VectorD wr_ji = e_particles.X(i) - e_particles.X(j);
			VectorT sr_ji = PointSet<d>::Rotate_To_TPlane(wr_ji, e_particles.E(i));
			VectorT grad_W = neighbor_params.kernel.template Grad<d - 1>(sr_ji, neighbor_params.sph_radius, KernelType::QUINTIC);
			a_ii_2 += e_particles.SA(j) * (temp[i] + e_particles.SA(i) * grad_W).dot(grad_W);
		}
		for (int k = 0; k < b_nbs.size(); k++) {
			int j = b_nbs[k];
			VectorD wr_ji = e_particles.X(i) - boundary.X(j);
			VectorT sr_ji = PointSet<d>::Rotate_To_TPlane(wr_ji, e_particles.E(i));
			VectorT grad_W = neighbor_params.kernel.template Grad<d - 1>(sr_ji, neighbor_params.sph_radius, KernelType::QUINTIC);
			a_ii_2 += boundary.SA(j) * (temp[i] + e_particles.SA(i) * grad_W).dot(grad_W);
		}

		//Compute Diagonal Elements a_ii
		a_ii_2 *= -1.;
		a_ii_2 *= dt * Ma * 1. / e_particles.H(i);
		a_ii[i] = a_ii_0 + a_ii_1 + a_ii_2;
		//Compute source term s_i
		real velocity_divergence = e_particles.Surface_Divergence(e_particles.X(i), e_particles.E(i), e_particles.LV(i), e_particles.LVRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC) +
			boundary.Surface_Divergence(e_particles.X(i), e_particles.E(i), e_particles.LV(i), boundary.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
		real gravity_divergence = Bs[i] * boundary.Surface_Divergence(e_particles.X(i), e_particles.E(i), g, boundary.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC) +
			grad_B.dot(PointSet<d>::Project_To_TPlane(g, e_particles.E(i)));
		real drag_divergence = e_particles.Surface_Divergence(e_particles.X(i), e_particles.E(i), vel_diff[i], vel_diff, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC) +
			boundary.Surface_Divergence(e_particles.X(i), e_particles.E(i), vel_diff[i], boundary.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
		drag_divergence = grad_inv_h[i].dot(PointSet<d>::Project_To_TPlane(vel_diff[i], e_particles.E(i))) + 1. / e_particles.H(i) * drag_divergence;
		s_i[i] = velocity_divergence - 1. / dt;
		if (!internal_force_only) {
			s_i[i] += dt * gravity_divergence + dt * Dr * drag_divergence;
		}
	}

	// Iterate
	int num_iters = 0;
	if (bursting) max_max_iter = 0;
	for (int l = 0; l < max_max_iter; l++) {
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
			VectorT grad2d = VectorT::Zero();
			//Compute grad p
			grad2d = e_particles.Surface_Gradient_Difference(i, e_particles.PRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
			VectorD grad3d;
			e_particles.surface.Unproject_To_World(grad2d, e_particles.E(i), grad3d);
			grad_p[i] = grad3d;
		}
#pragma omp parallel for
		for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
			//Compute LHS
			// zero order
			LHS[i] = -1. / (dt * e_particles.Gamma(i)) * e_particles.P(i);
			// first order
			LHS[i] += (dt * Ma * grad_inv_h[i]).dot(PointSet<d>::Project_To_TPlane(grad_p[i], e_particles.E(i)));
			// second order
			real laplace = e_particles.Surface_Divergence(e_particles.X(i), e_particles.E(i), grad_p[i], grad_p, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC) +
				boundary.Surface_Divergence(e_particles.X(i), e_particles.E(i), grad_p[i], boundary.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
			LHS[i] += dt * Ma * 1. / e_particles.H(i) * laplace;
			//Update pressure
			//e_particles.P(i) = std::max<real>(0, e_particles.P(i) + omega / a_ii[i] * (s_i[i] - LHS[i]));
			e_particles.P(i) = e_particles.P(i) + omega / a_ii[i] * (s_i[i] - LHS[i]);
			err[i] = abs(LHS[i] - s_i[i]);
		}

		num_iters++;
		real new_avg_err = AuxFunc::Mean<real>(err);
		real err_change = avg_err > 0. ? (avg_err - new_avg_err) / avg_err : 0.;
		if (verbose)
			std::cout << "[Tang Dynamics Solver] at iter: " << num_iters << " error: " << new_avg_err << " improvement: " << err_change * 100 << "%" << std::endl;
		if (avg_err > 0. && (abs(err_change) < 1.e-4 || err_change < 0.)) {
			if (verbose)
				std::cout << "[Tang Dynamics] Converged." << std::endl;
			break;
		}
		if (num_iters > max_iter && new_avg_err < 0.1) {
			break;
		}
		avg_err = new_avg_err;
	}
	if (verbose) {
		std::cout << "[Tang Dynamics] done. Num iters: " << num_iters << " final error: " << avg_err << std::endl;
	}

	//Advance the velocity
#pragma omp parallel for
	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
		//Compute pressure acceleration ap (ap3d)
		VectorT ap_2d = VectorT::Zero();
		VectorD ap_3d;

		ap_2d = -Ma / e_particles.H(i) * (e_particles.Surface_Gradient_Difference(i, e_particles.PRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC));

		VectorD increment = VectorD::Zero();

		e_particles.surface.Unproject_To_World(ap_2d, e_particles.E(i), ap_3d);
		increment += dt * ap_3d;
		if (!internal_force_only) {
			increment += dt * Dr / e_particles.H(i) * AuxFunc::Eliminate_Unit_Component(vel_diff[i], e_particles.Normal(i)) + dt * Bs[i] * AuxFunc::Eliminate_Unit_Component(g, e_particles.Normal(i));
		}
		if (is_rim_confidence[i] > 0.) { increment *= 0; }
		e_particles.LV_T(i) += increment;
	}

	if (verbose) {
		end_time = omp_get_wtime();
		std::cout << "[Tang Dynamics] done. Time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
	}
}

template class ParticleALEFilm<2>;
template class ParticleALEFilm<3>;
