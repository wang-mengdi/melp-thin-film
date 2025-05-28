//#include "SPHBubbleFilm.h"
//#include "RandomNumber.h"
//
//// apply black hole forces
//template<int d>
//void SPHBubbleFilm<d>::Evaporate(const real dt)
//{
//#pragma omp parallel for
//	for (int i = 0; i < l_particles.Size(); i++) {
//		real rate = .1;
//		real factor = 1 - dt * .1;
//		if (temp_func && chem_params.T_0 > 0.) factor *= l_particles.Temp(i) / chem_params.T_0;
//		l_particles.Vol(i) *= factor;
//		l_particles.H(i) *= factor;
//	}
//}
//
//// apply black hole forces
//template<int d>
//void SPHBubbleFilm<d>::Apply_Ext_Forces(const real dt)
//{
//	if (!ext_acc_func) return;
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) {
//		VectorD ext_acc = ext_acc_func(i);
//		ext_acc *= dt;
//		e_particles.V(i) += ext_acc;
//		e_particles.LV(i) += ext_acc;
//		e_particles.LV_T(i) += AuxFunc::Eliminate_Unit_Component(ext_acc, e_particles.Normal(i));
//	}
//}
//
//template<int d>
//Vector<real, d> SPHBubbleFilm<d>::Compute_COM(void)
//{
//	Array<VectorD> scaled_positions; scaled_positions.resize(e_particles.Size());
//	AuxFunc::Fill(scaled_positions, VectorD::Zero());
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) {
//		scaled_positions[i] = e_particles.SA(i) * e_particles.X(i);
//	}
//	VectorD avg_scaled_position = AuxFunc::Sum(scaled_positions);
//	avg_scaled_position /= fluid_area;
//	return avg_scaled_position;
//}
//
//
//// update max velocity
//template<int d>
//void SPHBubbleFilm<d>::Decay_Velocity(Array<VectorD>& vels, const real dt, const real tang_strength, const real norm_strength)
//{
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) {
//		VectorD norm_V = AuxFunc::Component_Along(vels[i], e_particles.Normal(i));
//		VectorD tang_V = vels[i] - norm_V;
//		norm_V *= exp(-norm_strength * dt);
//		tang_V *= exp(-tang_strength * dt);
//		vels[i] = norm_V + tang_V;
//	}
//}
//
//template<int d>
//void SPHBubbleFilm<d>::Update_E(void)
//{
//	e_particles.Update();
//	if (analytical_boundary.Available()) e_particles.Update_Phis(analytical_boundary);
//}
//
//template<int d>
//void SPHBubbleFilm<d>::Update_E_Positions(const real dt)
//{
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) {
//		e_particles.X(i) += e_particles.V(i) * dt;
//	}
//	if (camera_moving_along && !ext_acc_func) {
//#pragma omp parallel for
//		for (int i = 0; i < e_particles.Size(); i++) {
//			e_particles.X(i) += init_COM - curr_COM;
//		}
//	}
//	// enforce boundary condition
//	if (analytical_boundary.Available()) e_particles.Correct_Position_With_Analytical_Boundary(dt, simulation_scale, analytical_boundary);
//}
//
//template<int d>
//void SPHBubbleFilm<d>::Update_C_Curvature(void) {
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) {
//		std::function<real(const int)> z_func_i = [&](const int idx)->real {return e_particles.X(idx).dot(e_particles.Normal(i)); };
//		real KS = e_particles.Surface_Laplacian(i, z_func_i, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//		std::function<real(const int)> z_func_i_b = [&](const int idx)->real {return boundary.X(idx).dot(e_particles.Normal(i)); };
//		real boundary_KS = boundary.Surface_Laplacian(e_particles.X(i), e_particles.E(i), e_particles.X(i).dot(e_particles.Normal(i)), z_func_i_b, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//		real total_KS = KS + boundary_KS;
//		if (Eulerian_Near_Boundary(i, neighbor_params.sph_radius)) {
//			real max_abs_KS = 2. / simulation_scale;
//			if (abs(total_KS) > max_abs_KS) { total_KS *= max_abs_KS / abs(total_KS); }
//		}
//		e_particles.KS(i) = total_KS;
//	}
//	////smooth out KS
//	for (int i = 0; i < numeric_params.XSPH_KS_Passes; i++) {
//		XSPH_Smoothing<real>(0.999, e_particles.KSRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//	}
//
//}
//
//template<int d>
//void SPHBubbleFilm<d>::Update_E_Geometry(void) {
//	double begin_time, end_time;
//	if (verbose) {
//		std::cout << "[E Geometry] begin" << std::endl;
//		begin_time = omp_get_wtime();
//	}
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) {
//		// Compute the number density of each E particle
//		// Must remember to account for the boundary as well, so that the SA for particles near the boundary is not overly large
//		e_particles.NDen(i) = e_particles.Surface_Sum_Value(e_particles.X(i), e_particles.E(i), ones_e, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC)
//			+ boundary.Surface_Sum_Value(e_particles.X(i), e_particles.E(i), ones_b, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//		// SA = 1/NDen
//		e_particles.SA(i) = 1. / e_particles.NDen(i);
//		// H = LVol/SA
//		e_particles.H(i) = e_particles.LVol(i) / e_particles.SA(i);
//		// Clip to avoid singularity for things like 1/h
//		real e_min_h = 0.1 * fluid_volume / fluid_area;
//		e_particles.H(i) = std::max<real>(e_particles.H(i), e_min_h);
//		// Den = LM/SA
//		e_particles.Den(i) = e_particles.LM(i) / e_particles.SA(i);
//		// Gamma = Soap/SA
//		e_particles.Gamma(i) = e_particles.Soap(i) / e_particles.SA(i);
//		// ST = st_water - RT * gamma
//		e_particles.ST(i) = 7.275 * 0.01 - 8.314 * 298.15 * e_particles.Gamma(i);
//		// The mass of a e_particle is the summation of its own mass and the mass of the fluid contained in it
//		e_particles.M(i) = e_particles.LM(i) + e_particles.base_mass;
//		// Same with the the volume
//		e_particles.Vol(i) = e_particles.LVol(i) + e_particles.base_mass / e_particles.base_rho;
//	}
//
//	// if bursting, then near the rim, the estimation will be inaccurate
//	// hence we set the Gamma to be the average value.
//	if (bursting) {
//		avg_gamma = AuxFunc::Mean<real>(e_particles.GammaRef());
//#pragma omp parallel for
//		for (int i = 0; i < e_particles.Size(); i++) {
//			if (near_rim_confidence[i] > 0.)
//				e_particles.Gamma(i) = avg_gamma;
//		}
//		std::function<bool(const int)> near_rim = [&](const int idx)->bool {return (near_rim_confidence[idx]>0.);};
//		for (int pass = 0; pass < 3; pass++) {
//			XSPH_Smoothing(0.999, e_particles.GammaRef(), neighbor_params.kernel, 1. * neighbor_params.sph_radius, KernelType::QUINTIC, near_rim);
//		}
//	}
//
//	Update_C_Curvature();
//
//	fluid_area = AuxFunc::Sum(e_particles.SARef());
//	avg_l_sa = fluid_area / l_particles.Size();
//	avg_e_sa = fluid_area / e_particles.Size();
//
//	curr_COM = Compute_COM();
//	//std::cout << "Current COM: \n" << curr_COM << std::endl;
//
//	if (verbose) {
//		end_time = omp_get_wtime();
//		std::cout << "[E Geometry] done. Time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
//	}
//}
//
//template<int d>
//void SPHBubbleFilm<d>::Delete_Particles(void)
//{
//	real avg_NDen = AuxFunc::Mean<real>(e_particles.NDenRef());
//	real delete_threshold = 150;
//	// Determine which E particles should be deleted
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) {
//		if (e_particles.Size() < 0.02 * num_E_preburst) {
//			to_delete_e[i] = 1;
//		}
//		real scores_loc = abs(e_particles.KS(i));
//		if (scores_loc > delete_threshold
//			|| e_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.sph_radius).size() < 7
//			|| e_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), 3. * neighbor_params.sph_radius).size() < 200)// || e_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.sph_radius).size() < 10)
//			to_delete_e[i] = 1;
//	}
//
//	// Add deleted E particles to 3d solver
//	for (int i = 0; i < e_particles.Size(); i++) {
//		if (to_delete_e[i]) {
//			if (RandomFunc::Random_Real(0., 1.) < (1./(real)dilution_rate_3d))
//			fluid_3d->Add_Particle(e_particles.X(i), .6 * e_particles.V(i) + -4. * rim_normals[i].normalized() * simulation_scale, 2.5684e-08);
//		}
//	}
//
//	if (!fluid_3d->Is_Trivial()) {
//		fluid_3d->Update_Max_Velocity(); //need to update max velocity once the new particles are transferred over
//	}
//	int original_e = e_particles.Size();
//	int original_l = l_particles.Size();
//	int remaining_e = e_particles.Delete_Elements_Safe(to_delete_e);
//	if (remaining_e < original_e) {
//		Update_E();
//		if (!bursting) {
//			bursting = true;
//		}
//	}
//
//	// Determine which L particles should be deleted
//#pragma omp parallel for
//	for (int i = 0; i < l_particles.Size(); i++) {
//		Array<int> nbs = e_particles.nbs_searcher->Find_Neighbors(l_particles.X(i), 0.5 * neighbor_params.sph_radius);
//		if (nbs.size() < 1) to_delete_l[i] = 1;
//	}
//	int remaining_l = l_particles.Delete_Elements_Safe(to_delete_l);
//	if (remaining_l < original_l) {
//		Update_L();
//	}
//
//	if (verbose) {
//		std::cout << "[Delete Particles] Deleted: " << original_e - remaining_e << " E particles" << std::endl;
//		std::cout << "[Delete Particles] Deleted: " << original_l - remaining_l << " L particles" << std::endl;
//	}
//
//	if (remaining_l < original_l || remaining_e < original_e) {
//		Reinit_Arrays();
//	}
//}
//
//template<int d>
//void SPHBubbleFilm<d>::Update_Dynamics(const real dt) {
//	// this is to compute raw velocity (before solving)
//	Apply_Ext_Forces(dt);
//	// raw velocity computed
//	 
////	//Update tang velocity (LV_T)
//	Update_Tang_Dynamics(dt);
////	//Update norm velocity
//	Update_Norm_Dynamics(dt);
//	//smooth out e_v
//	XSPH_Smoothing_Tang_Norm(numeric_params.XSPH_V_Tang, numeric_params.XSPH_V_Norm, e_particles.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//	//additional smoothing of E.V
//	for (int i = 0; i < numeric_params.XSPH_V_Passes; i++) {
//		XSPH_Smoothing_Tang_Norm(0., 0.999, e_particles.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//	}
//}
//
//template<int d>
//void SPHBubbleFilm<d>::Advance_E(const real dt, const int iter) {
//	// Update position of each E particle
// 	Update_E_Positions(dt);
//	// Update E NB searcher
//	Update_E();
//	// Similar to pressure projection, but done after the new surface is formed
//	if (!bursting) Even_Out_E(dt, false, max_even_out_e_steps); //if bursting or fixed_surface then we don't even out E
//}
//
//
//template<int d>
//void SPHBubbleFilm<d>::Update_Capillary_Forces(void)
//{
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) {
//		real pressure_difference = e_particles.ST(i) * e_particles.KS(i);//capillary pressure difference
//		VectorD force = pressure_difference * avg_e_sa * e_particles.Normal(i);
//		e_particles.F(i) += force;
//	}
//}
//
//// Make L (Material Mass/Volume) Evenly Distributed
//// if "repacking" is set to false, then we don't apply the e_particles.V(i), since we assume that it's in the middle of simulation
//// and the velocity is applied in advance_E
//// if repacking is true, then we apply the tangential component of e_particles.V(i), so that the velocity is not forgotten each step
//template<int d>
//void SPHBubbleFilm<d>::Even_Out_E(const real _dt, const bool repacking, const int _max_steps)
//{
//	double begin_time, end_time;
//	if (verbose) {
//		std::cout << "[E redistribute] begin" << std::endl;
//		begin_time = omp_get_wtime();
//	}
//	real avg_err = -1.; // initial error is 0
//
//	int max_steps = _max_steps;
//	if (bursting) max_steps = 1;
//
//	int num_iter = 0;
//	for (int k = 0; k < max_steps; k++) {
//#pragma omp parallel for
//		for (int i = 0; i < e_particles.Size(); i++) {
//			e_particles.FauxV(i) = VectorD::Zero();
//		}
//		real dt = _dt;
//		if (repacking) { // if repacking, then the "dt" is solely dependent upon the max e velocity
//			dt = (cfl * neighbor_params.e_dx) / std::max<real>(cfl * neighbor_params.e_dx * frame_rate, max_vel_e);
//			// this below basically does what Advance_E does before Even_out_E
//			// but with velocity only acting along the normal direction
//#pragma omp parallel for
//			for (int i = 0; i < e_particles.Size(); i++) {
//				//we will perform an extra projection of velocity, since we are not technically inside simulation yet...!
//				e_particles.V(i) = AuxFunc::Eliminate_Unit_Component(e_particles.V(i), e_particles.Normal(i));
//				e_particles.X(i) += e_particles.V(i) * dt;
//			}
//			// enforce boundary condition
//			if (analytical_boundary.Available()) e_particles.Correct_Position_With_Analytical_Boundary(dt, simulation_scale, analytical_boundary);
//			Update_E();
//		}
//
//		real elas_coeff = .8 / e_evenness; //may want to modify 0.8 here
//		std::cout << "elas_coeff without clip: " << elas_coeff << std::endl;
//		real max_elas_coeff = 1./AuxFunc::Mean(e_particles.NDenRef());
//		elas_coeff = std::min<real>(elas_coeff, max_elas_coeff); // may want to modify 100 here
//		std::cout << "max allowed elas_coeff: " << max_elas_coeff << std::endl;
//		std::cout << "elas_coeff: " << elas_coeff << std::endl;
//		//elas_coeff *= 1. / 4.;
//		real solver_err = Update_Faux_Tang_Dynamics(dt, elas_coeff);
//		//real solver_err = 0;
//
//		Update_E_Max_Velocity();
//		if (repacking) Decay_Velocity(e_particles.VRef(), dt, 10., 10.);
//
//		// <-----Recompute NDen----->
//#pragma omp parallel for
//		for (int i = 0; i < e_particles.Size(); i++) {
//			// Compute the number density of each E particle
//			// Must remember to account for the boundary as well, so that the SA for particles near the boundary is not overly large
//			e_particles.NDen(i) = e_particles.Surface_Sum_Value(e_particles.X(i), e_particles.E(i), ones_e, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC)
//				+ boundary.Surface_Sum_Value(e_particles.X(i), e_particles.E(i), ones_b, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//			// SA = 1/NDen
//			e_particles.SA(i) = 1. / e_particles.NDen(i);
//		}
//		real new_avg_err = 1./e_particles.Size() * (AuxFunc::Max<real>(e_particles.NDenRef()) - AuxFunc::Min<real>(e_particles.NDenRef()));
//
//		num_iter++;
//		if (verbose || 1==1) 
//			std::cout << "[E Redistribute] At iter: " << num_iter << ", nden error: " << std::setprecision(2) << new_avg_err << std::endl;
//		if (avg_err >= 0.) {
//			real err_diff = abs((new_avg_err - avg_err)) / avg_err;
//			if (err_diff < 1.e-3) {
//				if (verbose)
//					std::cout << "[E Redistribute] Converged." << std::endl;
//				//break;
//			}
//		}
//		avg_err = new_avg_err;
//		std::cout << "[E redistribute] finished iter: " << k << " error is: " << avg_err << std::endl;
//		if (!repacking && avg_err <= 3.) break;
//	}
//	if (verbose) {
//		end_time = omp_get_wtime();
//		std::cout << "[E Redistribute] done. Num iters: " << num_iter << ", time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
//	}
//}
//
//template<int d>
//void SPHBubbleFilm<d>::Initial_Packing(const int e_steps, const int l_steps)
//{
//	//Even out E and L
//	Even_Out_E(1., true, e_steps);
//	//Clear velocities used in the process
//	AuxFunc::Fill(e_particles.VRef(), VectorD::Zero());
//}
//
//template<int d>
//real SPHBubbleFilm<d>::Compute_Enclosed_Volume(void) {
//	VectorD origin = VectorD::Zero();
//	real vol = 0.;
//	for (int i = 0; i < e_particles.Size(); i++) {
//		real height = fabs(e_particles.Normal(i).dot(e_particles.X(i) - origin)); // volume of the skewed cone with the origin
//		real cone_vol = real(1. / d) * e_particles.SA(i) * height;
//		if ((e_particles.X(i)-origin).dot(e_particles.Normal(i))>0.) {
//			vol += cone_vol;
//		}
//		else {
//			vol -= cone_vol;
//		}
//	}
//	//std::cout << "Current enclosed volume: \n" << vol << std::endl;
//	return vol;
//}
//
//template<int d>
//void SPHBubbleFilm<d>::Update_Enclosed_Forces(void)
//{
//	//compute enclosed volume
//	real vol = Compute_Enclosed_Volume();
//	if (verbose) std::cout << "[Enclosed] curr enclosed volume: " << vol << std::endl;
//	real pressure = (vol > 1e-8) ? (enclosed_amount)/vol : outer_pressure;
//	if (verbose) std::cout << "[Enclosed] curr enclosed pressure: " << pressure << std::endl;
//	VectorD center = AuxFunc::Mean<VectorD>(e_particles.XRef());
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) {
//		real sgn = 1;
//		if ((e_particles.X(i)-center).dot(e_particles.Normal(i)) < 0) sgn = -1;
//		real pressure_difference = sgn * (pressure - outer_pressure);
//		VectorD force = pressure_difference * avg_e_sa * e_particles.Normal(i);
//		e_particles.F(i) += force;
//	}
//}
//
//// MAIN TANGENTIAL ACCELERATION COMPUTATION, difference form
//template<int d>
//void SPHBubbleFilm<d>::Update_Tang_Dynamics(const real dt, const bool internal_force_only) {// Solve for Gamma is a fashion similar to Huang 2020
//	double begin_time, end_time;
//	if (verbose) {
//		std::cout << "[Tang Dynamics] begin" << std::endl;
//		begin_time = omp_get_wtime();
//	}
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
//		// Initialize P to Gamma
//		e_particles.P(i) = e_particles.Gamma(i); // Here we are just using P as a placeholder for the revised/updated Gamma
//	}
//	if (dt < 1.e-8) return; // if dt is 0. then return.
//
//	real Ma = 2. * 0.5 * chem_params.R * chem_params.T_0 / chem_params.rho_0; // Marangoni param
//	real Dr = chem_params.drag_C / chem_params.rho_0; // Drag param
//
//	//Begin Solve
//	//real omega = 0.2; //Jacobi Relaxation
//	real omega = 0.1; //try a more conservative omega??
//	//int max_iter = 2. * simulation_scale / neighbor_params.e_dx; //max_iter
//	int max_iter = 100.;
//	int max_max_iter = 1000.;
//	Array<VectorT> temp; temp.resize(e_particles.Size()); //temp = [sum over j] {SA_j * 1/2 * grad_ij} (for computing diagonal terms)
//	Array<real> Bs; Bs.resize(e_particles.Size()); //the bouyancy variable B
//	Array<VectorT> grad_inv_h; grad_inv_h.resize(e_particles.Size()); //grad of 1/h
//	Array<VectorD> vel_diff; vel_diff.resize(e_particles.Size()); //vel_diff
//	Array<real> a_ii; a_ii.resize(e_particles.Size()); //diagonal elements of matrix A
//	Array<real> s_i; s_i.resize(e_particles.Size()); //source term
//	Array<VectorD> grad_p; grad_p.resize(e_particles.Size()); //grad_p = Surface_Gradient_Symmetric(P)
//	Array<real> LHS; LHS.resize(e_particles.Size()); //Left hand side
//	Array<real> err; err.resize(e_particles.Size()); //error
//
//	real avg_err = 0.;
//
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
//		//Compute temp
//		//temp[i] = e_particles.Surface_Sum_Grad(e_particles.X(i), e_particles.E(i), e_particles.SARef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//		//temp[i] += boundary.Surface_Sum_Grad(e_particles.X(i), e_particles.E(i), boundary.SARef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//		//temp[i] *= -1;
//		temp[i] = VectorT::Zero();
//		//Compute B 
//		Bs[i] = 1 - chem_params.expansion * (e_particles.Temp(i) - chem_params.T_0);
//		//Compute Velocity difference
//		VectorD vel_air = (wind_func != nullptr) ? wind_func(i) : e_particles.LV(i);
//		VectorD diff = vel_air - e_particles.LV(i);
//		vel_diff[i] = AuxFunc::Eliminate_Unit_Component(diff, e_particles.Normal(i));
//	}
//
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
//		Array<int> e_nbs = e_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.sph_radius);
//		Array<int> b_nbs = boundary.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.sph_radius);
//		a_ii[i] = 0.;
//		//Compute a_ii
//		real a_ii_0 = 0.; // zero order term
//		real a_ii_1 = 0.; // first order term
//		real a_ii_2 = 0.; // second order term
//		// zero order
//		a_ii_0 = -1. / (dt * e_particles.Gamma(i));
//		// first order // boundary is already taken care of in "temp"
//		std::function<real(const int)> inv_h = [&](const int idx)->real {return 1. / e_particles.H(idx); };
//		grad_inv_h[i] = e_particles.Surface_Gradient_Difference(i, inv_h, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC); //No need to compensate for the difference form
//		VectorT grad_B = e_particles.Surface_Gradient_Difference(i, Bs, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//		a_ii_1 = temp[i].dot(dt * Ma * grad_inv_h[i]);
//
//		for (int k = 0; k < e_nbs.size(); k++) {
//			int j = e_nbs[k];
//			VectorD wr_ji = e_particles.X(i) - e_particles.X(j);
//			VectorT sr_ji = PointSet<d>::Rotate_To_TPlane(wr_ji, e_particles.E(i));
//			VectorT grad_W = neighbor_params.kernel.template Grad<d - 1>(sr_ji, neighbor_params.sph_radius, KernelType::QUINTIC);
//			a_ii_2 += e_particles.SA(j) * (temp[i] + e_particles.SA(i) * grad_W).dot(grad_W);
//		}
//		for (int k = 0; k < b_nbs.size(); k++) {
//			int j = b_nbs[k];
//			VectorD wr_ji = e_particles.X(i) - boundary.X(j);
//			VectorT sr_ji = PointSet<d>::Rotate_To_TPlane(wr_ji, e_particles.E(i));
//			VectorT grad_W = neighbor_params.kernel.template Grad<d - 1>(sr_ji, neighbor_params.sph_radius, KernelType::QUINTIC);
//			a_ii_2 += boundary.SA(j) * (temp[i] + e_particles.SA(i) * grad_W).dot(grad_W);
//		}
//
//		//Compute Diagonal Elements a_ii
//		a_ii_2 *= -1.;
//		a_ii_2 *= dt * Ma * 1. / e_particles.H(i);
//		a_ii[i] = a_ii_0 + a_ii_1 + a_ii_2;
//		//Compute source term s_i
//		real velocity_divergence = e_particles.Surface_Divergence(e_particles.X(i), e_particles.E(i), e_particles.LV(i), e_particles.LVRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC) +
//			boundary.Surface_Divergence(e_particles.X(i), e_particles.E(i), e_particles.LV(i), boundary.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//		real gravity_divergence = Bs[i] * boundary.Surface_Divergence(e_particles.X(i), e_particles.E(i), g, boundary.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC) +
//			grad_B.dot(PointSet<d>::Project_To_TPlane(g, e_particles.E(i)));
//		real drag_divergence = e_particles.Surface_Divergence(e_particles.X(i), e_particles.E(i), vel_diff[i], vel_diff, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC) +
//			boundary.Surface_Divergence(e_particles.X(i), e_particles.E(i), vel_diff[i], boundary.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//		drag_divergence = grad_inv_h[i].dot(PointSet<d>::Project_To_TPlane(vel_diff[i], e_particles.E(i))) + 1. / e_particles.H(i) * drag_divergence;
//		s_i[i] = velocity_divergence - 1. / dt;
//		if (!internal_force_only) {
//			s_i[i] += dt * gravity_divergence + dt * Dr * drag_divergence;
//		}
//	}
//
//	// Iterate
//	int num_iters = 0;
//	if (bursting) max_max_iter = 0;
//	for (int l = 0; l < max_max_iter; l++) {
//#pragma omp parallel for
//		for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
//			VectorT grad2d = VectorT::Zero();
//			//Compute grad p
//			grad2d = e_particles.Surface_Gradient_Difference(i, e_particles.PRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//			VectorD grad3d;
//			e_particles.surface.Unproject_To_World(grad2d, e_particles.E(i), grad3d);
//			grad_p[i] = grad3d;
//		}
//#pragma omp parallel for
//		for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
//			//Compute LHS
//			// zero order
//			LHS[i] = -1. / (dt * e_particles.Gamma(i)) * e_particles.P(i);
//			// first order
//			LHS[i] += (dt * Ma * grad_inv_h[i]).dot(PointSet<d>::Project_To_TPlane(grad_p[i], e_particles.E(i)));
//			// second order
//			real laplace = e_particles.Surface_Divergence(e_particles.X(i), e_particles.E(i), grad_p[i], grad_p, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC) +
//				boundary.Surface_Divergence(e_particles.X(i), e_particles.E(i), grad_p[i], boundary.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//			LHS[i] += dt * Ma * 1. / e_particles.H(i) * laplace;
//			//Update pressure
//			//e_particles.P(i) = std::max<real>(0, e_particles.P(i) + omega / a_ii[i] * (s_i[i] - LHS[i]));
//			e_particles.P(i) = e_particles.P(i) + omega / a_ii[i] * (s_i[i] - LHS[i]);
//			err[i] = abs(LHS[i] - s_i[i]);
//		}
//
//		num_iters++;
//		real new_avg_err = AuxFunc::Mean<real>(err);
//		real err_change = avg_err > 0. ? (avg_err - new_avg_err) / avg_err : 0.;
//		if (verbose)
//			std::cout << "[Tang Dynamics Solver] at iter: " << num_iters << " error: " << new_avg_err << " improvement: " << err_change * 100 << "%" << std::endl;
//		if (avg_err > 0. && (abs(err_change) < 1.e-4 || err_change < 0.)) {
//			if (verbose)
//				std::cout << "[Tang Dynamics] Converged." << std::endl;
//			break;
//		}
//		if (num_iters > max_iter && new_avg_err < 0.1) {
//			break;
//		}
//		avg_err = new_avg_err;
//	}
//	if (verbose) {
//		std::cout << "[Tang Dynamics] done. Num iters: " << num_iters << " final error: " << avg_err << std::endl;
//	}
//
//	//Advance the velocity
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
//		//Compute pressure acceleration ap (ap3d)
//		VectorT ap_2d = VectorT::Zero();
//		VectorD ap_3d;
//
//		ap_2d = -Ma / e_particles.H(i) * (e_particles.Surface_Gradient_Difference(i, e_particles.PRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC));
//
//		VectorD increment = VectorD::Zero();
//
//		e_particles.surface.Unproject_To_World(ap_2d, e_particles.E(i), ap_3d);
//		increment += dt * ap_3d;
//		if (!internal_force_only) {
//			increment += dt * Dr / e_particles.H(i) * AuxFunc::Eliminate_Unit_Component(vel_diff[i], e_particles.Normal(i)) + dt * Bs[i] * AuxFunc::Eliminate_Unit_Component(g, e_particles.Normal(i));
//		}
//		if (is_rim_confidence[i] > 0.) { increment *= 0; }
//		e_particles.LV_T(i) += increment;
//	}
//
//	if (verbose) {
//		end_time = omp_get_wtime();
//		std::cout << "[Tang Dynamics] done. Time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
//	}
//}
//
//template class SPHBubbleFilm<2>;
//template class SPHBubbleFilm<3>;
