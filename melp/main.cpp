#include <iostream>
#include <variant>
#include "ParseArgs.h"
#include "ParticleALEFilmDriver.h"
#include "ParticleALEFoamDriver.h"

#ifndef __Main_cpp__
#define __Main_cpp__

void Set_Threads(ParseArgs& parse_args) {
	int number_threads = parse_args.Get_Integer_Value("-tnum");
	omp_set_num_threads(number_threads);
	int max_threads = omp_get_max_threads();
	std::cout << "Set " << number_threads << " threads, run with " << max_threads << " cores\n";
}


int main(int argc,char* argv[])
{

	//RandomNumber unit_generator(0, 1);
	//std::cout << "[A] From MAIN: Getting a random number between 0, 1: " << unit_generator.Value() << std::endl;
	//std::cout << "[B] From MAIN: Getting another random number between 0, 1: " << RandomFunc::Random_Real(0., 1.) << std::endl;

	const int d = 3;

    ParseArgs parse_args;
    parse_args.Add_String_Argument("-o","output","output path");
    parse_args.Add_Integer_Argument("-s",400,"resolution");
	parse_args.Add_Integer_Argument("-test", 1, "test");
	parse_args.Add_Integer_Argument("-driver",1,"driver");
	parse_args.Add_Integer_Argument("-lf",100,"last frame");
	parse_args.Add_Double_Argument("-cfl",1,"cfl number");
	parse_args.Add_Integer_Argument("-tnum", 1, "number of threads");
	parse_args.Add_Integer_Argument("-verbose", 0, "verbose or not");
	parse_args.Add_Integer_Argument("-diagnosis", 0, "whether to print diagnosis");
	parse_args.Add_Integer_Argument("-iomode", 0, "if write essentials or write everything"); //iomode = 0 means verbose, iomode = 1 means essentials
	parse_args.Add_Integer_Argument("-timing", 0, "whether output timing information");
	parse_args.Add_Integer_Argument("-from", 0, "first farme");
	parse_args.Add_Integer_Argument("-stride", 50, "snapshot stride");
	parse_args.Add_String_Argument("-txt", "params.txt", "Description File");
	parse_args.Add_String_Argument("-e_init", "", "path to e intiialization");
	parse_args.Add_String_Argument("-l_init", "", "path to l intiialization");
	parse_args.Add_Integer_Argument("-e_repack_steps", 0, "how many e_repack steps");
	parse_args.Add_Integer_Argument("-l_repack_steps", 0, "how many l_repack steps");
    parse_args.Parse(argc,argv);

	Set_Threads(parse_args);

    std::string output_dir=parse_args.Get_String_Value("-o");
    const int scale=parse_args.Get_Integer_Value("-s");
	const int driver=parse_args.Get_Integer_Value("-driver");
	const int test=parse_args.Get_Integer_Value("-test");
	const int last_frame=parse_args.Get_Integer_Value("-lf");
	const int verbose = parse_args.Get_Integer_Value("-verbose");
	const int timing = parse_args.Get_Integer_Value("-timing");
	const int first_frame = parse_args.Get_Integer_Value("-from");
	const int iomode = parse_args.Get_Integer_Value("-iomode");
	std::string e_init = parse_args.Get_String_Value("-e_init");
	std::string l_init = parse_args.Get_String_Value("-l_init");
	const int e_repack_steps = parse_args.Get_Integer_Value("-e_repack_steps");
	const int l_repack_steps = parse_args.Get_Integer_Value("-l_repack_steps");
	const real cfl=parse_args.Get_Double_Value("-cfl");
	switch(driver){
	case 1: { // single bubble driver
		ParticleALEFilmDriver<d> driver;
		driver.scale = scale;
		driver.output_dir = output_dir;
		driver.test = test;
		driver.last_frame = last_frame;
		driver.snapshot_stride = parse_args.Get_Integer_Value("-stride");
		driver.cfl = cfl;
		driver.verbose = verbose;
		driver.first_frame = first_frame;
		driver.iomode = iomode;
		driver.param_file_name = parse_args.Get_String_Value("-txt");
		driver.e_init = e_init;
		driver.l_init = l_init;
		driver.e_repack_steps = e_repack_steps;
		driver.l_repack_steps = l_repack_steps;
		driver.Initialize();
		driver.Run();
		//driver.fluid.log_output.close();
	} break;
	case 2: { // multi bubble driver
		ParticleALEFoamDriver<d> driver;
		driver.scale = scale;
		driver.output_dir = output_dir;
		driver.test = test;
		driver.last_frame = last_frame;
		driver.snapshot_stride = parse_args.Get_Integer_Value("-stride");
		driver.cfl = cfl;
		driver.verbose = verbose;
		driver.first_frame = first_frame;
		//driver.iomode = iomode;
		//driver.param_file_name = parse_args.Get_String_Value("-txt");
		//driver.e_init = e_init;
		//driver.l_init = l_init;
		//driver.e_repack_steps = e_repack_steps;
		//driver.l_repack_steps = l_repack_steps;
		driver.Initialize();
		driver.Run();
		//driver.fluid.log_output.close();
	}break;
	}
}

#endif

//template<int d>
//real ParticleALEFilm<d>::Even_Out_E_Step(const real dt, const bool well_distributed) {
//	if (dt < 1.e-8) return -1.; // if dt is 0. then return. Otherwise would be numerically troublesome
//
//	// if E is initially well_distributed, we can afford a larger relaxation; if is not, we should use a smaller relaxation	
//	real omega = 0.05;// Relaxation
//	if (well_distributed) omega = 0.02;
//
//	//int max_iter = 2. * simulation_scale / neighbor_params.e_dx;
//	int max_iter = 100;
//	Array<Eigen::Matrix<real, d, d - 1>> vel_grads; vel_grads.resize(e_particles.Size());
//	Array<real> P; P.resize(e_particles.Size()); //temporary P
//	Array<VectorT> temp; temp.resize(e_particles.Size());
//	Array<real> a_ii; a_ii.resize(e_particles.Size());
//	Array<real> s_i; s_i.resize(e_particles.Size());
//	Array<VectorT> ap; ap.resize(e_particles.Size());
//	Array<VectorD> ap3d; ap3d.resize(e_particles.Size());
//	Array<real> lap; lap.resize(e_particles.Size());
//	Array<real> err; err.resize(e_particles.Size());
//	real avg_err = 0.;
//
//	//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
//		//Compute temp
//		temp[i] = e_particles.Surface_Sum_Grad(e_particles.X(i), e_particles.E(i), e_particles.SARef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//		temp[i] += boundary.Surface_Sum_Grad(e_particles.X(i), e_particles.E(i), boundary.SARef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//		//std::cout << "i is: " << i << std::endl;
//		//std::cout << "temp[i] is: \n" << temp[i] << std::endl;
//		//std::cout << "temp[i] norm is: \n" << temp[i].norm() << std::endl;
//		//std::cout << "position from center: " << e_particles.X(i).norm() / simulation_scale << std::endl;
//	}
//
//	//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
//		//Compute a_ii
//		a_ii[i] = 0.;
//		Array<int> nbs = e_particles.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.sph_radius);
//		for (int k = 0; k < nbs.size(); k++) {
//			int j = nbs[k];
//			VectorD wr_ji = e_particles.X(i) - e_particles.X(j);
//			VectorT sr_ji = PointSet<d>::Rotate_To_TPlane(wr_ji, e_particles.E(i));
//			VectorT grad_W = neighbor_params.kernel.Grad<d - 1>(sr_ji, neighbor_params.sph_radius, KernelType::QUINTIC);
//			a_ii[i] += 1. * (-temp[i] + e_particles.SA(i) * grad_W).dot(grad_W);
//			//a_ii[i] += 1. * (e_particles.SA(i) * grad_W).dot(grad_W);
//		}
//		real derrick = a_ii[i];
//		nbs = boundary.nbs_searcher->Find_Neighbors(e_particles.X(i), neighbor_params.sph_radius);
//		for (int k = 0; k < nbs.size(); k++) {
//			int j = nbs[k];
//			VectorD wr_ji = e_particles.X(i) - boundary.X(j);
//			VectorT sr_ji = PointSet<d>::Rotate_To_TPlane(wr_ji, e_particles.E(i));
//			VectorT grad_W = neighbor_params.kernel.Grad<d - 1>(sr_ji, neighbor_params.sph_radius, KernelType::QUINTIC);
//			a_ii[i] += 1. * (-temp[i] + e_particles.SA(i) * grad_W).dot(grad_W);
//			//a_ii[i] += 1. * (e_particles.SA(i) * grad_W).dot(grad_W);
//		}
//		real rose = a_ii[i] - derrick;
//		//Compute Diagonal Elements a_ii
//		a_ii[i] *= -1. * dt * dt;
//		//if (a_ii[i] > 0.) {
//		//	std::cout << "a_ii > 0! i is: " << i << std::endl;
//		//	std::cout << "position from center: " << e_particles.X(i).norm() / simulation_scale << std::endl;
//		//	std::cout << "position: \n" << e_particles.X(i) << std::endl;
//		//}
//		//if (i == 9) {
//		//if (a_ii[i] > 0.) {
//		//	std::cout << "a_ii > 0! i is: " << i << std::endl;
//		//	std::cout << "e contribution: " << derrick << " this should be positive!" << std::endl;
//		//	std::cout << "b contribution: " << rose << " this should also be positive!" << std::endl;
//		//}
//		//Compute source term s_i
//		real density_derivatives = e_particles.Surface_Number_Density_Derivative(e_particles.X(i), e_particles.E(i), e_particles.V(i), e_particles.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC) +
//			boundary.Surface_Number_Density_Derivative(e_particles.X(i), e_particles.E(i), e_particles.V(i), boundary.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//		s_i[i] = avg_NDen - (e_particles.NDen(i) + dt * density_derivatives);
//	}
//
//	// Iterate
//	for (int l = 0; l < max_iter; l++) {
//#pragma omp parallel for
//		for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
//			//Compute pressure acceleration ap (ap3d)
//			ap[i] = -1. * (e_particles.Surface_Gradient_Difference(i, P, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC));
//			VectorD ap_3d; // don't forget to unproject back to world coordinate
//			e_particles.surface.Unproject_To_World(ap[i], e_particles.E(i), ap_3d);
//			ap3d[i] = ap_3d;
//		}
//#pragma omp parallel for
//		for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
//			temp_LV[i] = 0.01 * ap3d[i];
//			//if (i == 9) std::cout << "pressure before: " << P[i] << std::endl;
//			//Compute Laplacian
//			real density_derivative_e = e_particles.Surface_Number_Density_Derivative(e_particles.X(i), e_particles.E(i), ap3d[i], ap3d, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//			real density_derivative_b = boundary.Surface_Number_Density_Derivative(e_particles.X(i), e_particles.E(i), ap3d[i], boundary.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//			real density_derivative = density_derivative_e + density_derivative_b;
//			lap[i] = dt * dt * density_derivative;
//			err[i] = abs(lap[i] - s_i[i]) / avg_NDen;
//			//if (i == 9) std::cout << "density derivative_e: " << density_derivative_e << std::endl;
//			//if (i == 9) std::cout << "density derivative_b: " << density_derivative_b << std::endl;
//			//if (i == 9) std::cout << "density derivative: " << density_derivative << std::endl;
//			//if (i == 9) std::cout << "dt2: " << dt * dt << std::endl;
//			//if (i == 9) std::cout << "LHS: " << lap[i] << std::endl;
//			//if (i == 9) std::cout << "RHS: " << s_i[i] << std::endl;
//			//if (i == 9) std::cout << "a_ii: " << a_ii[i] << std::endl;
//			//Update pressure
//			P[i] = std::max<real>(0, P[i] + omega / a_ii[i] * (s_i[i] - lap[i]));
//			//P[i] = P[i] + omega / a_ii[i] * (s_i[i] - lap[i]);
//			e_particles.P(i) = 1.E-3 * P[i];
//			//if (i == 9) std::cout << "pressure after: " << P[i] << std::endl;
//		}
//		real new_avg_err = AuxFunc::Mean<real>(err);
//		real err_change = avg_err > 0. ? (avg_err - new_avg_err) / avg_err : 0.;
//		if (verbose)
//			std::cout << "[E Redistribute Solver] At iter: " << l + 1 << " error: " << new_avg_err << " improvement: " << err_change * 100 << "%" << std::endl;
//		if (avg_err > 0. && abs(err_change) < 1.e-4) {
//			if (verbose)
//				std::cout << "[E Redistribute Solver] Converged." << std::endl;
//			break;
//		}
//		avg_err = new_avg_err;
//	}
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
//		vel_grads[i] = e_particles.Surface_Jacobian_Difference(i, e_particles.VRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC); //compute the velocity gradient
//	}
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
//		ap[i] = -1. * (e_particles.Surface_Gradient_Difference(i, P, neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC));
//		VectorD ap_3d; // don't forget to unproject back to world coordinate
//		e_particles.surface.Unproject_To_World(ap[i], e_particles.E(i), ap_3d);
//		e_particles.FauxV(i) += dt * ap_3d;
//	}
//	// same way with e.V, we mock the viscosity near the boundary
//	if (sticky_boundary) {
//		Boundary_Adhesion(e_particles.FauxVRef());
//	}
//
//	//smooth out e_faux_v
//	XSPH_Smoothing<VectorD>(numeric_params.XSPH_Even_E, e_particles.FauxVRef(), neighbor_params.kernel, neighbor_params.sph_radius, KernelType::QUINTIC);
//
//#pragma omp parallel for
//	for (int i = 0; i < e_particles.Size(); i++) { // for each E particle
//		VectorD displacement3 = e_particles.FauxV(i) * dt;
//		VectorT displacement2 = PointSet<d>::Project_To_TPlane(displacement3, e_particles.E(i));
//		e_particles.X(i) += displacement3;
//		VectorD velocity_change = vel_grads[i] * displacement2;
//		e_particles.V(i) += velocity_change;
//	}
//
//	if (analytical_boundary.Available()) e_particles.Correct_Position_With_Analytical_Boundary(dt, simulation_scale, analytical_boundary);
//	Update_E();
//
//	return avg_err;
//}
