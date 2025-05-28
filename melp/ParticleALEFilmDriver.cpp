#include "ParticleALEFilmDriver.h"
#include "MeshFunc.h"
#include "PointSetFunc.h"
#include "ArrayIO.h"
#include "FileInit.h"
#include "RenderFunc.h"
#include "EulerInitializer.h"

template<int d>
real ParticleALEFilmDriver<d>::Temperature(VectorD pos, real y_lowest, real y_highest, real strength, real perturbation_ratio)
{
	real y_range = y_highest - y_lowest;
	real coldest = 200;
	real hottest = 100 * strength + coldest;
	real temp_range = hottest - coldest;
	real x = pos[0];
	real y = pos[1];
	real z = pos[2];
	real base_temp = 300;
	if (y <= y_lowest) {
		base_temp = hottest;
	}
	else if (y >= y_highest) {
		base_temp = coldest;
	}

	real lerp = 1 - (y - y_lowest) / y_range; //1 at bottom, 0 at top
	real interp = lerp * lerp;
	base_temp = coldest + (temp_range) * interp;
	real x1 = 5 * x / fluid.simulation_scale;
	real y1 = 5 * y / fluid.simulation_scale;
	real z1 = 5 * z / fluid.simulation_scale;
	real temp_scale = 1. * strength * perturbation_ratio; 
	real temp1 = 0.9 * temp_scale;
	real temp2 = 1.5 * temp_scale;
	real temp3 = 1 * temp_scale;
	real temp4 = 2 * temp_scale;
	real temp5 = 3 * temp_scale;
	real temp6 = 0.6 * temp_scale;
	return base_temp + temp1 * sin(8 * x1 - 0.5) + temp2 * cos(15 * x1 + 1) + temp3 * sin(12 * y1 - 1.5)
			   + temp4 * cos(7 * y1 + 2) + temp5 * sin(10 * z1 - 2.5) + temp6 * cos(11 * z1 + 3);
}


template<int d>
void ParticleALEFilmDriver<d>::Case_1(void) {//sphere
	std::cout << "[Initialization] begin" << std::endl;
	double begin_time = omp_get_wtime();
	if constexpr (d == 3) {
		//---simulation params---//
		max_iter_per_frame = -1;
		cfl = 2.;
		frame_rate = 30.;
		real R = (real)1.;
		int fineness = 2;
		fluid.simulation_scale = R;
		fluid.simulation_fineness = fineness;
		fluid.verbose = verbose;

		real e_dx = 1.; real l_dx = 1.;
		e_dx = PointSetFunc::Initialize_Sphere_Points_Regular(VectorD::Zero(), R, 1200, fluid.e_particles);
		l_dx = PointSetFunc::Initialize_Sphere_Points_Regular(VectorD::Zero(), R, 20000, fluid.l_particles);

#pragma omp parallel for
		for (int i = 0; i < fluid.e_particles.Size(); i++) {
			real multiplier = Perlin_Noise(fluid.e_particles.X(i), 1234, 0.5 / R, .2, 2);
			fluid.e_particles.X(i) += R * multiplier * fluid.e_particles.Normal(i);
			fluid.e_particles.X(i)[0] *= 1.4;
		}
#pragma omp parallel for
		for (int i = 0; i < fluid.l_particles.Size(); i++) {
			real multiplier = Perlin_Noise(fluid.l_particles.X(i), 1234, 0.5 / R, .2, 2);
			fluid.l_particles.X(i) += R * multiplier * fluid.l_particles.Normal(i);
			fluid.l_particles.X(i)[0] *= 1.4;
		}

		fluid.e_only = fluid.l_particles.Size() < 1;

		if (verbose) {
			std::cout << "[Initialization] e_dx is: " << e_dx << std::endl;
			std::cout << "[Initialization] e_number is: " << fluid.e_particles.Size() << std::endl;
			std::cout << "[Initialization] l_dx is: " << l_dx << std::endl;
			std::cout << "[Initialization] l_number is: " << fluid.l_particles.Size() << std::endl;
		}

		//---configure NB params---//
		real e_np_on_h = 4;
		fluid.neighbor_params = NeighborParams<d>(e_dx, e_np_on_h, e_np_on_h);

		//---configure chem and numeric params---//
		fluid.chem_params = ChemicalParams();
		fluid.numeric_params = NumericalParams();

		//---initialize Lagrangian quantities---//
		real l_thickness = 5e-7; // typical thickness = 500 nm
		real e_thickness = 1e-8; // newton black film = 10 nm
		real l_gamma = 1e-7; // soap conc = 1e-7 mol/m^2
		real rho = fluid.chem_params.rho_0;
		real T = fluid.chem_params.T_0; //temperature
		real rho_base = 1.e6 * rho;
		real area = 4 * pi * R * R; // for sphere
		real l_volume = area * l_thickness;
		fluid.fluid_volume = l_volume;
		fluid.fluid_area = area;
		fluid.fluid_soap = area * l_gamma;

		fluid.e_particles.Set_Values(rho_base, rho_base * (area / fluid.e_particles.Size() * e_thickness));
		fluid.l_particles.Set_Values(area, l_thickness, l_gamma, rho, T);

		fluid.g = -.2 * VectorD::Unit(1) * 3.;
		fluid.camera_moving_along = true;
		fluid.Initialize(cfl, frame_rate, output_dir);

		if (first_frame < 1) {
			if (e_repack_steps >= 1) {
				fluid.Even_Out_E(1., true, e_repack_steps);
				Save_Snapshot(-1);
			}
			if (l_repack_steps >= 1) {
				fluid.Even_Out_L(l_repack_steps);
				Save_Snapshot(-1);
			}
		}

		fluid.enclosed = 1;
		fluid.enclosed_vol = fluid.Compute_Enclosed_Volume();
		if (verbose) std::cout << "[Initialization] init enclosed volume: " << fluid.enclosed_vol << std::endl;
		real young_laplace_pressure = AuxFunc::Max<real>(fluid.e_particles.STRef()) * 2 * 1./R * 100;
		fluid.enclosed_amount = fluid.enclosed_vol * (fluid.outer_pressure + young_laplace_pressure);

		// external forces
		std::function<real(const int)> temp_func = [&](const int idx)->real {
			return Temperature(fluid.l_particles.X(idx), -fluid.simulation_scale, fluid.simulation_scale, 2.);
		};
		fluid.temp_func = temp_func;

		// noise up the initial height
#pragma omp parallel for
		for (int i = 0; i < fluid.l_particles.Size(); i++) {
			real multiplier = 1. + Perlin_Noise(fluid.l_particles.X(i), 9199, 0.7/R, 0.5, 2);
			fluid.l_particles.Vol(i) *= multiplier;
			//fluid.l_particles.Soap(i) *= multiplier;
		}

		fluid.L2E();
		fluid.Update_E_Geometry();
		fluid.Update_L_Geometry(0.);

		if (verbose) {
			fluid.avg_h = AuxFunc::Mean<real>(fluid.e_particles.HRef());
			if (verbose) std::cout << "[Initialization] init avg h: " << fluid.avg_h << std::endl;
			Array<int> nbs;
			nbs = fluid.e_particles.nbs_searcher->Find_Neighbors(fluid.e_particles.X(0), fluid.neighbor_params.sph_radius);
			std::cout << "[Initialization] E finds E num nbs: " << nbs.size() << std::endl; 
			nbs = fluid.l_particles.nbs_searcher->Find_Neighbors(fluid.e_particles.X(0), fluid.neighbor_params.interp_radius);
			std::cout << "[Initialization] E finds L num nbs: " << nbs.size() << std::endl;
			if (!fluid.e_only) {
				nbs = fluid.e_particles.nbs_searcher->Find_Neighbors(fluid.l_particles.X(0), fluid.neighbor_params.interp_radius);
				std::cout << "[Initialization] L finds E num nbs: " << nbs.size() << std::endl;
			}
		}
	}
	double end_time = omp_get_wtime();
	std::cout << "[Initialization] done. Time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
}

template<int d>
void ParticleALEFilmDriver<d>::Write_Particle_Basics(const NAParticles<d>& particles, const std::string &frame_dir, const std::string& prefix)
{
	BinaryDataIO::Write_Vector_Array_3D<real, d>(frame_dir + "/" + prefix + "_x_bin", particles.XRef());
	//particles.Write_To_File_3d_Fast(frame_dir + "/" + prefix + "_tracker_points");
	RenderFunc::Write_Points_Float<d, real>(frame_dir + "/" + prefix + "_tracker_points", particles.XRef());
}

template<int d>
void ParticleALEFilmDriver<d>::Write_Output_Files(const int frame)
{
	namespace fs = std::filesystem;

	static double begin_time = omp_get_wtime();
	static double last_time = begin_time;
	static double now_time = last_time;

	Base::Write_Output_Files(frame);

	fs::path frame_path = fs::path(output_dir) / std::to_string(frame);

	//write snapshot first things first
	if (frame % snapshot_stride == (int)0) Save_Snapshot(frame);

	Array<VectorD> e_normals; e_normals.resize(fluid.e_particles.Size());
#pragma omp parallel for
	for (int i = 0; i < e_normals.size(); i++) e_normals[i] = fluid.e_particles.Normal(i);

	Array<VectorD> l_normals; l_normals.resize(fluid.l_particles.Size());
#pragma omp parallel for
	for (int i = 0; i < l_normals.size(); i++) l_normals[i] = fluid.l_particles.Normal(i);

	RenderFunc::Write_Vectors_Float<d, real>((frame_path / "e_tracker_circles").string(), fluid.e_particles.XRef(), e_normals, true);
	RenderFunc::Write_Vectors_Float<d, real>((frame_path / "l_tracker_circles").string(), fluid.l_particles.XRef(), l_normals, true);
	BinaryDataIO::Write_Array((frame_path / "l_h_bin").string(), fluid.l_particles.HRef());
	BinaryDataIO::Write_Array((frame_path / "l_sa_bin").string(), fluid.l_particles.SARef());
	if (!fluid.e_only) {
		BinaryDataIO::Write_Array((frame_path / "e_h_bin").string(), fluid.e_particles.HRef());
	}
	else {
		Array<real> e_h; e_h.resize(fluid.e_particles.Size());
		real avg_nden = AuxFunc::Mean<real>(fluid.e_particles.NDenRef());
#pragma omp parallel for
		for (int i = 0; i < e_h.size(); i++) e_h[i] = 3.e-7 * fluid.e_particles.NDen(i)/avg_nden;
		BinaryDataIO::Write_Array(frame_dir + "/e_h_bin", e_h);
	}
	//BinaryDataIO::Write_Array(frame_dir + "/e_h_bin", fluid.temp_scalar);
	BinaryDataIO::Write_Array(frame_dir + "/e_sa_bin", fluid.e_particles.SARef());

	RenderFunc::Write_Points_Float<d, real>(frame_dir + "/v_tracker_points", fluid_3d.particles.XRef());
	// That should be all that's needed for final rendering 
	if (iomode == 1) return;

	//tracker points
	RenderFunc::Write_Points_Float<d, real>(frame_dir + "/e_tracker_points", fluid.e_particles.XRef());
	RenderFunc::Write_Points_Float<d, real>(frame_dir + "/l_tracker_points", fluid.l_particles.XRef());

	//height and nden
	real h_scale = 0.3 * 1. / 5.e-7 * fluid.simulation_scale;
	real nden_scale = 1. / fluid.e_particles.Size() * fluid.simulation_scale;
	RenderFunc::Write_Customized_Segments_Float<d>(frame_dir + "/e_point_height", fluid.e_particles.XRef(), e_normals, fluid.e_particles.GammaRef(), h_scale);
	//RenderFunc::Write_Customized_Segments_Float<d>(frame_dir + "/e_point_height", fluid.e_particles.XRef(), e_normals, fluid.e_particles.GammaRef(), 3. * h_scale);
	RenderFunc::Write_Customized_Segments_Float<d>(frame_dir + "/e_point_nden", fluid.e_particles.XRef(), e_normals, fluid.e_particles.NDenRef(), nden_scale);
	//RenderFunc::Write_Customized_Segments_Float<d>(frame_dir + "/e_point_nden", fluid.e_particles.XRef(), e_normals, fluid.e_particles.KSRef(), 1.E-3 * h_scale);
	RenderFunc::Write_Customized_Segments_Float<d>(frame_dir + "/e_point_pressure", fluid.e_particles.XRef(), e_normals, fluid.temp_scalar, 1.e1);
	//RenderFunc::Write_Customized_Segments_Float<d>(frame_dir + "/e_point_pressure", fluid.e_particles.XRef(), e_normals, fluid.e_particles.PRef(), 1.E5);
	
	//velocities
	//RenderFunc::Write_Vectors_Float<d, real>(frame_dir + "/e_point_velocity", fluid.e_particles.XRef(), fluid.e_particles.BH_LVRef());
	RenderFunc::Write_Vectors_Float<d, real>(frame_dir + "/e_point_velocity", fluid.e_particles.XRef(), fluid.e_particles.VRef());
	//RenderFunc::Write_Vectors_Float<d, real>(frame_dir + "/e_point_force", fluid.e_particles.XRef(), fluid.temp_vector);
	//RenderFunc::Write_Vectors_Float<d, real>(frame_dir + "/e_point_force", fluid.e_particles.XRef(), fluid.e_particles.VRef());
	RenderFunc::Write_Vectors_Float<d, real>(frame_dir + "/e_point_force", fluid.e_particles.XRef(), fluid.e_particles.FRef());
	RenderFunc::Write_Vectors_Float<d, real>(frame_dir + "/l_point_velocity", fluid.e_particles.XRef(), fluid.e_particles.LV_TRef());

	//segment mesh
	PointSetFunc::Write_Local_Frames_To_File(frame_dir + "/e_segment_mesh", fluid.e_particles, 0.2 * fluid.simulation_scale);

	////boundary and 3d particles
	RenderFunc::Write_Points_Float<d, real>(frame_dir + "/boundary_points", fluid.boundary.XRef());
}

template<int d>
void ParticleALEFilmDriver<d>::Save_Snapshot(const int frame)
{
	namespace fs = std::filesystem;
	std::cout << "save snapshot for frame " << frame << "\n";// fluid.Numerical_Check_Force();
	std::string snapshot_dir = output_dir + "/snapshot/" + std::to_string(frame);
	if (!fs::exists(snapshot_dir)) {
		std::cout << "#     Create snapshot directory: " << snapshot_dir << std::endl;
		fs::create_directories(snapshot_dir);
	}
	//if (!File::Directory_Exists(snapshot_dir.c_str())) File::Create_Directory(snapshot_dir);
	fluid.e_particles.Save_Snapshot(snapshot_dir + "/e_particles.bin");
	fluid.l_particles.Save_Snapshot(snapshot_dir + "/l_particles.bin");
	fluid.boundary.Save_Snapshot(snapshot_dir + "/boundary_particles.bin");
	fluid_3d.particles.Save_Snapshot(snapshot_dir + "/3d_particles.bin");
}

template<int d>
void ParticleALEFilmDriver<d>::Load_Snapshot(const int frame)
{
	std::string snapshot_dir = output_dir + "/snapshot/" + std::to_string(frame);
	fluid.e_particles.Load_Snapshot(snapshot_dir + "/e_particles.bin");
	fluid.l_particles.Load_Snapshot(snapshot_dir + "/l_particles.bin");
	fluid.boundary.Load_Snapshot(snapshot_dir + "/boundary_particles.bin");
	fluid_3d.particles.Load_Snapshot(snapshot_dir + "/3d_particles.bin");
	if (fluid.e_particles.Size()) fluid.Update_E();
	if (fluid.l_particles.Size()) fluid.Update_L();
	if (fluid_3d.particles.Size()) { 
		fluid_3d.particles.Update(); 
		fluid.bursting = true;
	}
	fluid.Reinit_Arrays();
}

template<int d>
void ParticleALEFilmDriver<d>::Load_E_Initialization(const std::string e_path)
{
	std::string snapshot_dir = output_dir + "/../" + e_path + "/snapshot/" + std::to_string(-1);
	fluid.e_particles.Load_Snapshot(snapshot_dir + "/e_particles.bin");
}

template<int d>
void ParticleALEFilmDriver<d>::Load_L_Initialization(const std::string l_path)
{
	std::string snapshot_dir = output_dir + "/../" + l_path + "/snapshot/" + std::to_string(-1);
	fluid.l_particles.Load_Snapshot(snapshot_dir + "/l_particles.bin");
}


template class ParticleALEFilmDriver<2>;
template class ParticleALEFilmDriver<3>;
