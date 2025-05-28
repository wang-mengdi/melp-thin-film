#include "ParticleALEFoamDriver.h"
#include "MeshFunc.h"
#include "PointSetFunc.h"
#include "ArrayIO.h"
#include "FileInit.h"
#include "RenderFunc.h"
#include "EulerInitializer.h"


template<int d>
real ParticleALEFoamDriver<d>::Temperature(VectorD pos)
{
	real strength = .05;
	real perturbation_ratio = 1.;
	real y_highest = .5;
	real y_lowest = -.5;
	real y_range = y_highest-y_lowest;
	real coldest = 300;
	real hottest = 100 * 1. * strength + coldest;
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
	real interp = lerp;
	base_temp = coldest + (temp_range)*interp;
	/*return base_temp;*/
	real x1 = 4 * x / multi_fluid.simulation_scale;
	real y1 = 4 * y / multi_fluid.simulation_scale;
	real z1 = 4 * z / multi_fluid.simulation_scale;
	real temp_scale = 8 * strength * perturbation_ratio;
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
void ParticleALEFoamDriver<d>::Setup_Bubble(ParticleALEFilm<d>& fluid, const real radius, const VectorD center, const VectorD velocity, const real thickness_multiplier, const real fineness_multiplier, const std::string l_init, const std::string e_init) {//code to add a bubble
	if constexpr (d == 3) {
		real R = radius; //10cm radius
		int fineness = multi_fluid.simulation_fineness * fineness_multiplier;
		fluid.simulation_scale = R;
		fluid.simulation_fineness = fineness;
		fluid.verbose = verbose;
		int e_sub_num = 4 + (fineness - 1) + (R - 1);
		int l_sub_num = e_sub_num + 2;
		real base_dx = 0.0150084;
		VectorD ctr = center;
		real e_dx = 1.; real l_dx = 1.;
		e_dx = PointSetFunc::Initialize_Sphere_Points_Regular(center, R, 1200, fluid.e_particles);
		if (!multi_fluid.e_only) {
			l_dx = PointSetFunc::Initialize_Sphere_Points_Regular(center, R, 5000, fluid.l_particles);
		}
#pragma omp parallel for
		for (int i = 0; i < fluid.e_particles.Size(); i++) {
			real multiplier = ParticleALEFilmDriver<d>::Perlin_Noise(fluid.e_particles.X(i), 1234, 0.5 / R, .2, 2);
			fluid.e_particles.X(i) += R * multiplier * fluid.e_particles.Normal(i);
			//fluid.e_particles.X(i)[0] *= 1.4;
		}
#pragma omp parallel for
		for (int i = 0; i < fluid.l_particles.Size(); i++) {
			real multiplier = ParticleALEFilmDriver<d>::Perlin_Noise(fluid.l_particles.X(i), 1234, 0.5 / R, .2, 2);
			fluid.l_particles.X(i) += R * multiplier * fluid.l_particles.Normal(i);
			//fluid.l_particles.X(i)[0] *= 1.4;
		}

		real e_np_on_h = 4;
		real interp_np_on_h = 4;
		fluid.neighbor_params = NeighborParams<d>(e_dx, e_np_on_h, interp_np_on_h);

		//---configure chem and numeric params---//
		fluid.chem_params = ChemicalParams();
		fluid.chem_params.heat_capacity_param = 2.;
		fluid.numeric_params = NumericalParams();
		fluid.numeric_params.XSPH_V_Norm = 0.99;

		//---initialize Lagrangian quantities---//
		real l_thickness = 3e-7 * thickness_multiplier; // typical thickness = 500 nm
		real e_thickness = 1e-8; // newton black film = 10 nm
		real l_gamma = 1e-7 * thickness_multiplier;// *thickness_multiplier; // soap conc = 1e-7 mol/m^2
		real rho = fluid.chem_params.rho_0;
		real T = fluid.chem_params.T_0; //temperature
		real rho_base = 2.e5 * rho;
		real area = 4 * pi * R * R; // for sphere
		real l_volume = area * l_thickness;
		fluid.fluid_volume = l_volume;
		fluid.fluid_area = area;
		fluid.fluid_soap = area * l_gamma;
		fluid.l_geometry_sync = false;


		fluid.e_particles.Set_Values(rho_base, rho_base * (area / fluid.e_particles.Size() * e_thickness));
		if (!multi_fluid.e_only) {
			fluid.l_particles.Set_Values(area, l_thickness, l_gamma, rho, T);
		}

		fluid.g = multi_fluid.g;

#pragma omp parallel for
		for (int i = 0; i < fluid.e_particles.Size(); i++) {
			fluid.e_particles.V(i) = velocity;
		}

		if (!multi_fluid.e_only) {
#pragma omp parallel for
			for (int i = 0; i < fluid.l_particles.Size(); i++) {
				//std::cout << "this is called: " << i << std::endl;
				fluid.l_particles.V(i) = velocity;
			}

			// noise up the initial height
#pragma omp parallel for
			for (int i = 0; i < fluid.l_particles.Size(); i++) {
				real multiplier = 1. + ParticleALEFilmDriver<d>::Perlin_Noise(fluid.l_particles.X(i), 9199, 0.7 / R, 0.12, 2);
				fluid.l_particles.Vol(i) *= multiplier;
			}
		}

		fluid.Initialize(cfl, frame_rate, output_dir);

		fluid.enclosed = 1;
		fluid.enclosed_vol = fluid.Compute_Enclosed_Volume();
		if (verbose) std::cout << "[Foam Initialization] init enclosed volume: " << fluid.enclosed_vol << std::endl;
		real young_laplace_pressure = 100 * AuxFunc::Max<real>(fluid.e_particles.STRef()) * 2 * 1. / R;
		fluid.enclosed_amount = fluid.enclosed_vol * (multi_fluid.outer_pressure + young_laplace_pressure);
		if (verbose) std::cout << "[Foam Initialization] enclosed amount: " << fluid.enclosed_amount << std::endl;
		//fluid.enforce_momentum_conserve = true;

		if (verbose) {
			fluid.avg_h = AuxFunc::Mean<real>(fluid.e_particles.HRef());
			if (verbose) std::cout << "[Initialization] init avg h: " << fluid.avg_h << std::endl;
			Array<int> nbs;
			nbs = fluid.e_particles.nbs_searcher->Find_Neighbors(fluid.e_particles.X(0), fluid.neighbor_params.sph_radius);
			std::cout << "[Initialization] E finds E num nbs: " << nbs.size() << std::endl;
			if (!multi_fluid.e_only) {
				nbs = fluid.l_particles.nbs_searcher->Find_Neighbors(fluid.e_particles.X(0), fluid.neighbor_params.interp_radius);
				std::cout << "[Initialization] E finds L num nbs: " << nbs.size() << std::endl;
				nbs = fluid.e_particles.nbs_searcher->Find_Neighbors(fluid.l_particles.X(0), fluid.neighbor_params.interp_radius);
				std::cout << "[Initialization] L finds E num nbs: " << nbs.size() << std::endl;
			}
		}
	}
}


template<int d>
void ParticleALEFoamDriver<d>::Case_1(void) {//three bubbles
	std::cout << "[Foam Initialization] begin" << std::endl;
	double begin_time = omp_get_wtime();
	if constexpr (d == 3) {
		//---simulation params---//
		max_iter_per_frame = -1;
		cfl = .99;
		frame_rate = 30.;
		multi_fluid.simulation_fineness = 2;
		multi_fluid.simulation_scale = 1.;
		multi_fluid.g = VectorD::Zero();
		multi_fluid.outer_pressure = 0.;
		multi_fluid.e_dx = multi_fluid.simulation_scale; // set to be big, to be revised later
		multi_fluid.sph_radius = multi_fluid.simulation_scale; // set to be big, to be revised later
		int num_bubbles = 3;
		Array<real> Rs; 
		Rs.push_back(multi_fluid.simulation_scale); 
		Rs.push_back(multi_fluid.simulation_scale);
		Rs.push_back(multi_fluid.simulation_scale);
		Array<VectorD> ctrs; 
		ctrs.push_back(VectorD::Unit(1) * 1.1 * multi_fluid.simulation_scale); 
		ctrs.push_back(VectorD::Unit(0) * 2. * multi_fluid.simulation_scale);
		ctrs.push_back(-VectorD::Unit(1) * 1.1 * multi_fluid.simulation_scale); 
		Array<VectorD> vels; 
		//vels.push_back(VectorD::Zero()); vels.push_back(VectorD::Zero()); vels.push_back(VectorD::Zero());
		vels.push_back(0.3 * -VectorD::Unit(1) + 0.15 * VectorD::Unit(0));
		vels.push_back(-0.3 * VectorD::Unit(0));
		vels.push_back(0.3 * VectorD::Unit(1) + 0.15 * VectorD::Unit(0));
		Array<real> thicks; 
		thicks.push_back(2.5); 
		thicks.push_back(1.);
		thicks.push_back(1.5);
		for (int i = 0; i < num_bubbles; i++) {
			if (verbose) std::cout << "[Foam Initialization] Begin Initializing region: " << i << std::endl;
			std::shared_ptr<ParticleALEFilm<d>> fluid_ptr = std::make_shared<ParticleALEFilm<d>>();
			Setup_Bubble(*fluid_ptr, Rs[i], ctrs[i], vels[i], thicks[i]);
			multi_fluid.regions.push_back(fluid_ptr);
			// set the e_dx of multifluid to be the minimum of e_dx of all bubbles. Practically, all bubbles should have the same e_dx/fineness
			multi_fluid.e_dx = std::min<real>(multi_fluid.e_dx, fluid_ptr->neighbor_params.e_dx);
			multi_fluid.sph_radius = std::min<real>(multi_fluid.sph_radius, fluid_ptr->neighbor_params.sph_radius);
			if (verbose) std::cout << "[Foam Initialization] Done Initializing region: " << i << std::endl;
		}
	}

	multi_fluid.conserve_momentum_along.push_back(0); multi_fluid.conserve_momentum_along.push_back(1); multi_fluid.conserve_momentum_along.push_back(2);
	multi_fluid.Initialize(cfl, frame_rate, output_dir);

	for (int i = 0; i < multi_fluid.Num_Regions(); i++) {
		/*multi_fluid.region_modulus[i] = .5;*/
		multi_fluid.region_modulus[i] = 10.;
	}

	double end_time = omp_get_wtime();
	std::cout << "[Foam Initialization] done. Time taken: " << std::setprecision(2) << end_time - begin_time << std::endl;
}

template<int d>
void ParticleALEFoamDriver<d>::Write_Output_Files(const int frame)
{
	static double begin_time = omp_get_wtime();
	static double last_time = begin_time;
	static double now_time = last_time;

	Base::Write_Output_Files(frame);
	for (int i = 0; i < multi_fluid.Num_Regions(); i++) {
		//std::cout << "trying to write region: " << i << std::endl;
		Write_Output_Files_Region(frame, i);
	}

	////boundary and 3d particles
	RenderFunc::Write_Points_Float<d, real>(frame_dir + "/boundary_points", multi_fluid.boundary.XRef());
	RenderFunc::Write_Points_Float<d, real>(frame_dir + "/v_tracker_points", fluid_3d.particles.XRef());

	//Snapshot
	if (frame % snapshot_stride == (int)0) Save_Snapshot(frame);
}

template<int d>
void ParticleALEFoamDriver<d>::Write_Output_Files_Region(const int frame, const int idx)
{
	if (idx >= multi_fluid.Num_Regions()) return;

	if (multi_fluid.freeze[idx]) {
		std::string e_points_filename = frame_dir + "/e_tracker_points_" + std::to_string(idx);
		std::string e_tracker_circles_filename = frame_dir + "/e_tracker_circles_" + std::to_string(idx);
		std::ofstream output(e_points_filename, std::ios::binary); if (!output)return;
		File::Write_Binary(output, 0);
		output.close();
		std::ofstream output2(e_tracker_circles_filename, std::ios::binary); if (!output2)return;
		File::Write_Binary(output2, 0);
		output2.close();
		return;
	}


	//std::string region_dir = frame_dir + "/" + std::to_string(idx);
	//if (!File::Directory_Exists(region_dir.c_str()))File::Create_Directory(region_dir);
	ParticleALEFilm<d>& fluid = *multi_fluid.regions[idx];
 	Array<VectorD> e_normals; e_normals.resize(fluid.e_particles.Size());
#pragma omp parallel for
	for (int i = 0; i < e_normals.size(); i++) e_normals[i] = fluid.e_particles.Normal(i);

	Array<VectorD> l_normals; l_normals.resize(fluid.l_particles.Size());
#pragma omp parallel for
	for (int i = 0; i < l_normals.size(); i++) l_normals[i] = fluid.l_particles.Normal(i);

	RenderFunc::Write_Vectors_Float<d, real>(frame_dir + "/e_tracker_circles_" + std::to_string(idx), fluid.e_particles.XRef(), e_normals, true);
	RenderFunc::Write_Vectors_Float<d, real>(frame_dir + "/l_tracker_circles_" + std::to_string(idx), fluid.l_particles.XRef(), l_normals, true);
	
	//tracker points
	RenderFunc::Write_Points_Float<d, real>(frame_dir + "/e_tracker_points_" + std::to_string(idx), fluid.e_particles.XRef());
	RenderFunc::Write_Points_Float<d, real>(frame_dir + "/l_tracker_points_" + std::to_string(idx), fluid.l_particles.XRef());

	//BinaryDataIO::Write_Array(frame_dir + "/e_h_bin_" + std::to_string(idx), fluid.temp_scalar);
	BinaryDataIO::Write_Array(frame_dir + "/e_h_bin_" + std::to_string(idx), fluid.e_particles.HRef());
	//BinaryDataIO::Write_Array(frame_dir + "/e_sa_bin_" + std::to_string(idx), fluid.e_particles.SARef());
	//remember to uncomment this
	BinaryDataIO::Write_Array(frame_dir + "/l_h_bin_" + std::to_string(idx), fluid.l_particles.HRef());
	//remember to comment out these
	Array<real> l_h; l_h.resize(fluid.l_particles.Size());
#pragma omp parallel for
	for (int i = 0; i < l_h.size(); i++) l_h[i] = ((real)fluid.l_particles.I(i));
	BinaryDataIO::Write_Array(frame_dir + "/l_origin_bin_" + std::to_string(idx), l_h);
	// done comment
	//BinaryDataIO::Write_Array(frame_dir + "/l_sa_bin_" + std::to_string(idx), fluid.l_particles.SARef());
	//BinaryDataIO::Write_Array(frame_dir + "/l_origin_bin_" + std::to_string(idx), fluid.l_particles.IRef());

//#pragma omp parallel for
	//for (int i = 0; i < fluid.e_particles.Size(); i++) fluid.temp_vector[i] = 10000. * fluid.e_particles.F(i);
	//for (int i = 0; i < fluid.e_particles.Size(); i++) fluid.temp_vector[i] = 10000. * fluid.e_particles.H(i) * fluid.e_particles.Normal(i);
	RenderFunc::Write_Vectors_Float<d, real>(frame_dir + "/e_point_force_" + std::to_string(idx), fluid.e_particles.XRef(), fluid.temp_vector);

	//PointSetFunc::Write_Local_Frames_To_File(frame_dir + "/e_segment_mesh_" + std::to_string(idx), fluid.e_particles, 0.05 * fluid.simulation_scale);

}

template<int d>
void ParticleALEFoamDriver<d>::Save_Eigen_Matrix(std::string fileName, MatrixXd matrix)
{
	const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

	std::ofstream file(fileName);
	if (file.is_open())
	{
		file << matrix.format(CSVFormat);
		file.close();
	}
}

template<int d>
MatrixXd ParticleALEFoamDriver<d>::Load_Eigen_Matrix(std::string fileToOpen)
{
	std::vector<real> matrixEntries;
	std::ifstream matrixDataFile(fileToOpen);
	std::string matrixRowString;
	std::string matrixEntry;
	int matrixRowNumber = 0;
	while (std::getline(matrixDataFile, matrixRowString))
	{
		std::stringstream matrixRowStringStream(matrixRowString);

		while (std::getline(matrixRowStringStream, matrixEntry, ','))
		{
			matrixEntries.push_back(std::stod(matrixEntry));
		}
		matrixRowNumber++;
	}
	if (matrixEntries.size() > 0) {
		return Eigen::Map<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);
	}
	else {
		MatrixXd m;
		return m;
	}
}

template<int d>
void ParticleALEFoamDriver<d>::Save_Snapshot(const int frame)
{
	namespace fs = std::filesystem;
	std::cout << "save snapshot for frame " << frame << "\n";// fluid.Numerical_Check_Force();
	std::string snapshot_dir = output_dir + "/snapshot/" + std::to_string(frame);
	//if (!File::Directory_Exists(snapshot_dir.c_str())) File::Create_Directory(snapshot_dir);
	if (!fs::exists(snapshot_dir)) {
		std::cout << "#     Create snapshot directory: " << snapshot_dir << std::endl;
		fs::create_directories(snapshot_dir);
	}
	for (int i = 0; i < multi_fluid.Num_Regions(); i++) {
		//std::cout << "trying to save snapshot for: " << i << std::endl;
		if (multi_fluid.freeze[i]) continue;
		ParticleALEFilm<d>& fluid = *multi_fluid.regions[i];
		fluid.e_particles.Save_Snapshot(snapshot_dir + "/e_particles_" + std::to_string(i) + ".bin");
		fluid.l_particles.Save_Snapshot(snapshot_dir + "/l_particles_" + std::to_string(i) + ".bin");
		fluid.boundary.Save_Snapshot(snapshot_dir + "/boundary_particles_" + std::to_string(i) + ".bin");
	}
	fluid_3d.particles.Save_Snapshot(snapshot_dir + "/3d_particles.bin");
	Save_Eigen_Matrix(snapshot_dir + "/stresses.csv", multi_fluid.stresses);
}

template<int d>
void ParticleALEFoamDriver<d>::Load_Snapshot(const int frame)
{
	if (frame > 0) {
		if (multi_fluid.bubble_merge_func) {
			multi_fluid.bubble_merge_func(frame - 1);
		}
		if (multi_fluid.bubble_defreeze_func) {
			//std::cout << "running defreeze func!" << std::endl;
			multi_fluid.freeze = multi_fluid.bubble_defreeze_func(frame - 1);
		}
	}

	if (!multi_fluid.recompute_init_COM) {
		multi_fluid.Update_All_Geometry(0.);
		multi_fluid.init_COM = multi_fluid.curr_COM;
		//std::cout << "after load snapshot, init COM is: " << multi_fluid.init_COM << std::endl;
	}

	std::string snapshot_dir = output_dir + "/snapshot/" + std::to_string(frame);
	for (int i = 0; i < multi_fluid.Num_Regions(); i++) {
		if (multi_fluid.freeze[i]) continue;
		std::cout << "trying to load snapshot for: " << i << std::endl;
		ParticleALEFilm<d>& fluid = *multi_fluid.regions[i];
		//std::cout << "1" << std::endl;
		fluid.e_particles.Load_Snapshot(snapshot_dir + "/e_particles_" + std::to_string(i) + ".bin");
		//std::cout << "2" << std::endl;
		fluid.l_particles.Load_Snapshot(snapshot_dir + "/l_particles_" + std::to_string(i) + ".bin");
		//std::cout << "3" << std::endl;
		fluid.boundary.Load_Snapshot(snapshot_dir + "/boundary_particles_" + std::to_string(i) + ".bin");
		//std::cout << "4" << std::endl;
		if (fluid.e_particles.Size()) fluid.Update_E();
		//std::cout << "5" << std::endl;
		if (fluid.l_particles.Size()) fluid.Update_L();
		//std::cout << "6" << std::endl;
		fluid.Reinit_Arrays();
	}
	//std::cout << "7" << std::endl;
	fluid_3d.particles.Load_Snapshot(snapshot_dir + "/3d_particles.bin");
	if (fluid_3d.particles.Size()) {
		fluid_3d.particles.Update();
	}
	//std::cout << "got this far " << std::endl;
	multi_fluid.Reset_Master_List(); //this step is necessary because these std::vectors are actually addresses, we can't recover them
	//std::cout << "got that far " << std::endl;
	//std::cout << "try to load stresses " << std::endl;
	multi_fluid.stresses = Load_Eigen_Matrix(snapshot_dir + "/stresses.csv");
	//std::cout << "done loading stresses " << std::endl;


	if (multi_fluid.recompute_init_COM) {
		multi_fluid.Update_All_Geometry(0.);
		multi_fluid.init_COM = multi_fluid.curr_COM;
		//std::cout << "after load snapshot, init COM is: " << multi_fluid.init_COM << std::endl;
	}
}


template<int d>
void ParticleALEFoamDriver<d>::Load_E_Initialization(ParticleALEFilm<d>& fluid, const std::string e_path)
{
	std::string snapshot_dir = output_dir + "/../" + e_path + "/snapshot/" + std::to_string(-1);
	fluid.e_particles.Load_Snapshot(snapshot_dir + "/e_particles.bin");
}

template<int d>
void ParticleALEFoamDriver<d>::Load_L_Initialization(ParticleALEFilm<d>& fluid, const std::string l_path)
{
	std::string snapshot_dir = output_dir + "/../" + l_path + "/snapshot/" + std::to_string(-1);
	fluid.l_particles.Load_Snapshot(snapshot_dir + "/l_particles.bin");
}


template class ParticleALEFoamDriver<2>;
template class ParticleALEFoamDriver<3>;
