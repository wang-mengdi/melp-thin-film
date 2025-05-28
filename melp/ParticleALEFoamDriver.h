//////////////////////////////////////////////////////////////////////////
// Particle ALE film driver 
// Copyright (c) (2020-), Yitong Deng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __ParticleALEFoamDriver_h__
#define __ParticleALEFoamDriver_h__
#include <iostream>
#include <filesystem>
#include <fstream>
#include "Driver.h"
#include "PointSetFunc.h"
#include "AuxFunc.h"
#include "PerlinNoise.hpp"
#include "BoundaryParticles.h"
#include "ParticleALEFilm.h"
#include "ParticleALEFoam.h"
#include "Fluid3DSPH.h"
#include "ParticleALEFilmDriver.h"

template<int d> class ParticleALEFoamDriver : public Driver
{
	Typedef_VectorDii(d);
	Typedef_MatrixD(d);
	using Base=Driver;

public:
	using VectorT=Vector<real,d-1>;								////tangential vector type
	using VectorTi=Vector<int,d-1>;								////tangential vector int

	ParticleALEFoam<d> multi_fluid;
	Fluid3DSPH<d> fluid_3d;

	virtual void Advance_To_Target_Time(const real target_time)
	{
		bool done = false; 
		for (int substep = 1; !done; substep++) {
			real max_vel = 0.;
			if (multi_fluid.Is_Trivial()) {
				multi_fluid.Defreeze_Bubbles();
			}
			if (!multi_fluid.Is_Trivial()) { 
				max_vel = multi_fluid.max_vel; 
				if (verbose) std::cout << "[Driver]: multi_fluid max vel:" << multi_fluid.max_vel << std::endl;
			}
			if (!fluid_3d.Is_Trivial()) { 
				max_vel = std::max<real>(max_vel, fluid_3d.max_vel); 
				if (verbose) std::cout << "[Driver]: fluid 3d max vel:" << fluid_3d.max_vel << std::endl;
			}
			real vel = std::max(CFL() * multi_fluid.e_dx * frame_rate, max_vel);
			real dt = CFL() * multi_fluid.e_dx / vel;
			if (max_iter_per_frame > 0) {
				dt = std::max(dt, 1.0 / (frame_rate * max_iter_per_frame));
			}
			if (time + dt >= target_time) { dt = target_time - time; done = true; }
			else if (time + 2 * dt >= target_time) { dt = (real).5 * (target_time - time); }
			Advance_One_Time_Step(dt, time, substep-1);
			time += dt;
		}
	}

	virtual void Advance_One_Time_Step(const real dt,const real time, const int iter)
	{
		double begin_time = omp_get_wtime();
		std::cout << std::defaultfloat;

		multi_fluid.curr_frame = current_frame;
		// remember to also notify the individual fluids as well
		for (int i = 0; i < multi_fluid.Num_Regions(); i++) {
			ParticleALEFilm<d>& fluid = *multi_fluid.regions[i];
			fluid.curr_frame = current_frame;
		}

		if (!multi_fluid.Is_Trivial()) multi_fluid.Advance(dt, time, iter);
		if (!fluid_3d.Is_Trivial()) fluid_3d.Advance(dt, time);

		double end_time = omp_get_wtime();
		int iter_per_frame = (int)(1.0 / frame_rate / dt + 0.5);
		if (dt < 1e-8) iter_per_frame = 0;
		std::cout << "["
			<< std::setw(6) << iter_per_frame << " Iterations Per Frame]     ... "
			<< std::setw(7) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << end_time - begin_time << "s used, current time: " << time << std::endl;
		std::cout << std::defaultfloat;
	}
	
	virtual void Initialize()
	{
		multi_fluid.fluid_3d = &fluid_3d;
		cfl=.1; //default CFL
		switch (test) {
		case 1:	Case_1(); break; //Three Bubbles
		}
		if (first_frame > 0) {
			int last_saved = int(first_frame / snapshot_stride) * snapshot_stride;
			std::cout << "Run from snapshot frame " << last_saved << ".\n";
			current_frame = last_saved;
			time = Time_At_Frame(current_frame);
			Load_Snapshot(last_saved);
			std::cout << "Snapshot loaded.\n";
		}
	}

	void Case_1(void);

	real Temperature(VectorD pos);

	void Setup_Bubble(ParticleALEFilm<d>& fluid, const real radius, const VectorD center, const VectorD velocity, const real thickness_multiplier = 1., const real fineness_multiplier = 1., const std::string l_init = "", const std::string e_init = "");

	virtual void Write_Output_Files(const int frame);
	virtual void Write_Output_Files_Region(const int frame, const int idx);

	void Save_Snapshot(const int frame);
	void Load_Snapshot(const int frame);
	void Load_E_Initialization(ParticleALEFilm<d>& fluid, const std::string e_path);
	void Load_L_Initialization(ParticleALEFilm<d>& fluid, const std::string l_path);

	void Save_Eigen_Matrix(std::string fileName, MatrixXd matrix);
	MatrixXd Load_Eigen_Matrix(std::string fileToOpen);
};

	
#endif
