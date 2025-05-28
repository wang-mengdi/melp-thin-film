//////////////////////////////////////////////////////////////////////////
// Particle ALE film driver 
// Copyright (c) (2020-), Yitong Deng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __ParticleALEFilmDriver_h__
#define __ParticleALEFilmDriver_h__
#include "Driver.h"
#include "PointSetFunc.h"
//#include "TinyObjLoader.h"
#include "AuxFunc.h"
#include "PerlinNoise.hpp"
#include "BoundaryParticles.h"
#include "ParticleALEFilm.h"
#include "Fluid3DSPH.h"

template<int d> class ParticleALEFilmDriver : public Driver
{
	Typedef_VectorDii(d);
	Typedef_MatrixD(d);
	using Base=Driver;
public:
	using VectorT=Vector<real,d-1>;								////tangential vector type
	using VectorTi=Vector<int,d-1>;								////tangential vector int

	ParticleALEFilm<d> fluid;
	Fluid3DSPH<d> fluid_3d;

	int iomode = 0;
	int e_repack_steps = 0;
	int l_repack_steps = 0;
	std::string e_init;
	std::string l_init;
	std::string param_file_name;

	virtual void Run()
	{
		Timer timer;
		if (current_frame == 0) Write_Output_Files(current_frame);
		while (current_frame < last_frame) {
			if (fluid.bursting && fluid.e_particles.Size() > 0) {
				frame_rate = fluid.bursting_frame_rate;
			}
			else {
				frame_rate = fluid.frame_rate;
			}
			timer.Reset();
			Advance_To_Target_Time(time + 1./frame_rate);
			current_frame++;
			double frame_time = timer.Elapse(PhysicalUnits::s);
			double total_time = timer.Total_Elapsed(PhysicalUnits::s);
			Info("Frame time: {:.2f}s/{:.2f}s, ETA {:.2f}s", frame_time, total_time, total_time / current_frame * (last_frame - current_frame));

			Write_Output_Files(current_frame);
		}
	}

	virtual void Advance_To_Target_Time(const real target_time)
	{
		bool done = false;
		for (int substep = 1; !done; substep++) {
			if (fluid.bursting && fluid.e_particles.Size() > 0) {
				frame_rate = fluid.bursting_frame_rate;
			}
			else {
				frame_rate = fluid.frame_rate;
			}
			real max_vel = 0.;
			if (!fluid.Is_Trivial()) { 
				max_vel = fluid.max_vel; 
				if (verbose) std::cout << "[Driver]: fluid max vel:" << fluid.max_vel << std::endl;
				if (verbose) std::cout << "[Driver]: fluid E max vel:" << fluid.max_vel_e << std::endl;
				if (verbose) std::cout << "[Driver]: fluid L max vel:" << fluid.max_vel_l << std::endl;
				if (verbose) std::cout << "[Driver]: fluid 3D max vel:" << fluid.max_vel_3d << std::endl;
			}
			if (!fluid_3d.Is_Trivial()) { 
				max_vel = std::max<real>(max_vel, fluid_3d.max_vel); 
				if (verbose) std::cout << "[Driver]: fluid 3d max vel:" << fluid_3d.max_vel << std::endl;
			}
			real vel = std::max(CFL() * fluid.neighbor_params.e_dx * frame_rate, max_vel);
			real dt = CFL() * fluid.neighbor_params.e_dx / vel;
			if (max_iter_per_frame > 0) {
				dt = std::max(dt, 1.0 / (frame_rate * max_iter_per_frame));
			}
			if (time + dt >= target_time) { dt = target_time - time; done = true; 
			}
			else if (time + 2 * dt >= target_time) { dt = (real).5 * (target_time - time); 			
			}
			double begin_time = omp_get_wtime();
			std::cout << std::defaultfloat;
			Advance_One_Time_Step(dt, time, substep-1);
			double end_time = omp_get_wtime();
			int iter_per_frame = substep + (int)((target_time-time)/dt + 0.5);
			if (dt < 1e-8) iter_per_frame = 0;
			std::cout << "["
				<< std::setw(6) << iter_per_frame << " Iterations Per Frame]     ... "
				<< std::setw(7) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << end_time - begin_time << "s used, current time: " << time << std::endl;
			std::cout << std::defaultfloat;
			
			time += dt;
		}
	}

	virtual void Advance_One_Time_Step(const real dt,const real time, const int iter)
	{
		fluid.curr_frame = current_frame;

		if (!fluid_3d.Is_Trivial()) fluid_3d.Advance(dt, time);
		if (!fluid.Is_Trivial()) fluid.Advance(dt, time, iter);
	}
	
	virtual void Initialize()
	{
		fluid.fluid_3d = &fluid_3d;
		cfl=.1;
		switch (test) {
		case 1:	Case_1(); break; //Bubble
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

	static real Perlin_Noise(VectorD pos, std::uint32_t seed, real perlin_freq, real perlin_scale, std::uint32_t octaves=4);
	real Temperature(VectorD pos, real y_lowest, real y_highest, real strength = 1., real perturbation_ratio = 1.);
	void Case_1(void);

	void Write_Particle_Basics(const NAParticles<d>& particles, const std::string& frame_dir, const std::string& prefix);
	virtual void Write_Output_Files(const int frame);

	void Save_Snapshot(const int frame);
	void Load_Snapshot(const int frame);
	void Load_E_Initialization(const std::string e_path);
	void Load_L_Initialization(const std::string l_path);
};

template<int d>
real ParticleALEFilmDriver<d>::Perlin_Noise(VectorD pos, std::uint32_t seed, real perlin_freq, real perlin_scale, std::uint32_t octaves)
{
	/// 
	/// Return a perlin noise value.
	/// 
	/// pos : position of a particle.
	/// seed : A random seed.
	/// perlin_freq : The larger this parameter, the more frequent the noise.
	/// perlin_scale : The scale of offset from the original position. 
	/// octaves : Number of octaves in perlin noise, default value is 4.
	/// noise_value : Return value of perlin noise generator, range in [-1, 1].
	/// 
	const siv::PerlinNoise perlin(seed);
	real noise_value;
	if constexpr (d == 3) {
		noise_value = perlin.octave3D(pos[0] * perlin_freq, pos[1] * perlin_freq, pos[2] * perlin_freq, octaves);
	}
	else {
		noise_value = perlin.octave2D(pos[0] * perlin_freq, pos[1] * perlin_freq, octaves);
	}
	return noise_value * perlin_scale;
}
	
#endif
