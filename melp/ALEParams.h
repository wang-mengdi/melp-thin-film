//////////////////////////////////////////////////////////////////////////
// Parameters for particle ALE film
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __ALEParams_h__
#define __ALEParams_h__
#include "Kernels.h"
#include "FileInit.h"

//control neighbor search radius, kernels etc
template<int d>
class NeighborParams :public SimulatorParams {
public:
	real e_dx = 1.0;
	real sph_radius = 1.0;
	real interp_radius = 1.0;
	KernelSPH kernel;//to calculate values on Euler Particles (used by L particles, since they don't have reference frame)
	NeighborParams() :SimulatorParams("neighbor") {}
	void Initialize(real _e_dx, real ee, real el) {
		e_dx = _e_dx;
		sph_radius = ee * e_dx;
		interp_radius = el * e_dx;
	}
	NeighborParams(real _e_dx, real ee, real el) :SimulatorParams("neighbor") {
		Initialize(_e_dx, ee, el);
	}
};

class ChemicalParams : public SimulatorParams {
public:
	real R = 8.3144598; //gas constant
	real rho_0 = 997.; //density of water
	real T_0 = 300.; //room temperature
	real expansion = 0.3; //expansion param
	real drag_C = 1.e-3; //drag param, see paper 
	real heat_capacity_param = 1.;
	ChemicalParams() :SimulatorParams("chemical") {}
};

class NumericalParams : public SimulatorParams {
public:
	real Velocity_Decay = .01;
	real XSPH_V_Tang = 0.999;
	real XSPH_V_Norm = 0.1;
	int XSPH_KS_Passes = 0;
	int XSPH_V_Passes = 0;
	real XSPH_V_Passes_Strength = 0.999;
	NumericalParams() :SimulatorParams("numerical") {}
};

//dedicated for 3d sph
template<int d>
class SPH3DParams :public SimulatorParams {
public:
	real dx = 1.0;
	int np_on_h = 3;
	real radius = 1.0;
	real st_param = 3.;
	real repulsion_param = 10000.;
	real vis_param = 0.0002;
	real pressure_param = 0.0001;
	real total_force_multiplier = 1.;
	KernelSPH avg_world;//3D Kernel
	SPH3DParams() :SimulatorParams("sph_3d") {}
	SPH3DParams(real _dx, int _np) :SimulatorParams("neighbor") {
		dx = _dx;
		np_on_h = _np;
		radius = dx * np_on_h;
	}
};

#endif