#ifndef __Fluid3DSPH_h__
#define __Fluid3DSPH_h__
#include "SPH3DParticles.h"
#include "ImplicitShape.h"
#include "ALEParams.h"

template<int d>
class Fluid3DSPH
{
	Typedef_VectorDii(d);
	Typedef_MatrixD(d);

public:
	using VectorT = Vector<real, d - 1>;
	real max_vel = 10.;
	VectorD g = VectorD::Zero();
	SPH3DParticles<d> particles;
	ImplicitShape<d>* analytical_boundary;
	SPH3DParams<d> params;
	Array<real> ones;
	real simulation_scale = 1.;
	int mode = 0.; //0 means Akinci approach, suitable for simulating big droplets. 1 means custom approach, suitable for simulating a lot of small droplets
	real rest_rho = 1.;

public:

	void Initialize(void)
	{
		particles.Initialize();
		if (particles.Size() > 0) {
			Array<int> nbs = particles.nbs_searcher->Find_Neighbors(particles.X(0), params.radius);
			std::cout << "num nbs: " << nbs.size() << std::endl;
		}
		if (mode == 0) {
			//params.
		}
		if (mode == 2) {
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				particles.Rho(i) = particles.World_Sum_Value(particles.X(i), particles.MRef(), params.avg_world, params.radius, KernelType::QUINTIC);
			}
			std::cout << "Max Rho at init: " << AuxFunc::Max(particles.RhoRef()) << std::endl;
			std::cout << "Min Rho at init: " << AuxFunc::Min(particles.RhoRef()) << std::endl;
			rest_rho = AuxFunc::Max(particles.RhoRef());
		}
	}

	real Cohesion_Kernel(const real &r, const real &h) const{
		real alpha = 32.0 / (pi * pow(h, 9));
		if (0 <= r && 2 * r <= h) {//0~0.5h
			return alpha * 2 * pow((h - r) * r, 3) - pow(h, 6) / 64.0;
		}
		else if (2 * r > h && r <= h) {
			return alpha * pow((h - r) * r, 3);
		}
		else return 0;
	}

	real Repulsion_Kernel(const real& r, const real& h) const {
		if (r > h/20. && r <= h) {
			real relative_r = r / h;
			return 1. / (relative_r * relative_r);
		}
		else return 0;
	}

	real Adhesion_Kernel(const real& r, const real& h) const {
		real alpha = 0;
		if (r < h && 2 * r > h) {//0~0.5h
			alpha = 0.007 / pow(h, 3.25) * pow(-4*pow(r,2)/h + 6*r -2*h, 1. / 4.);
		}
		return alpha;
	}


	void Add_Particle(const VectorD& X, const VectorD& V, const real& m) {
		int k = particles.Add_Element();
		particles.X(k) = X;
		particles.V(k) = V;
		particles.M(k) = m;
	}

	bool Is_Trivial(void) {
		return particles.Is_Trivial();
	}

	void Update_Max_Velocity(void)
	{
		max_vel = 0.;
		for (int i = 0; i < particles.Size(); i++) {
			max_vel = std::max(max_vel, particles.V(i).norm());
		}
	}

	void Update_Density(void)
	{
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			particles.Rho(i) = particles.World_Sum_Value(particles.X(i), particles.MRef(), params.avg_world, params.radius, KernelType::QUINTIC);
			particles.Vol(i) = particles.M(i) / particles.Rho(i);
		}
	}

	void Update_Normals(void)
	{
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			particles.SN(i) = params.radius * particles.Gradient(i, ones, params.avg_world, params.radius, KernelType::QUINTIC);
			// note: this is not the REAL normal because it is not normalized, neither should it be
		}
	}

	void Update_Cohesion_Force(void)
	{
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			VectorD f_cohesion = VectorD::Zero();
			Array<int> nbs = particles.nbs_searcher->Find_Neighbors(particles.X(i), params.radius);
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				VectorD r_ij = particles.X(i) - particles.X(j);
				f_cohesion += -2./(particles.Rho(i) + particles.Rho(j)) * particles.M(i) * particles.M(j) * r_ij.normalized() * Cohesion_Kernel(r_ij.norm(), params.radius);
			}
			f_cohesion *= params.st_param;
			particles.F(i) += f_cohesion;
		}
	}


	void Update_Repulsion_Force(void)
	{
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			VectorD f_repulsion = VectorD::Zero();
			Array<int> nbs = particles.nbs_searcher->Find_Neighbors(particles.X(i), params.radius);
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				VectorD r_ij = particles.X(i) - particles.X(j);
				f_repulsion += particles.M(i) * particles.M(j) * r_ij.normalized() * Repulsion_Kernel(r_ij.norm(), params.radius);
			}
			f_repulsion *= params.repulsion_param;
			particles.F(i) += f_repulsion;
		}
	}

	void Update_Curvature_Force(void)
	{
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			VectorD f_curvature = VectorD::Zero();
			Array<int> nbs = particles.nbs_searcher->Find_Neighbors(particles.X(i), params.radius);
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				f_curvature += -2. / (particles.Rho(i) + particles.Rho(j)) * particles.M(i) * (particles.SN(i) - particles.SN(j));
			}
			f_curvature *= params.st_param;
			particles.F(i) += f_curvature;
		}
	}


	void Update_WC_Force(void)
	{
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			std::function<real(const int)> p_func_i = [&](const int idx)->real {
				return std::max(0., 100000 * params.pressure_param * (MathFunc::Quick_Pow(particles.Rho(idx)/rest_rho, 2) - 1.)); 
			};
			VectorD f_wc = -particles.Vol(i) * particles.Gradient_Difference(i, p_func_i, params.avg_world, params.radius, KernelType::QUINTIC);
			particles.F(i) += f_wc;
		}
	}

	void Update_Vis_Force(void)
	{
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			VectorD f_vis = VectorD::Zero();
			f_vis = particles.Laplacian(i, particles.VRef(), params.avg_world, params.radius, KernelType::QUINTIC);
			particles.F(i) += params.vis_param * particles.Vol(i) * f_vis;
		}
	}

	void Smooth_Velocity(real strength)
	{
		Array<VectorD> smoothed_V; smoothed_V.resize(particles.Size());
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			bool has_nbs;
			VectorD vel = particles.World_Avg_Value(particles.X(i), particles.VRef(), params.avg_world, params.radius, KernelType::QUINTIC, has_nbs);
			if (has_nbs) smoothed_V[i] = particles.V(i) * (1.-strength) + vel * strength;
			else smoothed_V[i] = particles.V(i);
		}
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			particles.V(i) = smoothed_V[i];
		}
	}

	void Update_Velocity(real dt) {
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			VectorD increment = params.total_force_multiplier * particles.F(i) / particles.M(i) * dt;
			real max_increment = 0.1 * params.dx;
			if (increment.norm() > max_increment /dt) {
				increment = increment.normalized() * max_increment /dt;
			}
			particles.V(i) += increment;
			particles.V(i) += g * dt;
		}
	}

	void XSPH_Smooth_Velocity(const real alpha)
	{
		Array<VectorD> smoothed_V; smoothed_V.resize(particles.Size());
		AuxFunc::Fill(smoothed_V, VectorD::Zero());
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			Array<int> nbs = particles.nbs_searcher->Find_Neighbors(particles.X(i), params.radius);
			VectorD sum = VectorD::Zero();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				VectorD wr_pj = particles.X(j) - particles.X(i);
				real w = params.avg_world.template Weight<d>(wr_pj, params.radius, KernelType::QUINTIC);
				sum += particles.Vol(j) * (particles.V(j) - particles.V(i)) * w;
			}
			smoothed_V[i] = particles.V(i) + alpha * sum;
		}
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			particles.V(i) = smoothed_V[i];
		}
	}

	void Update_Position(real dt) {
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			particles.X(i) += particles.V(i) * dt;
		}
		if (analytical_boundary->Available()) particles.Correct_Position_With_Analytical_Boundary(dt, simulation_scale, *analytical_boundary);
	}

	void Clear(void)
	{
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			particles.F(i) = VectorD::Zero();
		}
	}


	virtual void Advance(const real dt, const real time = 0)
	{	
		if (Is_Trivial()) return;
		std::cout << "3D SPH advance!" << std::endl;
		
		ones.resize(particles.Size());
		std::fill(ones.begin(), ones.end(), 1.);

		particles.Update(); //update nb searcher

		Clear();
		Update_Density();

		if (mode == 0) {
			Update_Normals();
			Update_Curvature_Force();
			Update_Cohesion_Force();
			Update_WC_Force();
			Update_Velocity(dt);
			Smooth_Velocity(0.9);
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				particles.V(i) *= exp(-1. * dt);
			}
		}
		else if (mode == 1) {
			Update_Cohesion_Force();
			Update_Repulsion_Force();
			Update_Velocity(dt);
			Smooth_Velocity(0.33);
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				particles.V(i) *= exp(-1. * dt);
			}
		}
		else if (mode == 2) {
			Update_WC_Force();
			Update_Velocity(dt);
			XSPH_Smooth_Velocity(0.05);
			//Smooth_Velocity(0.0);
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				particles.V(i) *= exp(-1. * dt);
			}
		}
		Update_Position(dt);

		Update_Max_Velocity();
		ArrayFunc::Numerical_Check(particles.XRef(), "3d particeles.x");
	}
};


#endif
