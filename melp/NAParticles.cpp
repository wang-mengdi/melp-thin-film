#include "NAParticles.h"



template<int d>
void NAParticles<d>::Save_Neighbor_Infos(const Array<VectorD>& pos, Array<NeighborInfo>& infos, real radius)
{
	infos.resize(pos.size());
#pragma omp parallel for
	for (int i = 0; i < pos.size(); i++) {
		VectorD x = pos[i];
		Array<int> nbs = nbs_searcher->Find_Neighbors(x, radius);
		infos[i].assign(radius, nbs);
	}
}

template<int d>
void NAParticles<d>::Update_Phis(const ImplicitShape<d>& analytical_boundary)
{
	if (!analytical_boundary.Available()) {
		std::cerr << "SSPHFrame<d>::Update_Phis error: analytical boundary not available\n";
	}
	phis.resize(Size());
	bnd_normals.resize(Size());
#pragma omp parallel for
	for (int i = 0; i < Size(); i++) {
		real phi; VectorD normal;
		analytical_boundary.Nearest_Boundary(X(i), phi, normal);
		phis[i] = -phi;
		bnd_normals[i] = normal;
	}
}

template<int d>
void NAParticles<d>::Correct_Velocity_With_Analytical_Boundary(const real dt, const real scale, const ImplicitShape<d>& analytical_boundary)
{
	if (!analytical_boundary.Available()) {
		std::cerr << "SSPHFrame<d>::Correct_Velocity_With_Analytical_Boundary error: analytical boundary not available\n";
	}
	Update_Phis(analytical_boundary);
#pragma omp parallel for
	for (int i = 0; i < Size(); i++) {

		real cushion = 0.01 * scale;
		real phi = phis[i];
		VectorD normal = bnd_normals[i];
		if (phi > -cushion) {
			//V(i) *= 0;
			VectorD normal_velocity = V(i).dot(normal) * normal;
			VectorD tangential_velocity = V(i) - normal_velocity;
			if (phi < 0) {
				V(i) = tangential_velocity + (-phi / cushion) * normal_velocity;
			}
			else {
				real normal_rate = (-phi) / dt;
				V(i) = tangential_velocity + normal_rate * normal;
			}
		}
	}
}

template<int d>
void NAParticles<d>::Correct_Position_With_Analytical_Boundary(const real dt, const real scale, const ImplicitShape<d>& analytical_boundary)
{
	if (!analytical_boundary.Available()) {
		std::cerr << "SSPHFrame<d>::Correct_Velocity_With_Analytical_Boundary error: analytical boundary not available\n";
	}
	Update_Phis(analytical_boundary);
#pragma omp parallel for
	for (int i = 0; i < Size(); i++) {
		VectorD vel = VectorD::Zero();
		real cushion = 0.01 * scale;
		real phi = phis[i];
		//std::cout << "phi: " << phi << std::endl;
		VectorD normal = bnd_normals[i];
		if (phi > -cushion) {
			VectorD normal_velocity = V(i).dot(normal) * normal;
			VectorD tangential_velocity = V(i) - normal_velocity;
			if (phi < 0) {
				vel = VectorD::Zero();
			}
			else {
				real normal_rate = phi;
				vel = normal_rate * normal;
			}
			X(i) += vel;
			//V(i) = VectorD::Zero();
			V(i) = tangential_velocity;
			VectorD reflected_normal_velocity = VectorD::Zero();
			if (dt >= 1.e-8) reflected_normal_velocity = vel / dt;
			if (normal_velocity.dot(normal) > 0) V(i) += normal_velocity; // give back the normal velocity if it is not into the surface
			else V(i) += reflected_normal_velocity;
		}
	}
}


template<int d>
void NAParticles<d>::Print_Fastest_Particle(void)
{
	int max_idx = -1; real max_vel = 0.0;
	for (int i = 0; i < Size(); i++) {
		if (V(i).norm() > max_vel) {
			max_idx = i;
			max_vel = V(i).norm();
		}
	}
	if (max_idx != -1) {
		Print_Dynamics(max_idx);
	}
}


template class NAParticles<2>;
template class NAParticles<3>;

NeighborInfo::NeighborInfo()
{
	radius = -1;
	data.clear();
}

void NeighborInfo::assign(real r, const Array<int>& nbs)
{
	radius = r;
	data = nbs;
}

int NeighborInfo::size(void)const
{
	return data.size();
}

int NeighborInfo::operator[](int k)const
{
	return data[k];
}
