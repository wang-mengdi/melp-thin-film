#include "EulerParticles.h"

template<int d>
void EulerParticles<d>::Diffuse_Scalar(Array<real>& f, const real& dt, const real& coeff, const KernelSPH& kernel, const real radius, const KernelType& kernel_type)
{
	Array<real> diff_v;
	int n = Size();
	diff_v.resize(n);
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		diff_v[i] = coeff * Surface_Laplacian<real>(i, f, kernel, radius, kernel_type);
	}
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		f[i] += diff_v[i] * dt;
	}
}


template<int d>
void EulerParticles<d>::Update_C_Metric_Tensors(const KernelSPH& kernel, const real radius)
{
	////update metric tensor particles.G
#pragma omp parallel for
	for (int i = 0; i < Size(); i++) {
		const auto& nbs = nbs_searcher->Find_Neighbors(X(i), radius);
		VectorT grad_f = VectorT::Zero();
		for (int k = 0; k < nbs.size(); k++) {
			int j = nbs[k];
			real S_j = SA(j);
			VectorT lr_ji = Surface_Relative_Vec(j, i);
			VectorT grad_W = kernel.Grad<d - 1>(lr_ji, radius, KernelType::QUINTIC);
			grad_f += S_j * (X(j) - X(i)).dot(this->Normal(i)) * grad_W;
		}
		this->G(i) = surface.Metric_Tensor(grad_f);
	}
}


template<int d>
void EulerParticles<d>::Update_H_Curvatures(const KernelSPH& kernel, const real radius)
{
#pragma omp parallel for
	for (int i = 0; i < Size(); i++) {
		KH(i) = Surface_Laplacian<real>(i, this->HRef(), kernel, radius, KernelType::QUINTIC);
	}
}


template<int d>
void EulerParticles<d>::Update_C_Curvatures(const KernelSPH& kernel, const real radius)
{
#pragma omp parallel for
	for (int i = 0; i < Size(); i++) {
		std::function<real(const int)> z_func_i = [&](const int idx)->real {return X(idx).dot(this->Normal(i)); };
		KS(i) = Surface_Laplacian<real>(i, z_func_i, kernel, radius, KernelType::QUINTIC);
	}
}


template class EulerParticles<2>;
template class EulerParticles<3>;
