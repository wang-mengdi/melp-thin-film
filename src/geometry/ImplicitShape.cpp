#include "ImplicitShape.h"

template<int d>
ImplicitShape<d> ImplicitShape<d>::operator-(void)const
{
	ImplicitShape<d> b = *this;
	for (auto i = 0; i < b.shapes.size(); i++) {
		b.shapes[i].second *= -1;
	}
	return b;
}

template<int d>
ImplicitShape<d>& ImplicitShape<d>::operator+(const ImplicitShape<d>& shape_b)
{
	shapes.insert(shapes.end(), shape_b.shapes.begin(), shape_b.shapes.end());
	return *this;
}

template<int d>
ImplicitShape<d>& ImplicitShape<d>::operator-(const ImplicitShape<d>& shape_b)
{
	return (*this) + (-shape_b);
}

template<int d>
ImplicitShape<d> ImplicitShape<d>::Bounding_Box(const VectorD& center, const VectorD& lengths)
{
	ImplicitShape<d> shape = ImplicitShape<d>::Empty();
	for (int axis = 0; axis < d; axis++) {
		VectorD e = VectorD::Unit(axis);
		VectorD offset = e * lengths[axis] * 0.5;
		auto upper_plane = std::make_shared<Plane<d>>(-e, center + offset);
		auto lower_plane = std::make_shared<Plane<d>>(e, center - offset);
		shape += upper_plane;
		shape += lower_plane;
	}
	return shape;
}

template<int d>
ImplicitShape<d> ImplicitShape<d>::Bounding_Sphere(const VectorD& center, const real radius)
{
	auto sphere = std::make_shared<Sphere<d>>(center, radius);
	return ImplicitShape<d>::Empty() - sphere;
}

template<int d>
bool ImplicitShape<d>::Nearest_Boundary(const VectorD& pos, real& phi, VectorD& normal) const
{
	phi = std::numeric_limits<real>::max();
	normal = VectorD::Zero();
	if (shapes.empty()) return false;
	int min_idx = -1;
	for (int i = 0; i < shapes.size(); i++) {
		auto obj = shapes[i].first;
		int sgn = shapes[i].second;
		real obj_phi = obj->Phi(pos) * sgn;
		if (obj_phi < phi) {
			phi = obj_phi;
			min_idx = i;
		}
	}
	normal = shapes[min_idx].first->Normal(pos) * shapes[min_idx].second;
	return true;
}

template<int d>
bool ImplicitShape<d>::Inside(const VectorD& pos)const {
	bool inside = true;
	for (auto pair : shapes) {
		auto obj = pair.first;
		auto sgn = pair.second;
		if (sgn == 1) inside = inside && obj->Inside(pos);
		else inside = inside && (!obj->Inside(pos));
	}
	return inside;
}

template<int d>
real ImplicitShape<d>::Phi(const VectorD& pos) const
{
	int my_sgn = 1;
	real abs_phi = std::numeric_limits<real>::max();
	for (int i = 0; i < shapes.size(); i++) {
		auto obj = shapes[i].first;
		int shape_sgn = shapes[i].second;
		real obj_phi = obj->Phi(pos) * shape_sgn;
		if (obj_phi < 0) my_sgn = -1;
		abs_phi = std::min(abs_phi, fabs(obj_phi));
	}
	return abs_phi * my_sgn;
}


template<int d>
Vector<real,d> ImplicitShape<d>::Normal(const VectorD& pos,real dx) const
{
	VectorD normal;
	for(int i=0;i<d;i++){
		VectorD pos0=pos-VectorD::Unit(i)*dx;real phi0=Phi(pos0);
		VectorD pos1=pos+VectorD::Unit(i)*dx;real phi1=Phi(pos1);
		normal[i]=(phi1-phi0)/((real)2*dx);}
	real length=normal.norm();
	if(length>(real)0)normal/=length;return normal;
}

template class ImplicitShape<2>;
template class ImplicitShape<3>;