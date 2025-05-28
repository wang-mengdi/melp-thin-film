//////////////////////////////////////////////////////////////////////////
// Analytical boundary for particle system
// Copyright (c) (2018-), Mengdi Wang, Yitong Deng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "Common.h"
#include "GeometryPrimitives.h"
#include "SPHParticles.h"

//Implicit shape composed of several geometries and complementary geometries.
//In other words, the combination of some ImplicitGeometry<d> (s).
//It's like: shape1 AND shape2 AND shape3 AND...
//Follow the common level-set convention
//phi>0: outside
//phi<0: inside
template<int d>
class ImplicitShape {
	Typedef_VectorDii(d);
public:
	//(shape,sgn) pairs
	//sgn==1 means the shape itself
	//sgn==-1 means the complementary shape
	Array<std::pair<std::shared_ptr<ImplicitGeometry<d> >, int> > shapes;//shape, sign
	
public:
	////Basic Info
	ImplicitShape(void) { shapes.clear(); }
	int Size(void) const { return shapes.size(); }
	bool Available() const { return Size() > 0; } // if some boundary is active

	////Building Functions
	void Add_Obstacle(std::shared_ptr<ImplicitGeometry<d> > obj) { shapes.push_back(std::make_pair(obj, 1)); }//outside of this
	void Add_Bounding(std::shared_ptr<ImplicitGeometry<d> > obj) { shapes.push_back(std::make_pair(obj, -1)); }//inside of this
	ImplicitShape<d> operator - (void)const;
	ImplicitShape<d>& operator += (std::shared_ptr<ImplicitGeometry<d> > obj) { Add_Obstacle(obj); return *this; }
	ImplicitShape<d>& operator + (std::shared_ptr<ImplicitGeometry<d> > obj) { Add_Obstacle(obj); return *this; }
	ImplicitShape<d>& operator - (std::shared_ptr<ImplicitGeometry<d> > obj) { Add_Bounding(obj); return *this; }
	ImplicitShape<d>& operator -= (std::shared_ptr<ImplicitGeometry<d> > obj) { Add_Bounding(obj); return *this; }
	ImplicitShape<d>& operator + (const ImplicitShape<d>& shape_b);
	ImplicitShape<d>& operator - (const ImplicitShape<d>& shape_b);

	////Some Basic Shapes
	static ImplicitShape<d> Empty(void) { return ImplicitShape<d>(); }
	static ImplicitShape<d> Bounding_Box(const VectorD& center, const VectorD& lengths);
	static ImplicitShape<d> Bounding_Sphere(const VectorD& center, const real radius);

	////Geometry computations
	//if the shape is empty, return false
	//the nearest boundary (or the minimum phi) is saved in dis and normal
	bool Nearest_Boundary(const VectorD& pos, real& phi, VectorD& normal) const;
	bool Inside(const VectorD& pos) const;
	real Phi(const VectorD& pos) const;
	VectorD Normal(const VectorD& pos,real dx=(real)1e-3)const;
};

//////////////////////////////////////////////////////////////////////////
//// implicit geometry wrapper to be added to implicit shape
//////////////////////////////////////////////////////////////////////////
template<class Type,int d> class ImplicitGeometryWrapper : public ImplicitGeometry<d>
{Typedef_VectorDii(d);
public:
	std::shared_ptr<Type> object;
	ImplicitGeometryWrapper(std::shared_ptr<Type> _object):object(_object){}
	ImplicitGeometryWrapper<Type,d>& operator=(const ImplicitGeometryWrapper<Type,d>& copy){object=copy.object;return *this;}
	ImplicitGeometryWrapper(const ImplicitGeometryWrapper<Type,d>& copy){*this=copy;}

	virtual real Phi(const VectorD& pos) const 
	{
		return object->Phi(pos);
	}

	virtual VectorD Normal(const VectorD& pos) const
	{
		return object->Normal(pos);
	}
};

