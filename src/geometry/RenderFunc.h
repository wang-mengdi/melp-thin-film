//////////////////////////////////////////////////////////////////////////
// Assemble data for opengl_viewer. Especially for Points.
// Copyright (c) (2018-),Mengdi Wang
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "GeometryParticles.h"
#include "FaceField.h"
#include "MacGrid.h"
#include "BoundaryCondition.h"
#include "Json.h"
#include "LevelSet.h"

//All these functions directly output to opengl_viewer. 
//So all 2D needs to be converted to 3D here.
namespace RenderFunc {
	////Eulerian
	//Write Field
	template<int d,class T> void Write_Field(std::string file_name,const Field<T,d>& field);

	//write a scalar field stored at face centers on a MAC grid to a vector field stored on cell centers
	template<int d> void Write_Face_Field_On_Cells(const std::string& file_name,const MacGrid<d>& mac_grid,const FaceField<real,d>& val);
	
	//write a scalar field stored at face centers on a MAC grid to a vector field with each face having an axis-aligned vector
	template<int d> void Write_Face_Field_On_Faces(const std::string& file_name,const MacGrid<d>& mac_grid,const FaceField<real,d>& face_field);

	////Lagrangian
	template<int d,class T> void Write_Points_Float(std::string file_name,const Array<Vector<T,d> >& X);
	//with as_ends==false,write X--X+V. otherwise,write segments X--V. 
	template<int d,class T> void Write_Vectors_Float(std::string file_name,const Array<Vector<T,d> >& X,const Array<Vector<T,d> >& V,bool as_ends = false);
	//write every scalar as a point. it's on the normal direction of particle i,and the distance between it and particle i is arr[i]*scale.
	template<int d> void Write_Scalars_As_Points_Float(std::string file_name,const GeometryParticles<d>& particles,const Array<real>& arr,const real scale = 1.0);
	//each segment originates from xs[i],with direction normals[i] and length len_arr[i]*scale
	template<int d> void Write_Customized_Segments_Float(std::string file_name,const Array<Vector<real,d> >& xs,const Array<Vector<real,d> >& normals,const Array<real>& len_arr,const real scale = 1.0);

	////boundary condition
	template<int d> void Write_Dirichlet_Boundary_Conditions(std::string file_name,const MacGrid<d>& mac_grid,const BoundaryConditionMacGrid<d>& bc);
	template<int d> void Write_Neumann_Boundary_Conditions(std::string file_name,const MacGrid<d>& mac_grid,const BoundaryConditionMacGrid<d>& bc);
	
	////Marching cubes
	template<int d> void Write_LevelSet_Mesh(std::string file_name,LevelSet<d>& levelset);
}