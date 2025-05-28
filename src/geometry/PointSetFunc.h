//////////////////////////////////////////////////////////////////////////
// Point set functions
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __PointSetFunc_h__
#define __PointSetFunc_h__
#include <random>
#include <fstream>
#include <iostream>
#include "Mesh.h"
#include "MeshFunc.h"
#include "AuxFunc.h"
#include "Constants.h"
#include "GeometryParticles.h"
#include "RandomNumber.h"
#include "MacGrid.h"
#include "GeometryPrimitives.h"
#include "RenderFunc.h"


namespace PointSetFunc
{
	using namespace AuxFunc;

	//////////////////////////////////////////////////////////////////////////
	////2D points initialization
	real Circle_Point_Cloud(const Vector2& center, const real R, const int num, Array<Vector2> &points);//return min dist
	real Initialize_Circle_Points(const Vector2& c, const real R, const int p_num, GeometryParticles<2>& particles);
	real Initialize_Oval_Points(const Vector2& c, const real R, const real a, const real b, const int p_num, GeometryParticles<2>& particles);
	int Initialize_Circle_Points(const Vector2& c, const real R, const real dx, GeometryParticles<2>& particles);
	void Initialize_Rectangle_Points(const int nx, const int ny, const real dx, const Vector2& start, GeometryParticles<2>& particles);
	void Initialize_Round_Corner_Rectangle_Points(const int nx, const int ny, const real dx, const real r, const Vector2& start, GeometryParticles<2>& particles);
	void Initialize_Segment_Points(const Vector2& v1, const Vector2& v2, int N, GeometryParticles<2>& particles);
	void Initialize_Curve_Points(const Vector2& c,const real R,const real theta,const int N,GeometryParticles<2>& particles);
	void Extend_Hexagon_Packing(Points<2>& points, const real &R, const real &min_dist);//there must be exactly 1 element in points, and generate all hexagon packing within distance R around it

	//////////////////////////////////////////////////////////////////////////
	////3D points initialization
	real Ring_Point_Cloud(const Vector3& center, const Vector3& direction, const real R, const int num, Array<Vector3>& points);
	real Initialize_Sphere_Points(const Vector3& c, const real R, const int sub, GeometryParticles<3>& particles, const real randomness = 0., std::function<bool(Vector3)> filter_func = nullptr);
	real Initialize_Sphere_Points_Regular(const Vector3& center, const real R, const int num_pts, GeometryParticles<3>& pts);
	void Initialize_Sphere_Points_Random(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles);
	void Initialize_Sphere_Points_Random_Noise(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles, std::function<real(Vector3)> noise_func);
	void Initialize_Half_Sphere_Points_Random(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles, std::function<real(Vector3)> noise_func);
	void Initialize_Circle_Points_Hexagon(GeometryParticles<3>& points, const real& R, const real& min_dist, const real& random);
	void Initialize_Circle_Points_Hexagon_Mask(GeometryParticles<3>& points, const real& R_out, const real& R_in, const real& min_dist, const real& random);
	std::vector<int> Initialize_Circle_Points_Grid(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles, bool with_boundary);
	void Initialize_Circle_Rim_Points(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles);
	std::vector<int> Initialize_Circle_Points_Grid_Random(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles, real randomness = 1., bool with_boundary = true);
	std::vector<int> Initialize_Circle_Points_Grid_Random_Mask(const Vector3& c, const real R_out, const real R_in, const real dx, GeometryParticles<3>& particles, real randomness = 1.);
	std::vector<int> Initialize_Circle_Points_Grid_Noise(const std::function<real(Vector3)> noise_func, const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles, real randomness = 1., bool with_boundary = true);
	std::vector<int> Initialize_Circle_Points_Random(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles);
	void Initialize_Circle_Points_Random2(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles);
	std::vector<int> Initialize_Catenoid_Points(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles);
	void Initialize_Catenoid_Points2(const Vector3& c, const real R, const real separation, const real dx, GeometryParticles<3>& particles, real random, std::vector<int>* is_boundary);
	std::vector<int> Initialize_Catenoid_Points_Random(const Vector3& c, const real R, const real separation, const real dx, GeometryParticles<3>& particles);
	void Initialize_Rectangle_Points(const int nx, const int ny, const real dx, const Vector3& start, const Vector3& normal, GeometryParticles<3>& particles);
	void Initialize_Box_Points(const int nx, const int ny, const int nz, const real dx, const Vector3& start, GeometryParticles<3>& particles, const real randomness = 0.);
	std::vector<int> Initialize_Lattice_Points(const Vector3& domain_center, const Vector2i& counts, Vector3 k1, Vector3 k2, const real dx, GeometryParticles<3>& particles);//nx*ny cells, not nodes
	std::vector<int> Initialize_Lattice_Points_Mask(const Vector3& domain_center, const Vector2i& counts, const int padding_num, Vector3 k1, Vector3 k2, const real dx, GeometryParticles<3>& particles);
	std::vector<int> Initialize_Lattice_Points2(const Vector3& domain_min, const Vector2i& counts, Vector3 k1, Vector3 k2, const real dx, GeometryParticles<3>& particles);
	std::vector<int> Initialize_Catenoid_Points(const Vector3& center, const real R, const int nr, const real height, const int nh, GeometryParticles<3>& particles);
	MacGrid<3> Extend_FCC_Packing(Points<3>& points, const real& R, const real& min_dist);//there must be exactly 1 element in points, and generate all FCC packing within distance R around it

	void Read_Points_From_CSV(std::string locations, std::string normals, GeometryParticles<3>& particles);

	//////////////////////////////////////////////////////////////////////////
	////Universal points initialization
	template<int d>	void Extend_Inside(Points<d>& points, const real& dx, const real& max_radius, const ImplicitGeometry<d>* geom = nullptr);
	template<int d> void Extend_Identical_Particles(Points<d>& points, const Array<Vector<real, d> >& positions);

	//////////////////////////////////////////////////////////////////////////
	////File IO
	template<int d> void Write_Local_Frames_To_File(const std::string file_name, const GeometryParticles<d>& points, const real scale = (real).02)
	{
		SegmentMesh<3> s3;
		for (int i = 0; i < points.Size(); i++) {
			if (points.I(i) == -1)continue;
			for (int j = 0; j < d; j++) {
				Vector<real, d> e = points.E(i).col(j);
				s3.Vertices().push_back(V<3>(points.X(i)));
				Vector<real, d> s_end = points.X(i) + e * scale;
				s3.Vertices().push_back(V<3>(s_end));
				int s = (int)s3.Vertices().size();
				s3.Elements().push_back(Vector2i(s - 2, s - 1));
			}
		}
		s3.Write_To_File_3d(file_name);
	}

	template<int d> void Write_Tracker_Circles_To_File(const std::string file_name, const GeometryParticles<d>& points)
	{
		int pn = points.Size(); Array<Vector<real, d> > normals(pn);
		for (int i = 0; i < pn; i++) {
			normals[i] = points.Normal(i);
			//normals[i] = Fit_Vector<real, d>(0, 0, 1).normalized();
		}
		RenderFunc::Write_Vectors_Float<d, real>(file_name, points.XRef(), normals, true);
	}

	template<int d> void Write_Tracker_Circles_To_File2(const std::string file_name, const Array<Vector<real, d>> normals, const Array<Vector<real, d>> locations)
	{
		RenderFunc::Write_Vectors_Float<d, real>(file_name, locations, normals);
	}
};
#endif
