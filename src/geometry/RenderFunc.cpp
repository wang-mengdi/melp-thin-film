//////////////////////////////////////////////////////////////////////////
// Assemble data for opengl_viewer. Especially for Points.
// Copyright(c)(2018-),Mengdi Wang
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "RenderFunc.h"
#include "AuxFunc.h"
#include "Interpolation.h"
#include "MarchingCubes.h"

namespace RenderFunc {

template<int d,class T> void Write_Field(std::string file_name,const Field<T,d>& field)
{
	if constexpr(d==1||d==3)
		File::Write_Binary_To_File(file_name,&field);
	else {
		////PBG: needs conversion when d==2,works for scalar field only
		Field<T,3> field3;
		Dim_Conversion<T,d,3>(field,field3);
		File::Write_Binary_To_File(file_name,field3);}
}

#define Inst_Helper(T)														\
template void Write_Field<2,T>(std::string,const Field<T,2>&);			\
template void Write_Field<3,T>(std::string,const Field<T,3>&);
Inst_Helper(float);	Inst_Helper(double);Inst_Helper(int);Inst_Helper(ushort);
#undef Inst_Helper

template<int d> void Write_Face_Field_On_Cells(const std::string& file_name,const MacGrid<d>& mac_grid,const FaceField<real,d>& val)
{
	Field<Vector<real,d>,d> grid_v;
	Field<real,d> grid_a;
	grid_v.Resize(mac_grid.grid.cell_counts,Vector<real,d>::Zero());
	grid_a.Resize(mac_grid.grid.cell_counts,0.0);

	Interpolation<d> intp(mac_grid);
	intp.Interpolate_Faces_To_Cells(val,grid_v);
	iterate_cell(iter,mac_grid.grid){
		Vector<int,d> cell=iter.Coord();
		Vector<real,d> cell_val=grid_v(cell);
		real a=0;for(int k=0;k<d;k++)a+=cell_val[k];
		grid_a(cell)=a/d;}
	Write_Field<d,real>(file_name,grid_a);
}

template void Write_Face_Field_On_Cells<2>(const std::string& file_name,const MacGrid<2>& mac_grid,const FaceField<real,2>& val);
template void Write_Face_Field_On_Cells<3>(const std::string& file_name,const MacGrid<3>& mac_grid,const FaceField<real,3>& val);

template<int d> void Write_Face_Field_On_Faces(const std::string& file_name,const MacGrid<d>& mac_grid,const FaceField<real,d>& face_field)
{Typedef_VectorDii(d);
	int n=(int)mac_grid.Number_Of_Faces();
	float* xf=new float[n*8];
	memset((void*)xf,0,sizeof(float)*n*8);

	int s=0;for(int axis=0;axis<d;axis++){
		int face_num=mac_grid.Number_Of_Faces(axis);
		if(axis>0)s+=mac_grid.Number_Of_Faces(axis-1);
		#pragma omp parallel for
		for(int f=0;f<face_num;f++){VectorDi face=mac_grid.Face_Coord(axis,f);int i=s+f;
			VectorD pos=mac_grid.Face_Center(axis,face);
			VectorD pos2=pos+VectorD::Unit(axis)*face_field(axis,face);
			if constexpr(d==2){
				xf[i*8+0]=(float)pos[0];
				xf[i*8+1]=(float)pos[1];

				xf[i*8+4]=(float)pos2[0];
				xf[i*8+5]=(float)pos2[1];}
			else if constexpr(d==3){
				xf[i*8+0]=(float)pos[0];
				xf[i*8+1]=(float)pos[1];
				xf[i*8+2]=(float)pos[2];

				xf[i*8+4]=(float)pos2[0];
				xf[i*8+5]=(float)pos2[1];
				xf[i*8+6]=(float)pos2[2];}}}

	std::ofstream output(file_name,std::ios::binary);if(!output)return;
	File::Write_Binary(output,n*8);
	File::Write_Binary_Array(output,xf,n*8);
	delete [] xf;
}

template void Write_Face_Field_On_Faces<2>(const std::string&,const MacGrid<2>&,const FaceField<real,2>&);
template void Write_Face_Field_On_Faces<3>(const std::string&,const MacGrid<3>&,const FaceField<real,3>&);

template<int d,class T> void Write_Points_Float(std::string file_name,const Array<Vector<T,d>>& X)
{
	int n=(int)X.size();
	float* xf=new float[n*4];
	memset((void*)xf,0,sizeof(float)*n*4);
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		if constexpr(d==2){
			xf[i*4]=(float)X[i][0];
			xf[i*4+1]=(float)X[i][1];}
		else if constexpr(d==3){
			xf[i*4]=(float)X[i][0];
			xf[i*4+1]=(float)X[i][1];
			xf[i*4+2]=(float)X[i][2];}}
	std::ofstream output(file_name,std::ios::binary);if(!output)return;
	File::Write_Binary(output,n*4);
	File::Write_Binary_Array(output,xf,n*4);
	delete[] xf;
}

template void Write_Points_Float<2,real>(std::string,const Array<Vector2>&);
template void Write_Points_Float<3,real>(std::string,const Array<Vector3>&);
template void Write_Points_Float<2,float>(std::string,const Array<Vector2f>&);
template void Write_Points_Float<3,float>(std::string,const Array<Vector3f>&);

template<int d,class T> void Write_Vectors_Float(std::string file_name,const Array<Vector<T,d> >& X,const Array<Vector<T,d> >& V,bool as_ends)
{
	int n=(int)X.size();
	float* xf=new float[n*8];
	memset((void*)xf,0,sizeof(float)*n*8);
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		if constexpr(d==2){
			xf[i*8]=(float)X[i][0];
			xf[i*8+1]=(float)X[i][1];
			if(as_ends){
				xf[i*8+4]=(float)(V[i][0]);
				xf[i*8+5]=(float)(V[i][1]);}
			else {
				xf[i*8+4]=(float)(X[i][0]+V[i][0]);
				xf[i*8+5]=(float)(X[i][1]+V[i][1]);}}

		else if constexpr(d==3){
			xf[i*8]=(float)X[i][0];
			xf[i*8+1]=(float)X[i][1];
			xf[i*8+2]=(float)X[i][2];
			if(as_ends){
				xf[i*8+4]=(float)(V[i][0]);
				xf[i*8+5]=(float)(V[i][1]);
				xf[i*8+6]=(float)(V[i][2]);}
			else {
				xf[i*8+4]=(float)(X[i][0]+V[i][0]);
				xf[i*8+5]=(float)(X[i][1]+V[i][1]);
				xf[i*8+6]=(float)(X[i][2]+V[i][2]);}}}

	std::ofstream output(file_name,std::ios::binary);if(!output)return;
	File::Write_Binary(output,n*8);
	File::Write_Binary_Array(output,xf,n*8);
	delete[] xf;
}

template void Write_Vectors_Float<2,double>(std::string,const Array<Vector2>&,const Array<Vector2>&,bool);
template void Write_Vectors_Float<3,double>(std::string,const Array<Vector3>&,const Array<Vector3>&,bool);
template void Write_Vectors_Float<2,float>(std::string,const Array<Vector2f>&,const Array<Vector2f>&,bool);
template void Write_Vectors_Float<3,float>(std::string,const Array<Vector3f>&,const Array<Vector3f>&,bool);

template<int d> void Write_Scalars_As_Points_Float(std::string file_name,const GeometryParticles<d>& points,const Array<real>& arr,const real scale)
{
	int pn=points.Size();
	if(arr.size()!=pn)AuxFunc::Crash_With_Info("RenderFunc::Write_Scalars_As_Points error: size not match");
	Array<Vector<real,d>> to_write(pn);
	#pragma omp parallel for
	for(int i=0;i<pn;i++){
		to_write[i]=points.X(i)+ points.Normal(i).normalized()*arr[i]*scale;}
	RenderFunc::Write_Points_Float<d,real>(file_name,to_write);
}

template void Write_Scalars_As_Points_Float<2>(std::string,const GeometryParticles<2>&,const Array<real>&,const real);
template void Write_Scalars_As_Points_Float<3>(std::string,const GeometryParticles<3>&,const Array<real>&,const real);

template<int d> void Write_Customized_Segments_Float(std::string file_name,const Array<Vector<real,d>>& xs,const Array<Vector<real,d>>& normals,const Array<real>& len_arr,const real scale)
{
	Array<Vector<real,d> > to_write;
	int pn=(int)len_arr.size();
	to_write.resize(pn);
	#pragma omp parallel for
	for(int i=0;i<pn;i++){
		to_write[i]=normals[i].normalized()*len_arr[i]*scale;}
	Write_Vectors_Float<d,real>(file_name,xs,to_write);//in Particles.h
}

template void Write_Customized_Segments_Float<2>(std::string,const Array<Vector2>&,const Array<Vector2>&,const Array<real>&,const real);
template void Write_Customized_Segments_Float<3>(std::string,const Array<Vector3>&,const Array<Vector3>&,const Array<real>&,const real);

template<int d> void Write_Dirichlet_Boundary_Conditions(std::string file_name,const MacGrid<d>& mac_grid,const BoundaryConditionMacGrid<d>& bc)
{
	Particles<d> particles;
	for(auto p:bc.psi_D_values){
		Vector<int,d> cell=mac_grid.grid.Cell_Coord(p.first);
		Vector<real,d> pos=mac_grid.grid.Center(cell);
		int i=particles.Add_Element();particles.X(i)=pos;}
	particles.Write_To_File_3d(file_name);
}

template void Write_Dirichlet_Boundary_Conditions<2>(std::string,const MacGrid<2>&,const BoundaryConditionMacGrid<2>&);
template void Write_Dirichlet_Boundary_Conditions<3>(std::string,const MacGrid<3>&,const BoundaryConditionMacGrid<3>&);

template<int d> void Write_Neumann_Boundary_Conditions(std::string file_name,const MacGrid<d>& mac_grid,const BoundaryConditionMacGrid<d>& bc)
{
	Particles<d> particles;
	for(auto p:bc.psi_N_values){
		int axis=p.first[0];
		Vector<int,d> face=mac_grid.Face_Coord(axis,p.first[1]);
		Vector<real,d> pos=mac_grid.Face_Center(axis,face);
		int i=particles.Add_Element();particles.X(i)=pos;}
	//File::Write_Binary_To_File(file_name,particles);
	particles.Write_To_File_3d(file_name);
}

template void Write_Neumann_Boundary_Conditions<2>(std::string,const MacGrid<2>&,const BoundaryConditionMacGrid<2>&);
template void Write_Neumann_Boundary_Conditions<3>(std::string,const MacGrid<3>&,const BoundaryConditionMacGrid<3>&);

template<int d> void Write_LevelSet_Mesh(std::string file_name,LevelSet<d>& levelset)
{
	MarchingCubes<d> marching_cubes(levelset);
	marching_cubes.Marching();
	(*marching_cubes.mesh).Write_To_File_3d(file_name);
}

template void Write_LevelSet_Mesh<2>(std::string,LevelSet<2>&);
template void Write_LevelSet_Mesh<3>(std::string,LevelSet<3>&);
}


