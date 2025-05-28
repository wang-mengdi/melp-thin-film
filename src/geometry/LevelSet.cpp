//////////////////////////////////////////////////////////////////////////
// Level set
// Copyright (c) (2018-), Bo Zhu, Xingyu Ni
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <numeric>
#include <set>
#include <queue>
#include <utility>
#include <iostream>
#include "LevelSet.h"
#include "Constants.h"
#include "MacGrid.h"
#include "Timer.h"

template<int d> LevelSet<d>::LevelSet(const Grid<d>& _grid):intp(nullptr)
{Initialize(_grid);}

template<int d> void LevelSet<d>::Initialize(const Grid<d>& _grid)
{
	grid=_grid;phi.Resize(grid.cell_counts);
	phi.Fill(std::numeric_limits<real>::max());
	intp.reset(new Interpolation<d>(grid));
}

template<int d> void LevelSet<d>::Set_By_Geom(ImplicitGeometry<d>& geom)
{
	grid.Exec_Each(
		[&](const VectorDi& cell) {
			phi(cell) = geom.Phi(grid.Center(cell));
		}
	);
}

template<int d>
void LevelSet<d>::Set_By_Shape(const ImplicitShape<d>& shape)
{
	grid.Exec_Each(
		[&](const VectorDi& cell) {
			phi(cell) = shape.Phi(grid.Center(cell));
		}
	);
}

template<int d> real LevelSet<d>::Phi(const VectorD& pos) const 
{return intp->Interpolate_Centers(phi,pos);}

template<int d> Vector<real,d> LevelSet<d>::Normal(const VectorD& pos) const
{
	VectorD normal;
	for(int i=0;i<d;i++)normal[i]=(Phi(pos+VectorD::Unit(i)*grid.dx)-Phi(pos-VectorD::Unit(i)*grid.dx))/((real)2*grid.dx);
	return normal.normalized();
}

template<int d> Vector<real,d> LevelSet<d>::Gradient(const VectorD& pos) const
{
	VectorD normal;
	for(int i=0;i<d;i++)normal[i]=(Phi(pos+VectorD::Unit(i)*grid.dx)-Phi(pos-VectorD::Unit(i)*grid.dx))/((real)2*grid.dx);
	return normal;
}

template<int d> void LevelSet<d>::Update_Normals(FaceField<real,d>& normals) const
{
	MacGrid<d> mac_grid(grid);
	for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
		#pragma omp parallel for
		for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
			const VectorD& pos=mac_grid.face_grids[axis].Node(face);
			normals(axis,face)=Normal(pos)[axis];}}
}

template<int d> void LevelSet<d>::Update_Normals(Field<VectorD,d>& normals) const
{
	int cell_num=grid.Number_Of_Cells();
	#pragma omp parallel for
	for(int i=0;i<cell_num;i++){const VectorDi& cell=grid.Cell_Coord(i);
		const VectorD& pos=grid.Center(cell);normals(cell)=Normal(pos);}
}

template<int d> real LevelSet<d>::Curvature(const VectorD& pos) const
{
	real one_over_dx=(real)1/grid.dx;real one_over_two_dx=(real).5*one_over_dx;real curvature=(real)0;
	for(int i=0;i<d;i++){
		VectorD normal_left=Normal(pos-VectorD::Unit(i)*grid.dx);
		VectorD normal_right=Normal(pos+VectorD::Unit(i)*grid.dx);
		curvature+=(normal_right[i]-normal_left[i])*one_over_two_dx;}
	return abs(curvature)<one_over_dx?curvature:(curvature<=(real)0?(real)-1:(real)1)*one_over_dx;

}

template<int d> Vector<real,d> LevelSet<d>::Closest_Point(const VectorD& pos,const real epsilon) const
{
	VectorD normal=Gradient(pos);normal.normalize();
	return pos-normal*(Phi(pos)+epsilon);
}

template<int d> Vector<real,d> LevelSet<d>::Closest_Point_With_Iterations(const VectorD& pos,const int max_iter/*=5*/) const
{
	VectorD intf_pos=pos;
	for(int i=0;i<max_iter;i++){
		intf_pos=Closest_Point(intf_pos);
		if(Phi(intf_pos)<(real)0)return intf_pos;}
	return intf_pos;
}

template<int d> real LevelSet<d>::Cell_Fraction(const VectorDi& cell) const
{
	return (real).5-AuxFunc::Clamp(phi(cell),-(real).5*grid.dx,(real).5*grid.dx)/grid.dx;
}

template<int d> real LevelSet<d>::Total_Volume() const
{
	real total_vol=(real)0;
	real cell_vol=std::pow(grid.dx,d);
	iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();
		total_vol+=Cell_Fraction(cell);}
	return total_vol*cell_vol;
}

//////////////////////////////////////////////////////////////////////////
////Fast marching method

template<int d> void LevelSet<d>::Fix_With_Complementary(const LevelSet<d>& negative_levelset)
{
	const int cell_num = grid.Number_Of_Cells();
	grid.Exec_Each(
		[&](const VectorDi& cell) {
			if (phi(cell) > 0) return;
			for (int j = 0;j < Grid<d>::Number_Of_Nb_C();j++) {
				VectorDi nb = grid.Nb_C(cell, j);
				if (!grid.Valid_Cell(nb)) continue;
				if (negative_levelset.Interface(negative_levelset.phi(cell), negative_levelset.phi(nb))) {
					phi(cell) = -negative_levelset.phi(cell);
					return;}}});
}

template<int d> void LevelSet<d>::Fast_Marching(const real band_width)
{
	Timer timer;
	timer.Reset();
	
	Field<real,d> tent(grid.cell_counts,band_width<0?std::numeric_limits<real>::max():band_width);
	Array<ushort> done(grid.cell_counts.prod(),0);
	using PRI=std::pair<real,int>;
	std::priority_queue<PRI,Array<PRI>,std::greater<PRI> > heaps[2];
	const int cell_num=grid.cell_counts.prod();
	//real far_from_intf_phi_val=grid.dx*(real)5;

	//////////////////////////////////////////////////////////////////////////
	////precondition
	//// find interface cells
	#pragma omp parallel for
	for(int i=0;i<cell_num;i++){
		const VectorDi cell=grid.Cell_Coord(i);
		//if(abs(phi(cell))>far_from_intf_phi_val)continue;		////ATTENTION: this might cause problem if the levelset is badly initialized

		for(int j=0;j<Grid<d>::Number_Of_Nb_C();j++){
			VectorDi nb=grid.Nb_C(cell,j);
			if(!grid.Valid_Cell(nb))continue;
			if(Interface(phi(cell),phi(nb))){done[i]=1;break;}}}
	if(verbose)timer.Elapse_And_Output_And_Reset("FMM Precond: find interface");

	//// calculate interface phi values
	#pragma omp parallel for
	for(int c=0;c<cell_num;c++){
		if(!done[c])continue;		////select interface cells
		const VectorDi cell=grid.Cell_Coord(c);

		VectorD correct_phi = VectorD::Ones() * std::numeric_limits<real>::max();
		VectorDi correct_axis = VectorDi::Zero();
		for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
			VectorDi nb = grid.Nb_C(cell, i);
			if (!grid.Valid_Cell(nb)) continue;
			const int nb_idx = grid.Cell_Index(nb);
			if (done[nb_idx] && Interface(phi(cell), phi(nb))) {
				real c_phi = Theta(phi(cell), phi(nb)) * grid.dx; // always non-negative
				int axis = grid.Nb_C_Axis(i);
				correct_axis[axis] = 1;
				correct_phi[axis] = std::min(correct_phi[axis], c_phi);}}
		if (correct_axis != VectorDi::Zero()) {
			real hmnc_mean = (real)0;
			for (int i = 0; i < d; i++) {
				if (correct_axis[i] == 0)continue;
				hmnc_mean += (real)1 / (correct_phi[i] * correct_phi[i]);}
			hmnc_mean = sqrt((real)1 / hmnc_mean);
			tent(cell) = hmnc_mean;}
		else {
			std::cerr << "Error: [Levelset] bad preconditioning" << std::endl;
			std::exit(1);}}
	if(verbose)timer.Elapse_And_Output_And_Reset("FMM Precond: calculate interface phi");

	//// initialize heap with front cells
	#pragma omp parallel for
	for(int i=0;i<cell_num;i++){
		const VectorDi cell=grid.Cell_Coord(i);
		bool is_front=false;
		for(int j=0;j<Grid<d>::Number_Of_Nb_C();j++){
			VectorDi nb=grid.Nb_C(cell,j);
			if(!grid.Valid_Cell(nb))continue;
			const int nb_idx=grid.Cell_Index(nb);
			if(done[nb_idx]){is_front=true;break;}}

		if(is_front){
			real temp=Solve_Eikonal(cell,tent,done);
			if(temp<tent(cell)){
				tent(cell)=temp;
				#pragma omp critical
				{heaps[Sign(phi(cell))>0?0:1].push(PRI(temp,i));}}}}
	if(verbose)timer.Elapse_And_Output_And_Reset("FMM: Build heap");

	//// heap traversing
	#pragma omp parallel for
	for(int h=0;h<2;h++){
		auto& heap=heaps[h];
		while(!heap.empty()){
			const real top_val=heap.top().first;
			const int cell_idx=heap.top().second;
			const VectorDi cell=grid.Cell_Coord(cell_idx);
			heap.pop();
			if(tent(cell)!=top_val)continue;
			done[cell_idx]=1;

			for(int i=0;i<Grid<d>::Number_Of_Nb_C();i++){
				VectorDi nb=grid.Nb_C(cell,i);
				if (!grid.Valid_Cell(nb))continue;
				const int nb_idx=grid.Cell_Index(nb);
				if(!done[nb_idx]){
					real temp=Solve_Eikonal(nb,tent,done);
					#pragma omp critical
					{if(temp<tent(nb)){tent(nb)=temp;heap.push(PRI(temp,nb_idx));}}}}}}

	if(verbose)timer.Elapse_And_Output_And_Reset("FMM: Traverse heap");

	#pragma omp parallel for
	for(int i=0;i<cell_num;i++){
		phi.array[i]=Sign(phi.array[i])*tent.array[i];}
	if(verbose)timer.Elapse_And_Output_And_Reset("FMM: Write to phi");
}

template<int d> bool LevelSet<d>::Solve_Quadratic(const real p1, const real p2, const real dx, real &rst)
{
	if (abs(p1) >= abs(p2) + dx) { rst = p2 + dx; return true; }
	else if (abs(p2) >= abs(p1) + dx) { rst = p1 + dx; return true; }
	else {
		real delta = (real)2 * dx * dx - pow(p1 - p2, 2);
		if (delta < (real)0) { std::cerr << "Error: [Levelset] negative delta in Solve_Quadratic_2" << std::endl; return false; }
		rst = (real).5 * (p1 + p2 + sqrt(delta)); return true;}
}

template<int d> bool LevelSet<d>::Solve_Quadratic(const real p1, const real p2, const real p3, const real dx, real &rst)
{
	real delta = pow(p1 + p2 + p3, 2) - (real)3 * (p1 * p1 + p2 * p2 + p3 * p3 - dx * dx);
	if (delta < (real)0) {
		int i = 0; real p_max = abs(p1); if (abs(p2) > p_max) { i = 1; p_max = abs(p2); }if (abs(p3) > p_max) { i = 2; p_max = abs(p3); }
		real q1, q2; if (i == 0) { q1 = p2; q2 = p3; }
		else if (i == 1) { q1 = p1; q2 = p3; }
		else { q1 = p1; q2 = p2; }
		return Solve_Quadratic(q1, q2, dx, rst);}
	rst = (real)one_third * (p1 + p2 + p3 + sqrt(delta)); return true;
}

template<int d> real LevelSet<d>::Solve_Eikonal(const VectorDi &cell,const Field<real,d>& tent,const Array<ushort>& done)
{
	// calculate correct phi from nb interface cells
	VectorD correct_phi = VectorD::Ones() * std::numeric_limits<real>::max();
	VectorDi correct_axis = VectorDi::Zero();
	for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
		VectorDi nb = grid.Nb_C(cell, i);
		if (!grid.Valid_Cell(nb)) continue;
		const int nb_idx = grid.Cell_Index(nb);
		if(done[nb_idx]){
			int axis = grid.Nb_C_Axis(i); correct_axis[axis] = 1;
			correct_phi[axis] = std::min(correct_phi[axis], tent(nb));}}
	// update phi on the cell
	real new_phi;
	int n = correct_axis.sum();
	switch (n) {
	case 1: {
		real c_phi;
		for (int i = 0; i < d; i++)
			if (correct_axis[i] != 0) { c_phi = correct_phi[i]; break; }
		new_phi = grid.dx + c_phi;
	} break;
	case 2: {
		real p[2];
		int j = 0;
		for (int i = 0; i < d; i++)
			if (correct_axis[i] != 0) p[j++] = correct_phi[i];
		Solve_Quadratic(p[0], p[1], grid.dx, new_phi);
	} break;
	case 3: {
		Solve_Quadratic(correct_phi[0], correct_phi[1], correct_phi[2], grid.dx, new_phi);
	} break;
	default: {
		std::cerr << "Error: [Levelset] bad solving Eikonal" << std::endl;
		std::exit(1);
	} break;}
	return new_phi;
}

template class LevelSet<2>;
template class LevelSet<3>;
