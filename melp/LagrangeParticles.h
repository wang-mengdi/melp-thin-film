//////////////////////////////////////////////////////////////////////////
// Lagrangian Particles
// Copyright (c) (2018-), Yitong Deng, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __LagrangeParticles_h__
#define __LagrangeParticles_h__
#include "NAParticles.h"



template<int d> class LagrangeParticles: public NAParticles<d>
{
	Typedef_VectorDii(d);
	Typedef_MatrixD(d);
	using Base= NAParticles<d>;
	using VectorT = Vector<real, d - 1>;
	using MatrixT=Matrix<real,d-1>;
public:
	using Base::nbs_searcher;

	using Base::Size;
	using Base::X;
	using Base::XRef;
	using Base::V;
	using Base::SA;
	using Base::Vol;
	using Base::M;
	using Base::Temp;
	using Base::Soap;
	using Base::WM;
public:
	Array<NeighborInfo> e_nbs;
	////Particle attributes
	
	LagrangeParticles() {
		New_Attributes();
		nbs_searcher = std::make_shared<NeighborKDTree<d> >();
		nbs_searcher->Update_Points(XRef());
		//std::cout << "LagrangeParticles initialized\n";
		//Points<d>::Print_Attributes();
	}

	Declare_Attribute(VectorD, Mom, mom);
	Declare_Attribute(MatrixT, AMom, amom); //Affine momentum
	Declare_Attribute(MatrixT, BM, bm); //B in APIC
	Declare_Attribute(MatrixT, DM, dm); //D in APIC
	Declare_Attribute(real, Is_BH, is_bh); //D in APIC

	//Declare_Attribute_Inherent_Func(mom, amom, bm, dm, blackhole);
	Declare_Attribute_Inherent_Func(mom, amom, bm, dm, is_bh);

	void Initialize() {
		Update();
	}

	//set: Vol,M
	void Set_Values(real area, real thickness, real gamma, real rho, real temp = 300.) {
		int n = Size();
#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			V(i) = VectorD::Zero();
			Vol(i) = area * thickness / n;
			SA(i) = area / n;
			Soap(i) = area * gamma / n;;
			M(i) = Vol(i) * rho;
			Mom(i) = VectorD::Zero();
			AMom(i) = MatrixT::Zero();
			BM(i) = MatrixT::Zero();
			DM(i) = MatrixT::Identity();
			Temp(i) = temp;
			Is_BH(i) = 0.;
		}
	}

	virtual void Update(void) {
		nbs_searcher->Build_Data(XRef());
	}

};



#endif
