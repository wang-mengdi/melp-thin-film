#ifndef __BoundaryParticles_h__
#define __BoundaryParticles_h__
#include "GeometryParticles.h"
#include "NeighborSearcher.h"


template<int d> class BoundaryParticles: public GeometryParticles<d>
{Typedef_VectorDii(d);Typedef_MatrixD(d);using Base=GeometryParticles<d>;
public:
	using Base::Size; using Base::Resize; using Base::Add_Element; using Base::Add_Elements; using Base::Join; using Base::Copy_Element_From; using Base::Delete_Elements; using Base::Print_Attributes; using Base::Save_Snapshot; using Base::Load_Snapshot;
	using Base::X;using Base::XRef;using Base::V;using Base::F;using Base::M;
	using Base::I;using Base::E;using Base::G;
	using MatrixT=Matrix<real,d-1>;					////tangential vector type
	using VectorT=Vector<real,d-1>;					////tangential vector type
	using VectorTi=Vector<int,d-1>;					////tangential vector int

	std::shared_ptr<NeighborSearcher<d>> nbs_searcher;

	BoundaryParticles() {
		//New_Attributes();

		nbs_searcher = std::make_shared<NeighborKDTree<d> >();
		nbs_searcher->Update_Points(XRef());
		//std::cout << "BoundaryParticles initialized\n";
		//Points<d>::Print_Attributes();
	}

	void Initialize() {
		Update();
	}

	virtual void Update(void) {
		nbs_searcher->Build_Data(XRef());
	}

	void Apply_Rotation(int idx, MatrixD R) {
		X(idx) = R * X(idx);
		V(idx) = R * V(idx);
		F(idx) = R * F(idx);
		E(idx) = R * E(idx);
	}

	void Apply_Translation(int idx, VectorD T) {
		X(idx) = X(idx) + T;
	}

};

#endif
