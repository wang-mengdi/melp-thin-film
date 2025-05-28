//////////////////////////////////////////////////////////////////////////
// SPH Kernels
// Copyright (c) (2018-),Xiangxin Kong, Mengdi Wang
// Please see simplex/docs/kernels-math-en.md for documentation.
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#pragma once
#include <iostream>
#include <cmath>
#include "Common.h"
#include "Constants.h"
#include "AuxFunc.h"
#include "Json.h"

enum class KernelType { NONE, POLY6 = 0, SPIKY, CUBIC, QUINTIC, GAUSSIAN };

NLOHMANN_JSON_SERIALIZE_ENUM(KernelType, {
	{KernelType::NONE, "NONE"},
	{KernelType::POLY6, "POLY6"},
	{KernelType::SPIKY, "SPIKY"},
	{KernelType::CUBIC, "CUBIC"},
	{KernelType::QUINTIC,"QUINTIC"},
	{KernelType::GAUSSIAN,"GAUSSIAN"},
	})

class UnitKernel {
public:
	Vector<real, 3> alpha;//0,1,2, -> d=1,2,3
	virtual real Weight(const int d, const real r)const = 0;
	virtual real Grad(const int d, const real r)const = 0;
};

class UnitPOLY6 : public UnitKernel {
public:
	UnitPOLY6() { alpha = Vector3(35.0 / 32.0, 4.0 / pi, 315.0 / (64.0 * pi)); }
	virtual real Weight(const int d, const real r)const { return r < 1 ? alpha[d - 1] * MathFunc::Power3(1 - r * r) : 0; }
	virtual real Grad(const int d, const real r)const { return r < 1 ? alpha[d - 1] * 6 * MathFunc::Power2(r * r - 1) * r : 0; }
};

class UnitSPIKY : public UnitKernel {
public:
	UnitSPIKY() { alpha = Vector3(2.0, 10.0 / pi, 15.0 / pi); }
	virtual real Weight(const int d, const real r)const { return r < 1 ? alpha[d - 1] * MathFunc::Power3(1 - r) : 0; }
	virtual real Grad(const int d, const real r)const { return r < 1 ? alpha[d - 1] * (-3) * MathFunc::Power2(1 - r) : 0; }
};

class UnitCUBIC :public UnitKernel {
public:
	UnitCUBIC() { alpha = Vector3(4.0 / 3.0, 40.0 / (7.0 * pi), 8.0 / pi); }
	virtual real Weight(const int d, const real r)const {
		if (0 <= r && r < 0.5) return alpha[d - 1] * ((real)6 * r * r * r - (real)6 * r * r + 1);
		else if (0.5 <= r && r < 1) return alpha[d - 1] * 2 * MathFunc::Power3(1 - r);
		else return 0;
	}
	virtual real Grad(const int d, const real r)const {
		if (0 <= r && r < 0.5) return alpha[d - 1] * 6.0 * r * (3.0 * r - 2.0);
		else if (0.5 <= r && r < 1) return alpha[d - 1] * (-6.0) * MathFunc::Power2(1.0 - r);
		else return 0;
	}
};

class UnitQUINTIC :public UnitKernel {
public:
	UnitQUINTIC() { alpha = Vector3(1.0 / 40.0, 63.0 / (478.0 * pi), 81.0 / (359.0 * pi)); }
	virtual real Weight(const int d, const real r)const {
		if (0 <= r && r < 1.0 / 3) return alpha[d - 1] * (MathFunc::Power5(3 - 3 * r) - 6 * MathFunc::Power5(2 - 3 * r) + 15 * MathFunc::Power5(1 - 3 * r));
		else if (1.0 / 3 <= r && r < 2.0 / 3) return alpha[d - 1] * (MathFunc::Power5(3 - 3 * r) - 6 * MathFunc::Power5(2 - 3 * r));
		else if (2.0 / 3 < r && r < 1.0) return alpha[d - 1] * (MathFunc::Power5(3 - 3 * r));
		else return 0;
	}
	virtual real Grad(const int d, const real r)const {
		if (0 <= r && r < 1.0 / 3)
			return alpha[d - 1] * (-(real)15 * MathFunc::Power4(3 - 3 * r) + (real)90 * MathFunc::Power4(2 - 3 * r) - (real)225 * MathFunc::Power4(1 - 3 * r));
		else if (1.0 / 3 <= r && r < 2.0 / 3)
			return alpha[d - 1] * (-(real)15 * MathFunc::Power4(3 - 3 * r) + (real)90 * MathFunc::Power4(2 - 3 * r));
		else if (2.0 / 3 < r && r < 1.0)
			return alpha[d - 1] * (-(real)15 * MathFunc::Power4(3 - 3 * r));
		else return 0;
	}
};

class UnitGAUSSIAN :public UnitKernel {
public:
	Vector3 beta;
	real trunc_num;//take sqrt(2)*sigma=1.0/trunc_num
	UnitGAUSSIAN(real _trunc_num = 3) :trunc_num(_trunc_num) {
		alpha = Vector3(trunc_num / sqrt(pi), MathFunc::Power2(trunc_num) / pi, MathFunc::Power3(trunc_num) / pow(pi, 1.5));
		beta = alpha * (-2) * MathFunc::Power2(trunc_num);
	}
	virtual real Weight(const int d, const real r)const { return alpha[d - 1] * exp(-MathFunc::Power2(trunc_num * r)); }
	virtual real Grad(const int d, const real r)const { return beta[d - 1] * r * exp(-MathFunc::Power2(trunc_num * r)); }
};

class KernelSPH {
public:
	static UnitPOLY6 poly6;
	static UnitSPIKY spiky;
	static UnitCUBIC cubic;
	static UnitQUINTIC quintic;
	static UnitGAUSSIAN gaussian;
	std::array<UnitKernel*, 5> kernels;
	
	
	//Truncated at h. Only non-zero in [0,h)
	real h;
	std::array<real, 5> h_pows_inv;//3d in maximum, so we must have h_pows[4]
	KernelType ref_type;
	KernelSPH(const real _h = 1.0, const KernelType _type = KernelType::SPIKY) :h(_h),ref_type(_type) {
		h_pows_inv[0] = 1;
		for (int i = 1; i < 5; i++) { h_pows_inv[i] = h_pows_inv[i - 1] / h; }
		kernels[(int)KernelType::POLY6] = &KernelSPH::poly6;
		kernels[(int)KernelType::SPIKY] = &KernelSPH::spiky;
		kernels[(int)KernelType::CUBIC] = &KernelSPH::cubic;
		kernels[(int)KernelType::QUINTIC] = &KernelSPH::quintic;
		kernels[(int)KernelType::GAUSSIAN] = &KernelSPH::gaussian;
	}
	KernelSPH(const json& j) :KernelSPH(j.at("h").get<const real>(), j.at("type").get<const KernelType>()) {}

	real Weight(int d, real r, real h1, KernelType kernel_type = KernelType::NONE)const;
	real Weight(int d, real r, KernelType kernel_type = KernelType::NONE)const;// const { return Weight(d, r, h, kernel_type); }
	template<int d> real Weight(const Vector<real,d> &r, real h1, const KernelType& kernel_type = KernelType::NONE)const{return Weight(d, r.norm(), h1, kernel_type); }
	template<int d> real Weight(const Vector<real, d>& r, const KernelType& kernel_type = KernelType::NONE)const { return Weight(d, r.norm(), kernel_type); }

	template<int d> real Grad_Norm(const real r, real h1, KernelType kernel_type = KernelType::NONE) const {
		if (kernel_type == KernelType::NONE) kernel_type = ref_type;
		return kernels[(int)kernel_type]->Grad(d, fabs(r / h1)) / MathFunc::Quick_Pow(h1, d + 1);
	}
	template<int d> real Grad_Norm(const real r, KernelType kernel_type = KernelType::NONE) const {
		if (kernel_type == KernelType::NONE) kernel_type = ref_type;
		return kernels[(int)kernel_type]->Grad(d, fabs(r / h)) * h_pows_inv[d + 1];
	}
	template<int d> Vector<real, d> Grad(const Vector<real, d>& r, real h1, KernelType kernel_type = KernelType::NONE)const{
		real r_len = r.norm();
		if (r_len == 0) return r;
		return Grad_Norm<d>(r_len, h1, kernel_type) / r_len * r;
	}
	template<int d> Vector<real, d> Grad(const Vector<real, d>& r, KernelType kernel_type = KernelType::NONE)const {
		real r_len = r.norm();
		if (r_len == 0) return r;
		real grad_coeff = Grad_Norm<d>(r_len, kernel_type) / r_len;
		return  grad_coeff * r;
	}
};
