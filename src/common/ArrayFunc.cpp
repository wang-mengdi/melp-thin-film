#include "ArrayFunc.h"

namespace ArrayFunc {
	real Max_Abs_Error(const real x0, const Array<real>& arr) {
		real err = 0;
		for (const real& a : arr) {
			err = std::max(err, abs(x0 - a));
		}
		return err;
	}
	real Max_Rel_Error(const real x0, const Array<real>& arr) {
		AuxFunc::Assert(x0 > 0, "MAx_Rel_Error: must x0>0");
		return abs(Max_Abs_Error(x0, arr) / x0);
	}
	real Upper_Rel_Error(const real x0, const Array<real>& arr)
	{
		real err = 0;
		for (const real& a : arr) {
			err = std::max(err, a - x0);
		}
		return err / x0;
	}
	real L2_Abs_Error(const real x0, const Array<real>& arr)
	{
		real sum = 0;
		for (auto a:arr) {
			sum += (x0 - a) * (x0 - a);
		}
		return sqrt(sum) / arr.size();
	}
	real L2_Rel_Error(const real x0, const Array<real>& arr)
	{
		AuxFunc::Assert(x0 > 0, "L2_Rel_Error: must x0>0");
		return abs(L2_Abs_Error(x0, arr) / x0);
	}
}