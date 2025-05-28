#include "SPHParticles.h"

template<int d>
void SPHParticles<d>::Update(void)
{
	nbs_searcher->Update_Points(this->XRef());
	nbs_searcher->Record_All_Neighbors(kernel.h);
}
template class SPHParticles<2>;
template class SPHParticles<3>;