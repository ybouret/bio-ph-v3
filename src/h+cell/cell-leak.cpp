#include "cell.hpp"

void Cell:: leak( chemical::solution &lambda, double t, double zeta, const chemical::solution &S, const chemical::solution &S_out)
{
    
    assert( lambda.has_same_components_than(S) );
    assert( lambda.has_same_components_than(S_out));
    
    chemical::solution::const_iterator i     = S.begin();
    chemical::solution::const_iterator i_out = S_out.begin();
    chemical::solution::iterator       j     = lambda.begin();
    
    assert(nsp == lambda.components );
    for(size_t k=nsp;k>0;--k,++j,++i,++i_out)
    {
        chemical::component     &cc    = *j;
        const chemical::species &sp    = *cc.spec;
        const double             zz    = sp.z * zeta;
        const Permeability      &perm  = sp.data.as<Permeability>();
        const double             ampl  = perm.SP(t,zeta) * perm.factor;
        const double             X     = (*i).concentration;
        const double             X_out = (*i_out).concentration;
        const double             J     = - ampl * Psi(zz) * (X_out - X*exp(zz));
        cc.concentration = J/volume;
    }
    
    
    
}
