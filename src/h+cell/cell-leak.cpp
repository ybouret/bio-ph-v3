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


#include "yocto/math/fcn/zfind.hpp"

double Cell:: compute_Em()
{
    numeric<double>::function f( this, & Cell::BiasedPassiveFlux);
    
    triplet<double> zeta = { -0.1, 0, 0.1 };
    triplet<double> flux = {f(zeta.a), 0, f(zeta.c) };
    
    while( flux.a*flux.c > 0 )
    {
        flux.a = f( zeta.a -= 0.01 );
        flux.c = f( zeta.c += 0.01 );
    }
    
    
    std::cerr << "found zeta=" << zeta << " / flux=" << flux << std::endl;
    
    zfind<double> solve(1e-7);
    solve.run(f, zeta, flux);
    //std::cerr << "zeta=" << zeta.b << std::endl;
    const double ans = zeta.b * __R__ * Temperature / __Faraday__;
    std::cerr << "Em=" << 1000*ans << " mV" << std::endl;
    (void)BiasedPassiveFlux(zeta.b);
    std::cerr << "lam=" << *sol_tmp << std::endl;
    return ans;
}

double Cell:: BiasedPassiveFlux(double zeta)
{
    chemical::solution &lam = *sol_tmp;
    leak(lam,0.0,zeta,*sol_ins,*sol_out);
    return lam.sum_zC() + 0.5 * lam[ "K+" ];
}