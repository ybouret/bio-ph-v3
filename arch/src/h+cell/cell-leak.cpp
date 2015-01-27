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
        const double             J     =  ampl * Psi(zz) * (X_out - X*exp(zz));
        cc.concentration = J/volume;
    }
    
    
    
}


#include "yocto/math/fcn/zfind.hpp"

void Cell:: compute_Em()
{
    draw_line();
    std::cerr << "-- Computing Em" << std::endl;
    draw_line();

    numeric<double>::function f( this, & Cell::BiasedPassiveFlux);
    zfind<double> solve(1e-7);
    
    const double zeta = solve(f,-0.1,0.1);
    std::cerr << "zeta=" << zeta << std::endl;
    Em = zeta * ( __R__ * Temperature) / __Faraday__;
    std::cerr << "Em=" << 1000*Em << " mV" << std::endl;
    (void)BiasedPassiveFlux(zeta);
    std::cerr << "lam=" << *sol_tmp << std::endl;
}

double Cell:: BiasedPassiveFlux(double zeta)
{
    chemical::solution &lam = *sol_tmp;
    leak(lam,0.0,zeta,*sol_ins,*sol_out);
    return lam.sum_zC() + 0.5 * lam[ "K+" ];
}


void Cell:: AdjustPermeabilities(double alpha)
{
    chemical::collection::iterator j     = lib.begin();
    chemical::solution::iterator   i     = sol_ins->begin();
    chemical::solution::iterator   i_out = sol_out->begin();
    
    for(size_t k=nsp;k>0;--k,++i,++i_out,++j)
    {
        chemical::species       &sp    = **j;
        const int                 z    = sp.z;
        Permeability            &perm  = sp.data.as<Permeability>();
        const double             X     = (*i).concentration;
        const double             X_out = (*i_out).concentration;
        perm.factor  = 1;
        
        //if( sp.name != "K+" ) continue;
        
        if(z>0)
        {
            if(X>X_out)
            {
                perm.factor = alpha;
                continue;
            }
            
            if(X<X_out)
            {
                perm.factor = 1.0 / alpha;
                continue;
            }
        }
        
        if(z<0)
        {
            
            if(X<X_out)
            {
                perm.factor = alpha;
                continue;
            }
            
            if(X>X_out)
            {
                perm.factor = 1.0 / alpha;
                continue;
            }
            
        }
    }

}


double Cell:: ScaledPassiveFlux(double alpha)
{
    AdjustPermeabilities(alpha);
    const double   zeta = (__Faraday__ * Em) / ( __R__ * Temperature);
    const double   ans  = BiasedPassiveFlux(zeta);
    
    // restore permeabilities
    chemical::collection::iterator j     = lib.begin();
    for(size_t k=1;k<=nsp;++k)
    {
        chemical::species &sp             = **j;
        sp.data.as<Permeability>().factor = 1;
    }
    
    return ans;
}


void Cell:: adjust_Em()
{
    draw_line();
    std::cerr << "Adjusting Permeabilities for Em=" << Em * 1000 << " mV" << std::endl;
    draw_line();
    numeric<double>::function func( this, & Cell:: ScaledPassiveFlux );
    zfind<double> solve( 1e-7 );
    const double  alpha = solve(func,0.5,1.5);
    std::cerr << "alpha=" << alpha << std::endl;
    AdjustPermeabilities(alpha);
    for( chemical::collection::iterator j = lib.begin(); j != lib.end(); ++j )
    {
        const chemical::species &sp = **j;
        std::cerr << "\tPermeability Factor for " << sp.name << " : " << sp.data.as<Permeability>().factor << std::endl;
    }
    
}

