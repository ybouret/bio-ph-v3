#include "cell.hpp"

void Cell:: adjust_effectors()
{
    draw_line();
    std::cerr << "-- Adjusting Effectors Vmax" << std::endl;
    draw_line();
    
    const double        zeta =( __Faraday__ * Em ) / (__R__ * Temperature);
    chemical::solution &lam  = *sol_tmp;
    leak(lam, 0.0, zeta, *sol_ins, *sol_out);
    std::cerr << "lambda=" << lam << std::endl;
    
    chemical::solution rate(lib);
    
    const double lam_K = lam["K+"];
    {
        chemical::effector &NaK = eff["NaK"];
        if(lam_K>0)
            throw exception("bad [K+] concentrations/potential");
        rate.ldz();
        NaK.call(rate,0.0,zeta,*sol_ins,*sol_out);
        std::cerr << "NaK=" << rate << std::endl;
        const double rate_K = rate["K+"];
        if(rate_K<=0)
            throw exception("bad NaK rate!!");
        NaK.factor = -lam_K/(rate_K);
        std::cerr << "\tNaK.factor=" << NaK.factor << std::endl;
    }
}

void Cell:: append( chemical::solution &rates, double t, double zeta, const chemical::solution &S, const chemical::solution &S_out)
{
    eff.append(rates, t, zeta, S, S_out);
}