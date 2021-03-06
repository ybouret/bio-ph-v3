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
    
    const double lam_Cl = lam["Cl-"];
    {
        chemical::effector &AE = eff["AE"];
        if(lam_Cl>0)
        {
            throw exception("bad [Cl-] concentrations/potential");
        }
        rate.ldz();
        AE.call(rate, 0.0, zeta, *sol_ins, *sol_out);
        std::cerr << "AE=" << rate << std::endl;
        const double rate_Cl = rate["Cl-"];
        if(rate_Cl<=0) throw exception("bad AE rate!");
        AE.factor = -lam_Cl/rate_Cl;
        std::cerr << "\tAE.factor=" << AE.factor << std::endl;
    }
    
    const double lam_Na = lam["Na+"];
    {
        chemical::effector &NHE = eff["NHE"];
        const double lam_NHE = -(1.5*lam_K+lam_Na);
        std::cerr << "lam_NHE=" << lam_NHE << std::endl;
        if( lam_NHE < 0 ) throw exception("bad leak for Na+/K+ with NHE");
        rate.ldz();
        NHE.call(rate,0.0,zeta,*sol_ins,*sol_out);
        std::cerr << "NHE=" << rate << std::endl;
        const double rate_Na = rate["Na+"];
        if( rate_Na<=0) throw exception("bad NHE rate!");
        NHE.factor = lam_NHE/rate_Na;
        std::cerr << "\tNHE.factor=" << NHE.factor << std::endl;
    }
    
}

void Cell:: append( chemical::solution &rates, double t, double zeta, const chemical::solution &S, const chemical::solution &S_out)
{
    eff.append(rates, t, zeta, S, S_out);
}