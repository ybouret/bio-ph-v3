#include "cell.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/kernel/tao.hpp"
#include "yocto/lua/lua-config.hpp"

Cell:: ~Cell() throw() {}

Cell:: Cell( lua_State *vm, const double t ) :
HCell(vm,t),
iNa( lib["Na+"]->indx ),
iK(  lib["K+"]->indx  ),
iCl( lib["Cl-"]->indx ),
iH(  lib["H+"]->indx  ),
leak_K(  eff["lambda_K"]  ),
leak_Na( eff["lambda_Na"] ),
leak_Cl( eff["lambda_Cl"] ),
NHE( eff["NHE"] ),
AE2( eff["AE2"] ),
NaK( eff["NaK"] )
{
}

#include "yocto/math/fcn/zfind.hpp"

double Cell:: SteadyStateZeta()
{
    tmx = 0.0;
    in  = inside;
    ComputeOutsideComposition(tmx);

    
    const double zeta_step = 0.1;
    double zeta = 0.1;
    numeric<double>::function F(this, & Cell:: ComputePassiveZeroFlux);
    do
    {
        const double Fm = F(-zeta);
        const double Fp = F( zeta);
        //std::cerr << "zeta=" << zeta << ", Em=" << Z2E * zeta << ", Fm=" << Fm << ", Fp=" << Fp << std::endl;
        if(Fm*Fp<0)
        {
            break;
        }
        zeta += zeta_step;
        if(zeta>ZETA_MAX)
            throw exception("Cannot find a resting potential!");
    }
    while(true);

    zfind<double> solve(0);
    zeta = solve(F,-zeta,zeta);
    return zeta;

}



double Cell:: ComputePassiveZeroFlux(double zeta)
{

    in[iZeta] = zeta;

    double ans = 0;
    leak_K.rate(rho, tmx, in, out, params);
    ans += 1.5 * rho[ iK ] * leak_K.pace;

    leak_Na.rate(rho, tmx, in, out, params);
    ans += rho[ iNa ] * leak_Na.pace;

    leak_Cl.rate(rho, tmx, in, out, params);
    ans -= rho[ iCl ] * leak_Cl.pace;
    
    return ans;
}

void  Cell:: SetSteadyStatePotential(double Em)
{
    in[iZeta] = Em * E2Z;
    tmx       = 0.0;
    numeric<double>::function F(this, & Cell:: ComputePassiveZeroKFlux);
    double factor = 1.1;
    const double factor_step = 0.1;
    do
    {
        const double Fm = F(1.0/factor);
        const double Fp = F(factor);
        if(Fm*Fp<0)
        {
            break;
        }
        factor += factor_step;
        if(factor>20)
            throw exception("Need too much expansion to get Em=%g", Em );
    }
    while(true);
    zfind<double> solve(0);
    const double pace = solve(F,1.0/factor,factor);
    std::cerr << "increased by " << pace << std::endl;
    leak_K.pace = pace;

    inside[iZeta] = in[iZeta];
}

double Cell:: ComputePassiveZeroKFlux(double pace)
{
    leak_K.pace = pace;
    return ComputePassiveZeroFlux(in[iZeta]);
}


void Cell:: Setup(double Em)
{
    ComputeOutsideComposition(0.0);
    in = inside;
    
    //__________________________________________________________________________
    //
    // Set potential => Lambda etc
    //__________________________________________________________________________
    SetSteadyStatePotential(Em);

    //__________________________________________________________________________
    //
    // compute the steady state K+ leak
    //__________________________________________________________________________
    tao::ld(rho,0);
    leak_K.rate(rho, tmx, in, out, params);
    const double lambda_K = rho[iK]*leak_K.pace;
    std::cerr << "lambda_K=" << lambda_K << std::endl;
    if(lambda_K>=0)
    {
        throw exception("Invalid Potassium Concentrations or Potential!");
    }

    //__________________________________________________________________________
    //
    // compute scaling for NaK
    //__________________________________________________________________________
    const double rho_NaK = -lambda_K;
    std::cerr << "rho_NaK=" << rho_NaK << std::endl;
    tao::ld(rho,0);
    NaK.rate(rho, tmx, in, out, params);
    const double raw_NaK = rho[iK];
    NaK.pace = rho_NaK/raw_NaK;


    //__________________________________________________________________________
    //
    // compute scaling for AE2/chloride
    //__________________________________________________________________________
    tao::ld(rho,0);
    leak_Cl.rate(rho, tmx, in, out, params);
    const double lambda_Cl = rho[iCl]*leak_Cl.pace;
    if(lambda_Cl>=0)
    {
        throw exception("Invalid Chloride Leak: check concentration/potential");
    }
    const double rho_AE2 = -lambda_Cl;
    std::cerr << "rho_AE2=" << rho_AE2 << std::endl;
    tao::ld(rho,0);
    AE2.rate(rho, tmx, in, out, params);
    const double raw_AE2 = rho[iCl];
    AE2.pace = rho_AE2/raw_AE2;

    //__________________________________________________________________________
    //
    // equilibriate for NHE
    //__________________________________________________________________________
    const double rho_NHE = rho_AE2;
    tao::ld(rho,0);
    NHE.rate(rho, tmx, in, out, params);
    const double raw_NHE = rho[iNa];
    NHE.pace = rho_NHE/raw_NHE;

    std::cerr << std::endl;
    eff.rate(rho, tmx, in, out, params);
    std::cerr << "rho=" << std::endl;
    lib.display(std::cerr, rho) << std::endl;;

    eqs.absorb(tmx, rho, in);
    std::cerr << "rho1=" << std::endl;
    lib.display(std::cerr, rho) << std::endl;

}


void Cell:: Rates( array<double> &dYdt, double t, const array<double> &Y )
{
    assert(Y.size()>=nvar);
    assert(dYdt.size()>=nvar);

    //assume Y is normalized ?

    //__________________________________________________________________________
    //
    // compute outside
    //__________________________________________________________________________
    ComputeOutsideComposition(t);

    //__________________________________________________________________________
    //
    //fetch data
    //__________________________________________________________________________
    //const double zeta    = Y[iZeta];
    const double volume  = Y[iVolume];
    const double activeS = Y[iActiveS];
    const double surface = Y[iSurface];

    //__________________________________________________________________________
    //
    //evaluate all surface fluxes
    //__________________________________________________________________________
    eff.rate(rho, t, Y, out, params);



    //__________________________________________________________________________
    //
    //evaluate concentration changes per unit of time
    //__________________________________________________________________________
    const double ratio = activeS/volume;

    for(size_t i=M;i>0;--i)
    {
        dYdt[i] = rho[i] * ratio;
    }

    //__________________________________________________________________________
    //
    // electric equation...
    //__________________________________________________________________________

    //TODO: rewrite with a change of surface
    const double SI_Volume   = volume  * 1e-18;
    //const double SI_surface  = surface * 1e-12;
    const double dQdt        = lib.charge(dYdt) * SI_Volume;
    const double Capa        = surface * Cm;
    const double dzeta       = E2Z * dQdt / Capa;

    dYdt[iZeta] = dzeta;

    //__________________________________________________________________________
    //
    //chemical absorption + constants variation
    //__________________________________________________________________________
    eqs.absorb(t, dYdt, Y);
    

}



