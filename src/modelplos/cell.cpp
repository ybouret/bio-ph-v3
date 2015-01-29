#include "cell.hpp"
#include "yocto/exception.hpp"

Cell:: ~Cell() throw() {}

Cell:: Cell( lua_State *vm, const double t ) :
HCell(vm,t),
iNa( lib["Na+"]->indx ),
iK(  lib["K+"]->indx  ),
iCl( lib["Cl-"]->indx ),
iH(  lib["H+"]->indx  )
{
}

#include "yocto/math/fcn/zfind.hpp"

double Cell:: SteadyStateZeta()
{
    tmx = 0.0;
    in  = inside0;
    ComputeOutsideComposition(tmx);

    //eff["lambda_K" ].pace = 1;
    //eff["lambda_Na"].pace = 1;
    //eff["lambda_Cl"].pace = 1;

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
    effector &lambda_K  = eff["lambda_K"];
    effector &lambda_Na = eff["lambda_Na"];
    effector &lambda_Cl = eff["lambda_Cl"];

    in[iZeta] = zeta;

    double ans = 0;
    lambda_K.rate(rho, tmx, in, out, params);
    ans += 1.5 * rho[ iK ] * lambda_K.pace;

    lambda_Na.rate(rho, tmx, in, out, params);
    ans += rho[ iNa ] * lambda_Na.pace;

    lambda_Cl.rate(rho, tmx, in, out, params);
    ans -= rho[ iCl ] * lambda_Cl.pace;
    
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
    eff["lambda_K"].pace = pace;
}

double Cell:: ComputePassiveZeroKFlux(double pace)
{
    eff["lambda_K"].pace = pace;
    return ComputePassiveZeroFlux(in[iZeta]);
}


void Cell:: Setup(double Em)
{
    // Set potential => Lambda etc
    SetSteadyStatePotential(Em);

    effector &NaK       = eff["NaK"];
    effector &lambda_K  = eff["lambda_K"];

    // compute the real potassium leak
    lambda_K.rate(rho, tmx, in, out, params);
    const double lambda_K_steady = rho[iK]*lambda_K.pace;
    std::cerr << "lambda_K=" << lambda_K_steady << std::endl;
    if(lambda_K_steady>=0)
        throw exception("Invalid Potassium Concentrations or Potential!");

    // scale NaK pace
    const double rho_NaK_steady = -lambda_K_steady/2;
    std::cerr << "rho_NaK=" << rho_NaK_steady << std::endl;
    NaK.pace = 1.0;
    NaK.rate(rho, tmx, in, out, params);
    const double rho_NaK = 0.5*rho[iK];
    NaK.pace = rho_NaK_steady/rho_NaK;
    std::cerr << "\tJmax_NaK=" << NaK.pace << std::endl;

    // scale HNE
    effector &lambda_Na = eff["lambda_Na"];
    lambda_Na.rate(rho, tmx, in, out, params);
    const double lambda_Na_steady = rho[iNa]*lambda_Na.pace;
    std::cerr << "lambda_Na=" << lambda_Na_steady << std::endl;
    const double rho_NHE_steady = 3*rho_NaK_steady - lambda_Na_steady;
    if(rho_NHE_steady<=0)
        throw exception("Invalid Potassium/Sodium or Potential: negative NHE activty");

    std::cerr << "rho_NHE=" << rho_NHE_steady << std::endl;
    
    effector &NHE = eff["NHE"];
    NHE.rate(rho, tmx, in, out, params);
    const double rho_NHE = rho[iNa];
    NHE.pace = rho_NHE_steady/rho_NHE;

    std::cerr << "\tJmax_NHE=" << NHE.pace << std::endl;

    // scale AE
    effector &lambda_Cl = eff["lambda_Cl"];
    lambda_Cl.rate(rho, tmx, in, out, params);
    const double rho_AE2_steady = -rho[iCl]*lambda_Cl.pace;
    std::cerr << "rho_AE2=" << rho_AE2_steady << std::endl;
    if(rho_AE2_steady<=0)
        throw exception("Invalid Chloride/Potential: negatice AE2 activity");

    effector &AE2 = eff["AE2"];
    AE2.rate(rho, tmx, in, out, params);
    const double rho_AE2 = in[iCl];
    AE2.pace = rho_AE2_steady/rho_AE2;
    std::cerr << "\tJmax_AE2=" << AE2.pace << std::endl;

    rho.make(M,0);
    eff.rate(rho, tmx, in, out, params);
    std::cerr << "rho=" << std::endl;
    lib.display(std::cerr, rho) << std::endl;
    std::cerr << "charge=" << lib.charge(rho) << std::endl;

    eqs.absorb(tmx, rho, in);
    std::cerr << "rho1=" << std::endl;
    lib.display(std::cerr, rho) << std::endl;

    std::cerr << "charge=" << lib.charge(rho) << std::endl;

}

