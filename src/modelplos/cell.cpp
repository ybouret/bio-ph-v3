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
    ans += 1.5 * rho[ lib["K+"]->indx ] * lambda_K.pace;

    lambda_Na.rate(rho, tmx, in, out, params);
    ans += rho[ lib["Na+"]->indx ] * lambda_Na.pace;

    lambda_Cl.rate(rho, tmx, in, out, params);
    ans -= rho[ lib["Cl-"]->indx ] * lambda_Cl.pace;
    
    return ans;
}

void  Cell:: SetSteadyStatePotential(double Em)
{
    in[iZeta] = Em * E2Z;

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
