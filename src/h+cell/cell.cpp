#include "cell.hpp"
#include "yocto/exception.hpp"
#include "yocto/lua/lua-config.hpp"

HCell:: ~HCell() throw()
{
}


const char *HCell:: PARAMETERS[] =
{
    "zeta",
    "volume",
    "surface",
    "activeS"
};

const size_t HCell:: NUM_PARAMS = sizeof(PARAMETERS)/sizeof(PARAMETERS[0]);

HCell:: HCell( lua_State   *vm,
              const double  t0) :
L(vm),
lib(L, "lib" ),
eqs(L, "eqs",lib),
N(eqs.N),
M(eqs.M),
params(lib,PARAMETERS,NUM_PARAMS),
iZeta( params["zeta"] ),
iVolume( params["volume"] ),
iSurface(params["surface"]),
iActiveS(params["activeS"]),
nvar(params.nvar),
eff(L, "eff" ),
inside(nvar,0),
outside(),
out(),
weights(),
in(nvar,0),
tmx(0),
rho(M,0),
Temperature(298),
E2Z(Y_FARADAY/(Y_R*Temperature)),
Z2E((Y_R*Temperature)/Y_FARADAY),
odeint( Lua::Config::Get<lua_Number>(L,"ftol")   ),
diff_h( Lua::Config::Get<lua_Number>(L,"diff_h") ),
diffeq( this, & HCell::Call )
{
    std::cerr << lib << std::endl;
    std::cerr << eqs << std::endl;
    std::cerr << "#var=" << params.nvar << std::endl;

    //__________________________________________________________________________
    //
    // Computing Initial Inside Composition
    //__________________________________________________________________________
    std::cerr << "Loading Initial Inside Composition" << std::endl;
    {
        boot ini;
        __lua::load(L,ini, "ini", lib);
        eqs.create(inside, ini, t0);
    }
    std::cerr << "inside=" << std::endl;
    lib.display(std::cerr,inside) << std::endl;
    in = inside;

    //__________________________________________________________________________
    //
    // Loading Outside possible solutions
    //__________________________________________________________________________
    std::cerr << "Loading Outside Composition(s)" << std::endl;
    __lua::load(L, outside, "out", lib, eqs, t0);
    const size_t nr = outside.rows;
    if(nr<1)
        throw exception("Need at least one outside solution!");

    for(size_t i=1;i<=outside.rows;++i)
    {
        std::cerr << "outside" << i << "=" << std::endl;
        lib.display(std::cerr,outside[i]) << std::endl;
    }



    out.make(M,0.0);


    //__________________________________________________________________________
    //
    // Loading other parameters
    //__________________________________________________________________________
    inside[iZeta]    = 0;
    inside[iVolume]  = Lua::Config::Get<lua_Number>(L,"volume");
    inside[iSurface] = Lua::Config::Get<lua_Number>(L,"surface");
    inside[iActiveS] = inside[iSurface];

    std::cerr << "volume =" << inside[iVolume]  << " mu^3" << std::endl;
    std::cerr << "surface=" << inside[iSurface] << " mu^2" << std::endl;

    //__________________________________________________________________________
    //
    // preparing integrator
    //__________________________________________________________________________
    odeint.start(nvar);
}


void HCell:: ComputeOutsideComposition(const double t)
{
    const size_t ns = outside.rows; assert(ns>0);
    weights.make(ns,0.0);
    lua_settop(L,0);
    lua_getglobal(L, "weights");

    //-- push time
    lua_pushnumber(L,t);

    //-- call function
    if( lua_pcall(L, 1, ns, 0))
    {
        const char *err = lua_tostring(L, -1);
        throw exception("weights(%g): %s",t,err);
    }

    //-- get weights
    for(size_t i=ns;i>0;--i)
    {
        if( !lua_isnumber(L, -1))
        {
            throw exception("invalid weight #%u", unsigned(i));
        }
        weights[i] = lua_tonumber(L, -1);
        lua_pop(L,1);
    }
    //std::cerr << "weights=" << weights << std::endl;
    eqs.mix(out, outside, weights, t);
    //std::cerr << "out=" << out << std::endl;
}

void HCell:: Call(array<double> &dYdt, double t, const array<double> &Y)
{
    Rates(dYdt,t,Y);
}

void HCell:: Step(array<double> &Y, double t0, double t1)
{
    double hh = diff_h;
    odeint(diffeq,Y,t0,t1,hh,&eqs.callback);
}


const double HCell:: ZETA_MAX = 5.0;

#include "yocto/math/fcn/zfind.hpp"

double HCell:: ComputeRestingZeta(const double t)
{
    tmx = t;
    const double zeta_step = 0.1;
    double zeta = 0.1;
    numeric<double>::function F(this, & HCell:: ComputeFluxes);
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
    inside[iZeta] = zeta;
    return zeta;
}

double HCell:: ComputeFluxes(double zeta)
{
    in[iZeta] = zeta;
    eff.rate(rho, tmx, in, out, params);
    return lib.charge(rho);
}

void HCell:: add_header( ios::ostream &fp ) const
{
    fp << " pH";
    fp << " Em";
    fp << " volume";
    fp << " surface";
    fp << " activeS";

    for(library::const_iterator i = lib.begin(); i != lib.end(); ++i)
    {
        const string &id = (*i)->name;
        fp << ' ' << id;
    }
}

void HCell:: add_values( ios::ostream &fp, const array<double> &Y ) const
{
    assert(Y.size()>=nvar);
    fp(" %.15g",lib.pH(Y));
    fp(" %.15g",Y[iZeta]*Z2E*1000);
    fp(" %.15g",Y[iVolume]);
    fp(" %.15g",Y[iSurface]);
    fp(" %.15g",Y[iActiveS]);

    for(size_t i=1;i<=M;++i)
    {
        fp(" %.15g", Y[i]);
    }

}

