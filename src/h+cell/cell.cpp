#include "cell.hpp"
#include "yocto/exception.hpp"

HCell:: ~HCell() throw()
{
}


const char *HCell:: PARAMETERS[] =
{
    "zeta",
    "volume",
    "surface"
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
nvar(params.nvar),
eff(L, "eff" ),
inside0(nvar,0),
outside(),
out()
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
        eqs.create(inside0, ini, t0);
    }
    std::cerr << "inside0=" << std::endl;
    lib.display(std::cerr,inside0) << std::endl;

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
}


void HCell:: ComputeOutsideComposition(const double t)
{
    const size_t ns = outside.rows; assert(ns>0);
    vector_t     weights(ns,0.0);
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

    // get weights
    for(size_t i=ns;i>0;--i)
    {
        if( !lua_isnumber(L, -1))
        {
            throw exception("invalid weight #%u", unsigned(i));
        }
        weights[i] = lua_tonumber(L, -1);
        lua_pop(L,1);
    }
    std::cerr << "weights=" << weights << std::endl;
    eqs.mix(out, outside, weights, t);
    std::cerr << "out=" << out << std::endl;
}


