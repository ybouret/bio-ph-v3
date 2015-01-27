#include "cell.hpp"

HCell:: ~HCell() throw()
{
}


const char *HCell:: PARAMETERS[] =
{
    "zeta",
    "volume"
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
S0(nvar,0)
{
    std::cerr << lib << std::endl;
    std::cerr << eqs << std::endl;
    std::cerr << "#var=" << params.nvar << std::endl;
    std::cerr << "Loading Initial Inside Composition" << std::endl;
    {
        boot ini;
        __lua::load(L,ini, "ini", lib);
        eqs.create(S0, ini, t0);
    }
    std::cerr << "S0=" << S0 << std::endl;

    std::cerr << "Loading Outside Composition(s)" << std::endl;
    __lua::load(L, out, "out", lib, eqs, t0);
    std::cerr << "Out=" << out << std::endl;
}



