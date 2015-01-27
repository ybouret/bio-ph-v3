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
              const double  t0,
              const string &libID,
              const string &eqsID,
              const string &effID,
              const string &iniID) :
L(vm),
lib(L, libID.c_str() ),
eqs(L, eqsID.c_str(),lib),
N(eqs.N),
M(eqs.M),
params(lib,PARAMETERS,NUM_PARAMS),
nvar(params.nvar),
eff(L, effID.c_str() ),
S0(nvar,0)
{
    std::cerr << lib << std::endl;
    std::cerr << eqs << std::endl;
    std::cerr << "#var=" << params.nvar << std::endl;
    {
        boot ini;
        __lua::load(L,ini, iniID, lib);
        eqs.create(S0, ini, t0);
    }
}



