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

HCell:: HCell( lua_State  *vm,
              const string &libID,
              const string &eqsID,
              const string &effID) :
L(vm),
lib(L, libID.c_str() ),
eqs(L, eqsID.c_str(),lib),
eff(L, effID.c_str() ),
params(lib,PARAMETERS,NUM_PARAMS)
{
    std::cerr << lib << std::endl;
    std::cerr << eqs << std::endl;
    std::cerr << "#var=" << params.nvar << std::endl;
}



