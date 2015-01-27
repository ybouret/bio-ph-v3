#include "cell.hpp"

HCell:: ~HCell() throw()
{
}

HCell:: HCell( lua_State *vm,
              const string &libID,
              const string &eqsID,
              const string &effID) :
L(vm),
lib(L, libID.c_str() ),
eqs(L, eqsID.c_str(),lib),
eff(L, effID.c_str() )
{
    std::cerr << lib << std::endl;
    std::cerr << eqs << std::endl;
}

