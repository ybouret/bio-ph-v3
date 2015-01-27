#include "cell.hpp"

HCell:: ~HCell() throw()
{
}

HCell:: HCell( lua_State *vm, const string &libID) :
L(vm),
lib(vm, libID.c_str() )
{
}

