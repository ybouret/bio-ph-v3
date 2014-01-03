#include "cell.hpp"
#include "yocto/exception.hpp"

Cell:: ~Cell() throw()
{
}


Cell:: Cell( lua_State *L ) :
lib(),
nsp(0),
eqs()
{
    //__________________________________________________________________________
    //
    // Loading Species
    //__________________________________________________________________________
    chemical::_lua::load(L, lib, "species");
    (size_t &)nsp = lib.size();
    if(nsp<=0)
        throw exception("no species");
    
    //__________________________________________________________________________
    //
    // Loading equilibria
    //__________________________________________________________________________
    chemical::_lua::load(L, lib, eqs, "equilibria");
    eqs.build_from(lib);
    
    
}

