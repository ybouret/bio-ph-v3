#include "cell.hpp"
#include "yocto/exception.hpp"

Cell:: ~Cell() throw()
{
    
}


Cell:: Cell( lua_State *L ) :
lib(),
nsp(0),
eqs(),
species_ctor_cb(this, & Cell::species_ctor_fn)
{
    const string code = "function SP_ZERO(t,zeta) return 0; end";
    const string SP_Z = "SP_ZERO";

    Lua::Config::DoString(L, code);
    const SP_Function SP_ZERO(L,SP_Z,true);

    //__________________________________________________________________________
    //
    // Loading Species
    //__________________________________________________________________________
    
    // parse data
    chemical::_lua::load(L, lib, "species", &species_ctor_cb);
    (size_t &)nsp = lib.size();
    if(nsp<=0)
        throw exception("no species");
    std::cerr << lib << std::endl;
    
    // complete data
    for( chemical::collection::iterator i=lib.begin();i!=lib.end();++i)
    {
        chemical::species &sp = **i;
        if( !sp.data.is_active() )
        {
            std::cerr << "-- Setting Zero Permeability for " << sp.name << std::endl;
            sp.data.build<SP_Function,SP_Function>( SP_ZERO );
        }
    }
    
    
    //__________________________________________________________________________
    //
    // Loading equilibria
    //__________________________________________________________________________
    chemical::_lua::load(L, lib, eqs, "equilibria");
    eqs.build_from(lib);
    
    
}

void Cell:: species_ctor_fn(lua_State *L, chemical::species &sp )
{
    std::cerr << "-- Parsing Surface Permeability for " << sp.name << std::endl;
    
    if( !lua_isstring(L, -1) )
        throw exception("Permeability for '%s' is not a function name", sp.name.c_str());
    
    const string fn = lua_tostring(L, -1);
    std::cerr << "--\t" << fn << std::endl;
    SP_Function SP(L,fn,true);
    std::cerr << "--\t\t" << fn << "(0,0.1)=" << SP(0,0.1) << std::endl;
    
    sp.data.build<SP_Function,SP_Function>( SP );
}