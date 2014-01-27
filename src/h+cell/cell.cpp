#include "cell.hpp"
#include "yocto/exception.hpp"

Cell:: ~Cell() throw()
{
    
}


Cell:: Cell( lua_State *L ) :
lib(),
nsp(0),
eqs(),
init_ins(),
species_ctor_cb(this, & Cell::species_ctor_fn)
{
    const string code = "function SP_ZERO(t,zeta) return 0; end";
    
    Lua::Config::DoString(L, code);
    
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
    
    // create solution
    sol_ins.reset( new chemical::solution(lib) );
    
    //__________________________________________________________________________
    //
    // Loading equilibria
    //__________________________________________________________________________
    chemical::_lua::load(L, lib, eqs, "equilibria");
    eqs.build_from(lib);
    
    //__________________________________________________________________________
    //
    // Loading intializers
    //__________________________________________________________________________
    chemical::_lua::load(L, init_ins, "inside", lib);
    
    init_ins(eqs,lib, 0.0 );
    sol_ins->load( eqs.C );
    
    std::cerr << "S=" << *sol_ins << std::endl;
    std::cerr << "pH=" << sol_ins->pH() << std::endl;
    
}

void Cell:: species_ctor_fn(lua_State *L, chemical::species &sp )
{
    const string SP_Z = "SP_ZERO";
    const SP_Function SP_ZERO(L,SP_Z,true);

    std::cerr << "-- Parsing Surface Permeability for " << sp.name << std::endl;
    
    if( lua_isstring(L,-1))
    {
        const string fn = lua_tostring(L, -1);
        std::cerr << "--\t" << fn << std::endl;
        SP_Function SP(L,fn,true);
        std::cerr << "--\t\t" << fn << "(0,0.1)=" << SP(0,0.1) << std::endl;
        
        sp.data.make<SP_Function>( SP );
        return;
    }
    
    if( lua_isnil(L, -1) )
    {
        std::cerr << "-- setting to zero" << std::endl;
        sp.data.make<SP_Function>( SP_ZERO );
        return;
    }
    
    throw exception("Permeability for '%s' is not a string/nil", sp.name.c_str());
    
}

