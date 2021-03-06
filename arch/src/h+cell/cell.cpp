#include "cell.hpp"
#include "yocto/exception.hpp"


Permeability:: ~Permeability() throw()
{
}

Permeability:: Permeability( lua_State *L, const string &fn ) :
SP(L,fn,true),
factor(1)
{
}


Cell:: ~Cell() throw()
{
    
}


#define __CELL(FIELD) FIELD( Lua::Config::Get<lua_Number>(L,#FIELD) )

Cell:: Cell(const string &filename) :
VM(filename),
lib(),
nsp(0),
nvar(0),
idxE(0),
eqs(),
__CELL(surface),
__CELL(volume),
__CELL(Cm),        // Farad/microns^2
Capa(surface*Cm),
Em(0),
init_ins(),
dt_ini(1e-4),
ctrl(dt_ini),
odeint(1e-7),
X(),
drvs(this, & Cell::ComputeFields  ),
cb(  this, & Cell::NormalizeState ),
species_ctor_cb(this, & Cell::species_ctor_fn)
{
    const string code = "function SP_ZERO(t,zeta) return 0; end";
    Lua::Config::DoString(L, code);
    
    draw_line();
    std::cerr << "-- Cell Description" << std::endl;
    draw_line();
    std::cerr << "\tsurface = " << surface << " mu^2" << std::endl;
    std::cerr << "\tvolume  = " << volume  << " mu^3" << std::endl;
    std::cerr << "\tcapa    = " << Capa    << std::endl;
    
    //__________________________________________________________________________
    //
    // Loading Species
    //__________________________________________________________________________
    std::cerr << "-- \tReading Species" << std::endl;
    // parse data
    chemical::_lua::load(L, lib, "species", &species_ctor_cb);
    (size_t &)nsp = lib.size();
    if(nsp<=0)
        throw exception("no species");
    (size_t&)nvar = nsp + extra_nvar;
    (size_t&)idxE = nsp + 1;
    std::cerr << lib << std::endl;
    
    //__________________________________________________________________________
    //
    // create solver
    //__________________________________________________________________________
    odeint.start(nvar);
    X.make(nvar,0);
    
    //__________________________________________________________________________
    //
    // create solutions
    //__________________________________________________________________________
    sol_ins.reset( new chemical::solution(lib) );
    sol_out.reset( new chemical::solution(lib) );
    sol_tmp.reset( new chemical::solution(lib) );
    
    //__________________________________________________________________________
    //
    // Loading equilibria
    //__________________________________________________________________________
    std::cerr << "-- \tReading Equilibria" << std::endl;
    chemical::_lua::load(L, lib, eqs, "equilibria");
    eqs.build_from(lib);
    
    //__________________________________________________________________________
    //
    // Loading inside initializer
    //__________________________________________________________________________
    std::cerr << "-- \tReading Inside Initializer" << std::endl;
    chemical::_lua::load(L, init_ins, "inside", lib);
    
    
    //__________________________________________________________________________
    //
    // Loading effectors
    //__________________________________________________________________________
    std::cerr << "-- \tReading Effectors" << std::endl;
    chemical::_lua::load(L, eff, "effectors");
    
    //__________________________________________________________________________
    //
    // Loading outside initializer/solution
    //__________________________________________________________________________
    std::cerr << "-- \tReading Outside Compositions" << std::endl;
    lua_getglobal(L, "outside" );
    if( !lua_istable(L,-1) )
        throw exception("'outside' must be a table");
    
    const unsigned nout = lua_rawlen(L,-1);
    if(nout<=0)
        throw exception("need at least one name in 'outside'");
    
    vector<string> init_names;
    for(unsigned i=1;i<=nout;++i)
    {
        lua_rawgeti(L, -1, i);
        if( !lua_isstring(L,-1) )
            throw exception("outside #%u is not a string",i);
        
        const string s = lua_tostring(L,-1);
        init_names.push_back(s);
        lua_pop(L,1);
    }
    
    lua_settop(L,0);
    std::cerr << "loading " << init_names << std::endl;
    for(size_t i=1;i<=nout;++i)
    {
        std::cerr << "-- Using '" << init_names[i] << "'" << std::endl;
        chemical::boot::loader init_out;
        chemical::_lua::load(L, init_out, init_names[i], lib);
        std::cerr << "-- Computing Composition" << std::endl;
        init_out(eqs,lib,0.0);
        std::cerr <<  "-- Done" << std::endl;
        sol_out->load(eqs.C);
        std::cerr << init_names[i] << "=" << *sol_out << std::endl;
        std::cerr << "\tpH=" << sol_out->pH() << std::endl;
        out_mix.push_back( *sol_out );
    }
    
    weights.make(nout,0);
    
    draw_line();
    std::cerr << "-- Initializing Media" << std::endl;
    draw_line();
    initialize(0.0);
    
}

void Cell:: initialize(double t)
{
    init_ins(eqs,lib,t);
    sol_ins->load( eqs.C );
    
    std::cerr << "S=" << *sol_ins << std::endl;
    std::cerr << "pH=" << sol_ins->pH() << std::endl;

    compute_out(t);
    std::cerr << "S_out=" << *sol_out << std::endl;
    std::cerr << "pH_out=" << sol_out->pH() << std::endl;

    ctrl = dt_ini;
}


void Cell:: compute_out(double t)
{
    lua_getglobal(L, "weights");
    if(!lua_isfunction(L, -1))
        throw exception("invalid weights: expecting a function of time");

    lua_pushnumber(L, t);
    const size_t nout = out_mix.size(); assert(nout>0);
    if( lua_pcall(L, 1, nout, 0) )
        throw exception("weights: %s", lua_tostring(L, -1));
    for(size_t i=nout; i>0; --i )
    {
        if( !lua_isnumber(L, -1) )
            throw exception("invalid weight#%u", unsigned(i) );
        weights[i] = lua_tonumber(L, -1);
        lua_pop(L,1);
    }
    
    //std::cerr << "weights=" << weights << std::endl;
    sol_out->mix(eqs,out_mix, weights, t);
    
}


void Cell:: species_ctor_fn(lua_State *L, chemical::species &sp )
{
    const string SP_Z = "SP_ZERO";
    
    std::cerr << "-- Parsing Surface Permeability for " << sp.name << std::endl;
    
    if( lua_isstring(L,-1))
    {
        const string fn = lua_tostring(L, -1);
        std::cerr << "--\t" << fn << std::endl;
        sp.data.build<Permeability,lua_State*,string>( L, fn );
        SP_Function &SP = sp.data.as<Permeability>().SP;
        std::cerr << "--\t\t" << fn << "(0,0.1)=" << SP(0,0.1) << std::endl;
        
        return;
    }
    
    if( lua_isnil(L, -1) )
    {
        std::cerr << "-- setting to zero" << std::endl;
        sp.data.build<Permeability,lua_State*,string>( L, SP_Z );
        return;
    }
    
    throw exception("Permeability for '%s' is not a string/nil", sp.name.c_str());
    
}


void Cell:: save_header( ios::ostream &fp ) const
{
    fp("#t Em pH "); sol_ins->write_header(fp); fp("\n");
}

void Cell:: save_values( double t, ios::ostream &fp ) const
{
    fp("%e",t);
    fp(" %e", Em);
    fp(" %e", sol_ins->pH() );
    sol_ins->write_values(fp);
    fp("\n");
}
