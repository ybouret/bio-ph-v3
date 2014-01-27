#ifndef H_CELL_INCLUDED
#define H_CELL_INCLUDED 1

#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/lua/lua-maths.hpp"

#include "yocto/chemical/lua/io.hpp"
#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "yocto/ptr/auto.hpp"

using namespace yocto;
using namespace math;

typedef Lua::Function<double>        SP_Function;  //!< surface permeability

class Cell
{
public:
    chemical::collection  lib;     //!< library
    const size_t          nsp;     //!< #species
    chemical::equilibria  eqs;     //!< equilibria
    
    chemical::boot::loader init_ins; //!< inside initializer
    chemical::boot::loader init_out; //!< outside initializer
    
    auto_ptr<chemical::solution> sol_ins; //!< one solution inside
    auto_ptr<chemical::solution> sol_out; //!< one solution outside
    vector<chemical::solution>   out_mix; //!< result from mixing
    
    explicit Cell( lua_State *L );
    virtual ~Cell() throw();
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
    chemical::species_ctor species_ctor_cb;
    void                   species_ctor_fn(lua_State *L,chemical::species &sp);
};


#endif
