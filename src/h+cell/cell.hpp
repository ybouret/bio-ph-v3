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

typedef auto_ptr<chemical::solution> solution_ptr;
typedef Lua::Function<double>        SP_Function;  //!< surface permeability

class Cell
{
public:
    chemical::collection  lib;     //!< library
    const size_t          nsp;     //!< #species
    chemical::equilibria  eqs;     //!< equilibria
    
    chemical::initializer init_ins; //!< inside initializer
    chemical::initializer init_out; //!< outside initializer
    
    solution_ptr          sol_ins; //!< inside solution
    solution_ptr          sol_out; //!< inside solution
    
    explicit Cell( lua_State *L );
    virtual ~Cell() throw();
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
    chemical::species_ctor species_ctor_cb;
    void                   species_ctor_fn(lua_State *L,chemical::species &sp);
};


#endif
