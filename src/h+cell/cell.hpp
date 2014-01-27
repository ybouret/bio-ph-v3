#ifndef H_CELL_INCLUDED
#define H_CELL_INCLUDED 1

#include "vm.hpp"

#include "yocto/chemical/lua/io.hpp"
#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "yocto/ptr/auto.hpp"

using namespace yocto;
using namespace math;

typedef Lua::Function<double>        SP_Function;  //!< surface permeability

class Permeability
{
public:
    mutable SP_Function SP;     //!<
    double              factor; //!< default is 1
    
    Permeability(lua_State *L, const string &fn);
    ~Permeability() throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Permeability);
};


class Cell : public VM
{
public:
    chemical::collection  lib;     //!< library
    const size_t          nsp;     //!< #species
    chemical::equilibria  eqs;     //!< equilibria
    
    const double          surface; //!< cell surface
    const double          volume;  //!< cell volume
    
    
    chemical::boot::loader init_ins; //!< inside initializer
    
    auto_ptr<chemical::solution> sol_ins; //!< one solution inside
    auto_ptr<chemical::solution> sol_out; //!< one solution outside
    vector<chemical::solution>   out_mix; //!< result from mixing
    vector<double>               weights; //!< storing from lua function 'weights'
    
    
    explicit Cell(const string &filename);
    virtual ~Cell() throw();
    
    void compute_out( double t ); //! mix according to weights
    void initialize(double t);
    
    //! collect the passive leaks
    void leak( chemical::solution &lambda, double t, double zeta, const chemical::solution &S, const chemical::solution &S_out);
    
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
    chemical::species_ctor species_ctor_cb;
    void                   species_ctor_fn(lua_State *L,chemical::species &sp);
};


#endif
