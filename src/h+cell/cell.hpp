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
    double                Em;      //!< current potential
    
    chemical::boot::loader init_ins; //!< inside initializer
    
    auto_ptr<chemical::solution> sol_ins; //!< one solution inside
    auto_ptr<chemical::solution> sol_out; //!< one solution outside
    vector<chemical::solution>   out_mix; //!< result from mixing
    vector<double>               weights; //!< storing from lua function 'weights'
    auto_ptr<chemical::solution> sol_tmp; //!< store temporary
    
    explicit Cell(const string &filename);
    virtual ~Cell() throw();
    
    //! compute the new out solution from the mix, according to the weights
    void compute_out( double t );
    
    //! compute inside and outside composition
    void initialize(double t);
    
    //! collect the passive leaks
    void leak( chemical::solution &lambda, double t, double zeta, const chemical::solution &S, const chemical::solution &S_out);
    
    //! collect Em@t=0 from the current permeabilities
    void compute_Em();

    //! adjust permeabilities to match Em
    void adjust_Em();
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
    chemical::species_ctor species_ctor_cb;
    void                   species_ctor_fn(lua_State *L,chemical::species &sp);
    double BiasedPassiveFlux(double zeta);
    double ScaledPassiveFlux(double alpha);
    void   AdjustPermeabilities(double alpha);
    
};


#endif
