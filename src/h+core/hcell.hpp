#ifndef H_CELL_INCLUDED
#define H_CELL_INCLUDED 1

#include "luavm.hpp"
#include "yocto/physics/constants.hpp"
#include "yocto/ios/ostream.hpp"
#include "yocto/memory/pooled.hpp"

//! prototype minimal cell
/**
 - concentrations:      mol/L
 - volume:              mu^3 = 1e-18 m^2 = 1e-15 L
 - surface:             mu^2
 - potential:           mV
 - surface capacitance: muF/cm^2, 1e-14 F/mu^2
 */
class HCell
{
public:
    static const char  *PARAMS_REG[]; //!< built-in parameters: zeta, V, S
    static const size_t PARAMS_NUM;   //!< #PARAMS_REG
    typedef vector<string,memory::pooled::allocator> vector_s;
    //!
    /**
     - Loading species from           "lib"
     - Loading equilibria from        "eqs"
     - Loading effectors from         "eff"
     - Loading inside0   from         "ini"
     - Loading outside solutions from "out"
     - The global initial volume and surface must be available
     */
    explicit HCell(lua_State   *vm,
                   const double t0,
                   const char  *extra_params_reg[],
                   const size_t extra_params_num);

    virtual ~HCell() throw();

    lua_State         *L;           //!< internal virtual machine

    //__________________________________________________________________________
    //
    // Chemical Part
    //__________________________________________________________________________
    __lua::Library     lib;         //!< the library
    __lua::Equilibria  eqs;         //!< the global chemical system
    const size_t      &N;           //!< #eqs
    const size_t      &M;           //!< #species
    vector_s           params_reg;  //!< all the species
    parameters         params;      //!< extra parameters, initially 0
    const size_t      &nvar;        //!< params.size
    __lua::Effectors   eff;         //!< effectors
    vector_t           inside;      //!< initial inside concentration (+extra vars)
    double             tmx;         //!< current time
    vector_t           in;          //!< current inside
    matrix_t           outside;     //!< possible outside solutions
    vector_t           out;         //!< resulting from mix, using the lua "weights" function
    vector_t           weights;     //!< to store weights

    //__________________________________________________________________________
    //
    // Physical Part
    //__________________________________________________________________________
    const double       T;           //!< temperature
    const double       E2Z;         //!< zeta  = E2Z * Em, E2Z=F/(1000*R*T)
    const double       Z2E;         //!< Em/mV = Z2E * zeta
    const double       Cm;          //!< surface capacitance

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(HCell);
    const array<string> & fill_params_reg(const char  *extra_params_reg[],
                                          const size_t extra_params_num);
};


#endif

