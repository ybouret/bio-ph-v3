#ifndef H_CELL_INCLUDED
#define H_CELL_INCLUDED 1

#include "luavm.hpp"
#include "yocto/physics/constants.hpp"
#include "yocto/ios/ostream.hpp"
#include "yocto/memory/pooled.hpp"

//! prototype minimal cell
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
    __lua::Library     lib;         //!< the library
    __lua::Equilibria  eqs;         //!< the global chemical system
    const size_t      &N;           //!< #eqs
    const size_t      &M;           //!< #species
    vector_s           params_reg;  //!< all the species
    parameters         params;      //!< extra parameters
    const size_t      &nvar;        //!< params.size
    __lua::Effectors   eff;         //!< effectors
    vector_t           inside;      //!< initial inside concentration (+extra vars)
    vector_t           in;          //!< current inside
    matrix_t           outside;     //!< possible outside solutions
    vector_t           out;         //!< resulting from mix, using the lua "weights" function
    vector_t           weights;     //!< to store weights


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(HCell);
    const array<string> & fill_params_reg(const char  *extra_params_reg[],
                                          const size_t extra_params_num);
};


#endif

