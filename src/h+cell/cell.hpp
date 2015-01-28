#ifndef H_CELL_INCLUDED
#define H_CELL_INCLUDED 1

#include "yocto/chemical/lua/io.hpp"

using namespace yocto;
using namespace math;
using namespace chemical;

class HCell
{
public:
    static const char  *PARAMETERS[];
    static const size_t NUM_PARAMS;

    //!
    /**
     - Loading species from           "lib"
     - Loading equilibria from        "eqs"
     - Loading effectors from         "eff"
     - Loading inside0   from         "ini"
     - Loading outside solutions from "out"
     */
    explicit HCell(lua_State    *vm,
                   const double  t0);
    virtual ~HCell() throw();

    lua_State        *L;
    __lua::Library    lib;     //!< the library
    __lua::Equilibria eqs;     //!< the global chemical system
    const size_t      N;       //!< #eqs
    const size_t      M;       //!< #species
    parameters        params;  //!< extra parameters
    const size_t      nvar;    //!< #nvar for all
    __lua::Effectors  eff;     //!< effectors
    vector_t          inside0; //!< initial inside concentration (+extra vars)
    matrix_t          outside; //!< possible outside solutions
    vector_t          out;     //!< resulting from mix, using the lua "weights" function

    //! using outside solutions...
    void  ComputeOutsideComposition(const double t);

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(HCell);
};


#endif
