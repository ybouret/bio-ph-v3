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

    explicit HCell(lua_State    *vm,
                   const double  t0,
                   const string &libID,
                   const string &eqsID,
                   const string &effID,
                   const string &iniID);
    virtual ~HCell() throw();

    lua_State        *L;
    __lua::Library    lib;    //!< the library
    __lua::Equilibria eqs;    //!< the global chemical system
    const size_t      N;      //!< #eqs
    const size_t      M;      //!< #species
    parameters        params; //!< extra parameters
    const size_t      nvar;   //!< #nvar for all
    __lua::Effectors  eff;    //!< effectors
    vector_t          S0;     //!< initial inside concentration
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(HCell);
};


#endif
