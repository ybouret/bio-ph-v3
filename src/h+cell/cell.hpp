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

    explicit HCell(lua_State *vm,
                   const string &libID,
                   const string &eqsID,
                   const string &effID);
    virtual ~HCell() throw();

    lua_State        *L;
    __lua::Library    lib;
    __lua::Equilibria eqs;
    __lua::Effectors  eff;
    parameters        params;

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(HCell);
};


#endif
