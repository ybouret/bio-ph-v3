#include "../h+cell/cell.hpp"

#include "yocto/lua/lua-config.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/kernel/tao.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/round.hpp"
#include "yocto/program.hpp"

class Cell : public HCell
{
public:
    explicit Cell( lua_State *vm, const double t) :
    HCell(vm,t)
    {

    }

    virtual ~Cell() throw()
    {
    }

    virtual void   Rates( array<double> &dYdt, double t, const array<double> &Y )
    {

    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
};

YOCTO_PROGRAM_START()
{
    //--------------------------------------------------------------------------
    //
    // setting up Lua
    //
    //--------------------------------------------------------------------------
    LuaVM VM;
    lua_State *L = VM();
    if(argc>1)
    {
        Lua::Config::DoFile(L,argv[1]);
    }
    for(int i=2;i<argc;++i)
    {
        Lua::Config::DoString(L, argv[i]);
    }

    Cell cell(L,0.0);

    cell.ComputeOutsideComposition(0);

}
YOCTO_PROGRAM_END()
