#include "hcell.hpp"
#include "yocto/program.hpp"
#include "yocto/lua/lua-config.hpp"

YOCTO_PROGRAM_START()
{
    LuaVM vm;
    if(argc>1)
    {
        Lua::Config::DoFile(vm(),argv[1]);
    }

    HCell cell(vm(),0.0,NULL,0);
}
YOCTO_PROGRAM_END()

