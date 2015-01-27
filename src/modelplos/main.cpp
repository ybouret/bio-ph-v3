#include "../h+cell/cell.hpp"

#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"

#include "yocto/exception.hpp"

int main(int argc, char *argv[])
{
    try
    {
        //----------------------------------------------------------------------
        //
        // setting up Lua
        //
        //----------------------------------------------------------------------

        Lua::State VM;
        lua_State *L = VM();
        if(argc>1)
        {
            Lua::Config::DoFile(L,argv[1]);
        }
        for(int i=2;i<argc;++i)
        {
            Lua::Config::DoString(L, argv[i]);
        }

        string libID = "lib";
        string eqsID = "eqs";
        string effID = "eff";
        string iniID = "ini";

        HCell cell(L,0.0,libID,eqsID,effID,iniID);

        return 0;
    }
    catch( const exception &e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "Unhandled exception..." << std::endl;
    }
    return 1;
}