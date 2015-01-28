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
        __lua::register_functions(L);
        if(argc>1)
        {
            Lua::Config::DoFile(L,argv[1]);
        }
        for(int i=2;i<argc;++i)
        {
            Lua::Config::DoString(L, argv[i]);
        }



        HCell cell(L,0.0);

        cell.ComputeOutsideComposition(0.0);
        const double zm = cell.ComputeRestingZeta(0.0);

        std::cerr << "Resting Zeta=" << zm << std::endl;
        std::cerr << "Em=" << zm * cell.Z2E*1000.0 << " mV" << std::endl;

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