#include "cell.hpp"

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
        const double zm0 = cell.SteadyStateZeta();
        std::cerr << "Zm0=" << zm0 << std::endl;
        std::cerr << "Em0=" << zm0 * cell.Z2E*1000.0 << " mV" << std::endl;
        cell.SetSteadyStatePotential(-40.0e-3);
        const double zm1 = cell.SteadyStateZeta();
        std::cerr << "Em1=" << zm1 * cell.Z2E*1000.0 << " mV" << std::endl;

        cell.Setup(-60.0e-3);

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