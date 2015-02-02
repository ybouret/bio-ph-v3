#include "cell.hpp"

#include "yocto/lua/lua-config.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/kernel/tao.hpp"
#include "yocto/ios/ocstream.hpp"


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

        vector_t Y(cell.nvar,0);
        tao::set(Y,cell.inside);
        std::cerr << "zeta   =" << Y[cell.iZeta]    << ", Em=" << Y[cell.iZeta]*cell.Z2E*1000.0 << std::endl;
        std::cerr << "volume =" << Y[cell.iVolume]  << std::endl;
        std::cerr << "surface=" << Y[cell.iSurface] << std::endl;
        std::cerr << "activeS=" << Y[cell.iActiveS] << std::endl;

        std::cerr << "Y=" << Y << std::endl;

        std::cerr << "OsmIn ="  << cell.lib.osmolarity(Y)        << std::endl;
        std::cerr << "OsmOut=" << cell.lib.osmolarity(cell.out) << std::endl;

        static const char wheel[] = "|/-\\";
        size_t       count= 0;
        const double dt   = 0.5;

        {
            ios::ocstream fp("output.dat",false);
            fp("#t ");
            cell.add_header(fp);
            fp("\n");
        }

        for(double t=0;t<=1800;t+=dt)
        {
            std::cerr << '[' << wheel[ count++ % sizeof(wheel) ] << ']' << "\tt=" << t << "         " << '\r';
            std::cerr.flush();
            cell.Step(Y,t,t+dt);
            ios::ocstream fp("output.dat",true);
            fp("%g",t+dt);
            cell.add_values(fp,Y);
            fp("\n");

        }
        std::cerr << std::endl;

        std::cerr << "#calls=" << cell.ncalls << std::endl;

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