#include "cell.hpp"

#include "yocto/lua/lua-config.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/kernel/tao.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/round.hpp"

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

        //----------------------------------------------------------------------
        //
        // setting up the cell
        //
        //----------------------------------------------------------------------
        Cell cell(L,0.0);
        const double zm0 = cell.SteadyStateZeta();
        std::cerr << "Zm0=" << zm0 << std::endl;
        std::cerr << "Em0=" << zm0 * cell.Z2E*1000.0 << " mV" << std::endl;
        cell.SetSteadyStatePotential(-40.0e-3);
        const double zm1 = cell.SteadyStateZeta();
        std::cerr << "Em1=" << zm1 * cell.Z2E*1000.0 << " mV" << std::endl;

        cell.Setup(-60.0e-3);
        cell.eff.rescale_pace(50);

        //----------------------------------------------------------------------
        //
        // phase space
        //
        //----------------------------------------------------------------------
        vector_t Y(cell.nvar,0);
        tao::set(Y,cell.inside);
        std::cerr << "zeta   =" << Y[cell.iZeta]    << ", Em=" << Y[cell.iZeta]*cell.Z2E*1000.0 << std::endl;
        std::cerr << "volume =" << Y[cell.iVolume]  << std::endl;
        std::cerr << "surface=" << Y[cell.iSurface] << std::endl;
        std::cerr << "activeS=" << Y[cell.iActiveS] << std::endl;

        std::cerr << "Y=" << Y << std::endl;

        std::cerr << "OsmIn ="  << cell.lib.osmolarity(Y)        << std::endl;
        std::cerr << "OsmOut=" << cell.lib.osmolarity(cell.out) << std::endl;

        //----------------------------------------------------------------------
        //
        // timings
        //
        //----------------------------------------------------------------------
        double       dt      = Lua::Config::Get<lua_Number>(L,"dt");
        const double t_run   = Lua::Config::Get<lua_Number>(L,"t_run");
        double       dt_save = Lua::Config::Get<lua_Number>(L,"dt_save");
        const size_t every   = simulation_save_every(dt, dt_save);
        const size_t niter   = simulation_iter(t_run, dt, every);

        std::cerr << "dt=" << dt << std::endl;
        std::cerr << "savering every " << every << std::endl;
        std::cerr << "t_run=" << t_run << std::endl;
        std::cerr << "niter=" << niter << std::endl;

        static const char wheel[] = "|/-\\";

        {
            ios::ocstream fp("output.dat",false);
            fp("#t ");
            cell.add_header(fp);
            fp("\n");
            fp("0 ");
            cell.add_values(fp, Y);
            fp("\n");
        }

        for(size_t i=1;i<=niter;++i)
        {
            const double t0 = (i-1)*dt;
            const double t  =  i   *dt;
            std::cerr.flush();
            cell.Step(Y,t0,t);
            if(0==(i%every))
            {
                std::cerr << '[' << wheel[ i % (sizeof(wheel)-1) ] << ']' << "\tt=" << t << "         " << '\r';
                ios::ocstream fp("output.dat",true);
                fp("%g",t);
                cell.add_values(fp,Y);
                fp("\n");
            }

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