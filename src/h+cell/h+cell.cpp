#include "cell.hpp"
#include "yocto/exception.hpp"
#include "yocto/fs/vfs.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/round.hpp"
#include "yocto/code/utils.hpp"



int main(int argc, char *argv[])
{
    static const char *progname = vfs::get_base_name(argv[0]);
    try
    {
        if(false)
        {
            ios::ocstream fp("phi.dat",false);
            for(double x=-0.01;x<=0.01;x+=0.0001)
            {
                fp("%e %e\n", x, Psi(x));
            }
        }
        
        if(argc<=1)
            throw exception("usage: %s config.lua ...", progname);
        
        const string cfgfile = argv[1];
        
        
        
        Cell    cell(cfgfile);
        
        cell.compute_Em();
        std::cerr << "Em=" << cell.Em*1000 << std::endl;
        cell.leak( *cell.sol_tmp, 0.0, (cell.Em*__Faraday__)/(__R__*Temperature), * cell.sol_ins, * cell.sol_out);
        std::cerr << "lambda0=" << *cell.sol_tmp << std::endl;
        
        
        cell.Em = -60e-3;
        cell.adjust_Em();
        
        
        cell.adjust_effectors();
        double        t     = 0;
        const double  t_run   = Lua::Config::Get<lua_Number>(cell.L,"t_run");
        double        dt      = Lua::Config::Get<lua_Number>(cell.L,"dt");
        double        dt_save = Lua::Config::Get<lua_Number>(cell.L,"save_every");
        size_t        every   = 0;
        
        simulation_times(dt, dt_save, every);
        size_t        num_iter = max_of<size_t>(every, ceil(t_run / dt) );
        while( 0 != (num_iter%every) ) ++num_iter;
        
        ios::ocstream fp("toto.dat",false);
        
        cell.initialize(t);
        cell.leak( *cell.sol_tmp, t, (cell.Em*__Faraday__)/(__R__*Temperature), * cell.sol_ins, * cell.sol_out);
        std::cerr << "lambda1=" << *cell.sol_tmp << std::endl;
        
        
        Lua::Function<double> P_CO2(cell.L,"P_CO2");
        {
            ios::ocstream fp("co2.dat",false);
            fp( "0 %g\n", P_CO2(0));
        }
        cell.save_header(fp);
        cell.save_values(t, fp);
        
        size_t num_out=0;
        for(size_t iter=1;iter<=num_iter;++iter)
        {
            t               = (iter-1) * dt;
            const double t1 = iter*dt;
            cell.step(t,t1);
            t=t1;
            if(0 == (iter%every) )
            {
                cell.save_values(t,fp);
                {
                    ios::ocstream fp("co2.dat",true);
                    fp( "%g %g\n", t, P_CO2(t1));
                }
                std::cerr << ".";
                ++num_out;
                if(0==(num_out%32)) std::cerr << std::endl;
                std::cerr.flush();
                
            }
        }
        std::cerr << std::endl;
        
        
        return 0;
    }
    catch(const exception &e)
    {
        
        std::cerr << "*** in " << progname << std::endl;
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "*** in " << progname << std::endl;
        std::cerr << "unhandled exception!" << std::endl;
    }
    return 1;
}