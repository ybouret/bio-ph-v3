#include "cell.hpp"
#include "yocto/exception.hpp"
#include "yocto/fs/vfs.hpp"
#include "yocto/ios/ocstream.hpp"

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
        double        t  = 0;
        const double  dt = 1;
        ios::ocstream fp("toto.dat",false);
        
        cell.initialize(t);
        cell.leak( *cell.sol_tmp, t, (cell.Em*__Faraday__)/(__R__*Temperature), * cell.sol_ins, * cell.sol_out);
        std::cerr << "lambda1=" << *cell.sol_tmp << std::endl;
        
        
        cell.save_header(fp);
        cell.save_values(t, fp);
        while(t<1000)
        {
            const double t1 = t+dt;
            cell.step(t,t1);
            //std::cerr << "ctrl=" << cell.ctrl << std::endl;
            t=t1;
            cell.save_values(t,fp);
        }
        
        
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