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
        
        
        chemical::solution leak(cell.lib);
        
        cell.leak(leak, 0.0, 0.0,  * cell.sol_ins, * cell.sol_out );
        
        std::cerr << "leak=" << leak << std::endl;
        
        cell.compute_Em();
        std::cerr << "Em=" << cell.Em*1000 << std::endl;
        
        cell.Em = -60e-3;
        cell.adjust_Em();
        
        
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