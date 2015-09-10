#include "hcell.hpp"
#include "yocto/program.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/math/kernel/tao.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/round.hpp"


class Cell : public HCell
{
public:
    Cell(lua_State *vm) : HCell(vm,0.0,0,0) {}
    virtual ~Cell() throw() {}

    //! come with a legal composition
    virtual void Rates( array<double> &dYdt, double t, const array<double> &Y )
    {
        assert(Y.size()==nvar);
        // outside
        ComputeOutsideComposition(t);

        // get fluxes: moles/m^2/s
        eff.rate(rho,t,Y,out,params);

        // convert into dCdt
        const double S = Y[iSurface]*1e-12; // m^2
        const double V = Y[iVolume] *1e-15; // L
        const double J2C = S/V;
        for(size_t i=M;i>0;--i)
        {
            dYdt[i] = rho[i] * J2C;
        }



        // chemical reactions
        eqs.absorb(t,dYdt,Y);


        //std::cerr << "dYdt=" << dYdt << std::endl;
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
};

YOCTO_PROGRAM_START()
{
    LuaVM vm;
    lua_State *L = vm();
    if(argc>1)
    {
        Lua::Config::DoFile(L,argv[1]);
    }

    Cell cell(L);
    std::cerr << "Surf_exp=" << Lua::Config::Get<lua_Number>(L, "Surf_exp") << std::endl;
    std::cerr << "surface =" << Lua::Config::Get<lua_Number>(L, "surface")  << std::endl;


    const double zr = cell.ComputeRestingZeta(0);
    std::cerr << "zr=" << zr << std::endl;
    vector_t Y = cell.inside;
    Y[cell.iZeta] += 0.01;
    {
        ios::ocstream fp("sim.dat",false);
        cell.add_header(fp << "#t") << "\n";
        fp("%g", 0.0); cell.add_values(fp,Y) << "\n";
    }


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
    size_t            iw      = 0;

    for(size_t i=1;i<=niter;++i)
    {
        const double t0 = (i-1)*dt;
        const double t  =  i   *dt;
        std::cerr.flush();
        cell.Step(Y,t0,t);
        if(0==(i%every))
        {
            std::cerr << '[' << wheel[ iw++ % (sizeof(wheel)-1) ] << ']' << "\tt=" << t << "         " << '\r';
            ios::ocstream fp("sim.dat",true);
            fp("%g",t);
            cell.add_values(fp,Y) << "\n";
        }

    }
    std::cerr << std::endl;

    std::cerr << "#calls=" << cell.ncalls << std::endl;


}
YOCTO_PROGRAM_END()

