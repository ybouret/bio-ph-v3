#include "hcell.hpp"
#include "yocto/program.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/math/kernel/tao.hpp"
#include "yocto/ios/ocstream.hpp"


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
        const double S = Y[iSurface]*1e-6; // m^2
        const double V = Y[iVolume]*1e-15; // L
        const double J2C = S/V;
        for(size_t i=M;i>0;--i)
        {
            rho[i] *= J2C;
        }

        // chemical reactions
        eqs.absorb(t,rho,Y);

        for(size_t i=M;i>0;--i)
        {
            dYdt[i] = rho[i];
        }
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
};

YOCTO_PROGRAM_START()
{
    LuaVM vm;
    if(argc>1)
    {
        Lua::Config::DoFile(vm(),argv[1]);
    }

    Cell cell(vm());
    std::cerr << "Surf_exp=" << Lua::Config::Get<lua_Number>(vm(), "Surf_exp") << std::endl;
    std::cerr << "surface =" << Lua::Config::Get<lua_Number>(vm(), "surface")  << std::endl;

    vector_t Y = cell.inside;

    {
        ios::ocstream fp("sim.dat",false);
        cell.add_header(fp << "#t") << "\n";
        fp("%g ", 0.0); cell.add_values(fp,Y) << "\n";
    }
    

}
YOCTO_PROGRAM_END()

