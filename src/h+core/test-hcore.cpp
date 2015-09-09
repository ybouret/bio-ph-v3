#include "hcell.hpp"
#include "yocto/program.hpp"
#include "yocto/lua/lua-config.hpp"



class Cell : public HCell
{
public:
    Cell(lua_State *vm) : HCell(vm,0,0,0) {}
    virtual ~Cell() throw() {}

    virtual void Rates( array<double> &dYdt, double t, const array<double> &Y )
    {
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

    cell.ComputeOutsideComposition(0);
    cell.ComputeFluxes(0);

    cell.lib.display(std::cerr,cell.rho) << std::endl;
}
YOCTO_PROGRAM_END()

