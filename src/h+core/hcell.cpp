#include "hcell.hpp"
#include "yocto/exception.hpp"
#include "yocto/lua/lua-config.hpp"

const char * HCell:: PARAMS_REG[] = { "zeta", "V", "S" };
const size_t HCell:: PARAMS_NUM   = sizeof(HCell::PARAMS_REG)/sizeof(HCell::PARAMS_REG[0]);

HCell:: ~HCell() throw() {}

HCell:: HCell(lua_State   *vm,
              const double t0,
              const char  *extra_params_reg[],
              const size_t extra_params_num) :
L(vm),
lib(vm,"lib"),
eqs(vm,"eqs",lib),
N(eqs.N),
M(eqs.M),
params_reg(),
params(lib,fill_params_reg(extra_params_reg, extra_params_num) ),
nvar(params.nvar),
eff(vm,"eff"),
inside(nvar,0.0),
tmx(t0),
in(nvar,0),
outside(),
out(M,0.0),
weights(),
T( Lua::Config::Get<lua_Number>(L,"T" ) ),
E2Z( Y_FARADAY / (1000.0*Y_R*T) ),
Z2E(1.0/E2Z)
{

    params_reg.release();
    std::cerr << lib << std::endl;
    std::cerr << eqs << std::endl;
    //std::cerr << eff << std::endl;

    //__________________________________________________________________________
    //
    // Computing Initial Inside Composition
    //__________________________________________________________________________
    std::cerr << "Loading Initial Inside Composition" << std::endl;
    {
        boot ini;
        __lua::load(L,ini, "ini", lib);
        eqs.create(inside, ini, t0);
    }
    std::cerr << "inside=" << std::endl;
    lib.display(std::cerr,inside) << std::endl;

    //__________________________________________________________________________
    //
    // Loading Outside possible solutions
    //__________________________________________________________________________
    std::cerr << "Loading Outside Composition(s)" << std::endl;
    __lua::load(L, outside, "out", lib, eqs, t0);
    const size_t nr = outside.rows;
    if(nr<1)
        throw exception("Need at least one outside solution!");

    for(size_t i=1;i<=outside.rows;++i)
    {
        std::cerr << "outside" << i << "=" << std::endl;
        lib.display(std::cerr,outside[i]) << std::endl;
    }



    //__________________________________________________________________________
    //
    // end initializing
    //__________________________________________________________________________
    in = inside;
    
}

const array<string> &HCell:: fill_params_reg(const char **extra_params_reg,
                                             const size_t extra_params_num)
{
    params_reg.free();
    params_reg.reserve(extra_params_num+PARAMS_NUM);

    for(size_t i=0;i<PARAMS_NUM;++i)
    {
        const string p(PARAMS_REG[i]);
        params_reg.push_back(p);
    }

    for(size_t i=0;i<extra_params_num;++i)
    {
        assert(extra_params_reg);
        assert(extra_params_reg[i]);
        const string p(extra_params_reg[i]);
        params_reg.push_back(p);
    }

    return params_reg;
}