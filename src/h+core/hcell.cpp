#include "hcell.hpp"
#include "yocto/exception.hpp"
#include "yocto/lua/lua-config.hpp"

const char * HCell:: PARAMS_NAMES_REG[] = { "zeta", "V", "S" };
const char * HCell:: PARAMS_LOADS_REG[] = { "zeta0", "volume", "surface" };
const size_t HCell:: PARAMS_EXTRA_NUM   = sizeof(HCell::PARAMS_NAMES_REG)/sizeof(HCell::PARAMS_NAMES_REG[0]);

HCell:: ~HCell() throw() {}


#define HCELL_LUA_GET_(FIELD_NAME) Lua::Config::Get<lua_Number>(L,FIELD_NAME)
#define HCELL_LUA_GET(FIELD)  FIELD(  HCELL_LUA_GET_(#FIELD) )

HCell:: HCell(lua_State   *vm,
              const double t0,
              const char  *extra_params_reg[],
              const size_t extra_params_num) :
L(vm),
lib(vm,"lib"),
eqs(vm,"eqs",lib),
N(eqs.N),
M(eqs.M),
//params_reg(),
params(lib, PARAMS_NAMES_REG, PARAMS_LOADS_REG, PARAMS_EXTRA_NUM),
nvar(params.count),
eff(vm,"eff"),
inside(nvar,0.0),
tmx(t0),
in(nvar,0),
outside(),
out(M,0.0),
weights(),

HCELL_LUA_GET(T),
E2Z( Y_FARADAY / (1000.0*Y_R*T) ),
Z2E(1.0/E2Z),
HCELL_LUA_GET(Cm),

odeint( HCELL_LUA_GET_("ftol") ),
HCELL_LUA_GET(diff_h),
diffeq(this, & HCell::Call),
ncalls(0)
{

    //params_reg.release();
    std::cerr << lib << std::endl;
    std::cerr << eqs << std::endl;
    std::cerr << eff << std::endl;

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
    __lua::load(L,inside,params);
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

void HCell:: Call( array<double> &dYdt, double t, const array<double> &Y )
{
    ++ncalls;
}


#if 0
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
#endif


void HCell:: add_header( ios::ostream &fp ) const
{
    fp << " pH";
    fp << " Em";
    fp << " volume";
    fp << " surface";
    fp << " deltaOsm";

    for(library::const_iterator i = lib.begin(); i != lib.end(); ++i)
    {
        const string &id = (*i)->name;
        fp << ' ' << id;
    }
}

void HCell:: add_values( ios::ostream &fp, const array<double> &Y ) const
{
    assert(Y.size()>=nvar);
    fp(" %.15g",lib.pH(Y));
    fp(" %.15g",Y[M+1]*Z2E);
    fp(" %.15g",Y[M+2]);
    fp(" %.15g",Y[M+3]);
    fp(" %.15g",lib.osmolarity(Y)-lib.osmolarity(out));
    for(size_t i=1;i<=M;++i)
    {
        fp(" %.15g", Y[i]);
    }
    
}

