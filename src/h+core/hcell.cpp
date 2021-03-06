#include "hcell.hpp"
#include "yocto/exception.hpp"
#include "yocto/lua/lua-config.hpp"


////////////////////////////////////////////////////////////////////////////////
//
// HVariables
//
////////////////////////////////////////////////////////////////////////////////
const char * HVariables:: NAMES_REG[] = { "zeta", "V", "S" };
const char * HVariables:: LOADS_REG[] = { "zeta0", "volume", "surface" };
const size_t HVariables:: EXTRA_NUM   = sizeof(HVariables::NAMES_REG)/sizeof(HVariables::NAMES_REG[0]);
HVariables:: ~HVariables() throw() {}

HVariables:: HVariables() :
names(EXTRA_NUM,as_capacity),
loads(EXTRA_NUM,as_capacity)
{
    for(size_t i=0;i<EXTRA_NUM;++i)
    {
        const string name = NAMES_REG[i];
        const string load = LOADS_REG[i];
        names.push_back(name);
        loads.push_back(load);
    }
}

void HVariables:: ReleaseContent() throw()
{
    loads.release();
    names.release();
}


////////////////////////////////////////////////////////////////////////////////
//
// HCells
//
////////////////////////////////////////////////////////////////////////////////




HCell:: ~HCell() throw() {}


#define HCELL_LUA_GET_(FIELD_NAME) Lua::Config::Get<lua_Number>(L,FIELD_NAME)
#define HCELL_LUA_GET(FIELD)  FIELD(  HCELL_LUA_GET_(#FIELD) )

HCell:: HCell(lua_State   *vm,
              const double t0,
              const char  *extra_params_reg[],
              const size_t extra_params_num) :
HVariables(),
L(vm),
lib(vm,"lib"),
eqs(vm,"eqs",lib),
N(eqs.N),
M(eqs.M),
params(lib,names,loads),
nvar(params.count),
iZeta(params["zeta"]),
iVolume(params["V"]),
iSurface(params["S"]),
eff(vm,"eff"),
inside(nvar,0.0),
tmx(t0),
rho(M,0.0),
in(nvar,0),
outside(),
out(M,0.0),
weights(),

HCELL_LUA_GET(T),
E2Z( Y_FARADAY / (Y_R*T) ),
Z2E(1.0/E2Z),
HCELL_LUA_GET(Cm),

odeint( HCELL_LUA_GET_("ftol") ),
HCELL_LUA_GET(diff_h),
diffeq(this, & HCell::Call),
ncalls(0)
{

    ReleaseContent();

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
    lib.display(std::cerr,inside) << std::endl;
    std::cerr << "inside=" << inside << std::endl;
    
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
    std::cerr << "in =" << in     << std::endl;
    ComputeOutsideComposition(t0);
    std::cerr << "out=" << out    << std::endl;
    odeint.start(nvar);

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

////////////////////////////////////////////////////////////////////////////////
//
//...
//
////////////////////////////////////////////////////////////////////////////////

void HCell:: ComputeOutsideComposition(const double t)
{
    const size_t ns = outside.rows; assert(ns>0);
    weights.make(ns,0.0);
    lua_settop(L,0);
    lua_getglobal(L, "weights");

    //-- push time
    lua_pushnumber(L,t);

    //-- call function
    if( lua_pcall(L, 1, ns, 0))
    {
        const char *err = lua_tostring(L, -1);
        throw exception("weights(%g): %s",t,err);
    }

    //-- get weights
    for(size_t i=ns;i>0;--i)
    {
        if( !lua_isnumber(L, -1))
        {
            throw exception("invalid weight #%u", unsigned(i));
        }
        weights[i] = lua_tonumber(L, -1);
        lua_pop(L,1);
    }
    //std::cerr << "weights=" << weights << std::endl;
    eqs.mix(out, outside, weights, t);
    //std::cerr << "out=" << out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
//
// Output functions
//
////////////////////////////////////////////////////////////////////////////////
ios::ostream & HCell:: add_header( ios::ostream &fp ) const
{
    size_t iCol = 1; // time column
    fp << " pH"; ++iCol;
    fp << " Em"; ++iCol;
    fp << " volume"; ++iCol;
    fp << " surface"; ++iCol;
    fp << " deltaOsm"; ++iCol;

    for(library::const_iterator i = lib.begin(); i != lib.end(); ++i)
    {
        const string &id = (*i)->name;
        fp << ' ' << id;
        ++iCol;
        std::cerr << "\t\t[" << id << "]@column " << iCol << std::endl;
    }
    return fp;
}

ios::ostream & HCell:: add_header_csv( ios::ostream &fp ) const
{
    size_t iCol = 1; // time column
    fp << ",\"pH\""; ++iCol;
    fp << ",\"Em\""; ++iCol;
    fp << ",\"volume\""; ++iCol;
    fp << ",\"surface\""; ++iCol;
    fp << ",\"deltaOsm\""; ++iCol;

    for(library::const_iterator i = lib.begin(); i != lib.end(); ++i)
    {
        const string &id = (*i)->name;
        fp << ',' << '"' << id << '"';
        ++iCol;
        std::cerr << "\t\t[" << id << "]@column " << iCol << std::endl;
    }
    return fp;
}

ios::ostream & HCell:: add_values( ios::ostream &fp, const array<double> &Y ) const
{
    assert(Y.size()>=nvar);
    fp(" %.15g",lib.pH(Y));
    fp(" %.15g",Y[iZeta]*Z2E*1000.0);
    fp(" %.15g",Y[iVolume]);
    fp(" %.15g",Y[iSurface]);
    fp(" %.15g",lib.osmolarity(Y)-lib.osmolarity(out));
    for(size_t i=1;i<=M;++i)
    {
        fp(" %.15g", Y[i]);
    }
    return fp;
}

ios::ostream & HCell:: add_values_csv( ios::ostream &fp, const array<double> &Y ) const
{
    assert(Y.size()>=nvar);
    fp(",%.15g",lib.pH(Y));
    fp(",%.15g",Y[iZeta]*Z2E*1000.0);
    fp(",%.15g",Y[iVolume]);
    fp(",%.15g",Y[iSurface]);
    fp(",%.15g",lib.osmolarity(Y)-lib.osmolarity(out));
    for(size_t i=1;i<=M;++i)
    {
        fp(",%.15g", Y[i]);
    }
    return fp;
}


////////////////////////////////////////////////////////////////////////////////
//
// Computing Resting Zeta functions
//
////////////////////////////////////////////////////////////////////////////////
double HCell:: ComputeVolumicChargeRate(double zeta)
{
    in[iZeta] = zeta;
    // compute fluxes in moles/m^2/s
    eff.rate(rho, tmx, in, out, params);

    // convert in dC/dt
    const double V = in[iVolume]*1e-15; // L
    const double S = in[iSurface]*1e-12; // meters
    const double J2C = S/V;
    for(size_t i=1;i<=M;++i)
    {
        rho[i] *= J2C;
    }

    const double zSurf = in[iSurface];
    const double Cs    = (1e-14*Cm);
    const double Capa  = zSurf * Cs;
    const double dzCdt = V*lib.charge(rho);
    const double dQdt  = Y_FARADAY*dzCdt;
    const double dZdt  = (Y_FARADAY*dQdt)/(Y_R*T)/Capa;
    
    
    return dZdt;
}

const double HCell:: ZETA_MAX = 5.0;

#include "yocto/math/fcn/zfind.hpp"

double HCell:: ComputeRestingZeta(const double t)
{
    tmx = t;
    ComputeOutsideComposition(t);
    for(size_t i=1;i<=nvar;++i)
    {
        in[i] = inside[i];
    }

    const double zeta_step = 0.1;
    double zeta = 0.1;
    numeric<double>::function F(this, & HCell:: ComputeVolumicChargeRate);

    do
    {
        const double Fm = F(-zeta);
        const double Fp = F( zeta);
        if(Fm*Fp<0)
        {
            break;
        }
        zeta += zeta_step;
        if(zeta>ZETA_MAX)
            throw exception("Cannot find a resting potential!");
    }
    while(true);

    zfind<double> solve(0);
    zeta = solve(F,-zeta,zeta);
    inside[iZeta] = zeta;
    std::cerr << "zeta=" << zeta << "-> Em=" << Z2E * zeta << " mV, dZdt=" << ComputeVolumicChargeRate(zeta) << std::endl;
    return zeta;
}



////////////////////////////////////////////////////////////////////////////////
//
// differential engine
//
////////////////////////////////////////////////////////////////////////////////
void HCell:: Step(array<double> &Y, double t0, double t1)
{
    double hh = diff_h;
    odeint(diffeq,Y,t0,t1,hh,&eqs.callback);
}

void HCell:: Call( array<double> &dYdt, double t, const array<double> &Y )
{
    ++ncalls;
    Rates(dYdt, t, Y);

    // compute charges
    const double V   = Y[iVolume]*1e-15; //!< in Litters
    const double Smu = Y[iSurface];
    // TODO: this is a generic behavior ?
    const double Cs    = (1e-14*Cm); //!< in F/mu^2
    const double Capa  = Smu * Cs;
    const double dzCdt = V*lib.charge(dYdt);
    const double dQdt  = Y_FARADAY*dzCdt;
    const double dZdt  = (Y_FARADAY*dQdt)/(Y_R*T)/Capa;

    dYdt[iZeta] = dZdt;
}
