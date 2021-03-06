#include "vm.hpp"


const double __R__       = 8.314472;
double       Temperature = 273.15+25;
const double __Faraday__ = 96485.3415;




double Psi( double U ) throw()
{
	static const double c[] = { 1.0, 0.5, 1.0/12, 1.0/720, 1.0/30240};
	if( Fabs(U) < 1e-3 )
	{
        const double U2 = U*U;
        const double U4 = U2*U2;
        const double U6 = U4*U2;
		return c[0] - c[1] * U + c[2] * U2 - c[3] * U4  + c[4] * U6;
	}
	else
	{
		return U/(Exp(U)-1.0);
	}
}


VM:: ~VM() throw() {}


static inline int SurfaceAndVolume( lua_State *L )
{
    lua_Number r[3];
	if( lua_gettop(L) != 3 )
	{
		lua_pushstring(L, "incorrect #arguments: needs a,b,c");
		lua_error(L);
	}
	for( int i=1; i <= 3; ++i )
	{
		if( !lua_isnumber(L, i) )
		{
			lua_pushstring(L, "invalid argument, not a number");
			lua_error(L);
		}
		r[i-1] = lua_tonumber(L, i);
	}
	static const double p = 1.6075;
	const double V   = (4.0*numeric<double>::pi/3.0) * r[0] * r[1] * r[2];
	const double tmp = Pow( r[0]*r[1], p ) + Pow( r[0]*r[2], p ) + Pow( r[1]*r[2], p );
	const double S   = (4.0*numeric<double>::pi) * Pow( tmp/3.0, 1.0/p);
    lua_pushnumber(L,S);
	lua_pushnumber(L,V);
	return 2;
    
}


VM:: VM( const string &filename) :
Lua::State(),
L(0),
line_length(63)
{
    draw_line();
    std::cerr << "-- Registering Virtual Machine" << std::endl;
    draw_line();
    Lua::State &vm = *this;
    L = vm();
    lua_register(L, "SurfaceAndVolume", SurfaceAndVolume );
    draw_line();
    std::cerr << "-- Loading Resources from '" << filename << "'" << std::endl;
    Lua::Config::DoFile(L, filename);
    draw_line();
}

void VM:: draw_line() const
{
    for(size_t i=line_length;i>0;--i) std::cerr << "-";
    std::cerr << std::endl;
}

