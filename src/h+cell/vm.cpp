#include "vm.hpp"


const double __R__       = 8.314472;
double       Temperature = 273.15+25;
const double __Faraday__ = 96485.3415;




double Psi( double U ) throw()
{
	static const double c[4] = { 1.0, 0.5, 1.0/12, 1.0/720 };
	if( Fabs(U) < 1e-3 )
	{
		return c[0] - c[1] * U + c[2] * U*U + c[3] * U*U*U;
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
L(0)
{
   
    std::cerr << "-- Registering Virtual Machine" << std::endl;
    Lua::State &vm = *this;
    L = vm();
    lua_register(L, "SurfaceAndVolume", SurfaceAndVolume );
    std::cerr << "-- Loading Resources from '" << filename << "'" << std::endl;
    Lua::Config::DoFile(L, filename);
}