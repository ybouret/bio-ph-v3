
#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/lua/lua-maths.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/types.hpp"

using namespace yocto;
using namespace math;

extern double       Temperature;
extern const double __Faraday__;
extern const double __R__;
double Psi(double ) throw();

class VM : public Lua::State
{
public:
    explicit VM( const string &filename);
    virtual ~VM() throw();
    
    lua_State *L;
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(VM);
};
