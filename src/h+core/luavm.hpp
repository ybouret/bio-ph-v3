#ifndef H_LUAVM_INCLUDED
#define H_LUAVM_INCLUDED 1

#include "yocto/chemical/lua/io.hpp"
#include "yocto/lua/lua-state.hpp"

using namespace yocto;
using namespace math;
using namespace chemical;

//! LuaVM with extended built-in functions
class LuaVM : public Lua::State
{
public:
    explicit LuaVM();
    virtual ~LuaVM() throw();

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(LuaVM);
};

#endif

