#include "luavm.hpp"

LuaVM:: ~LuaVM() throw()
{
}

LuaVM:: LuaVM() : Lua::State()
{
    lua_State *L = (*this)();
    __lua::register_functions(L);
}
