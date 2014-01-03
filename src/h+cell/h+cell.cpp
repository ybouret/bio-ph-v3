#include "cell.hpp"
#include "yocto/exception.hpp"
#include "yocto/fs/vfs.hpp"

int main(int argc, char *argv[])
{
    static const char *progname = vfs::get_base_name(argv[0]);
    try
    {
        if(argc<=1)
            throw exception("usage: %s config.lua ...", progname);
        
        Lua::State VM;
        lua_State *L = VM();
        Cell cell(L);
        return 0;
    }
    catch(const exception &e)
    {
        
        std::cerr << "*** in " << progname << std::endl;
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "*** in " << progname << std::endl;
        std::cerr << "unhandled exception!" << std::endl;
    }
    return 1;
}