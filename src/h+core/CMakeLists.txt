
ADD_LIBRARY(h+core STATIC
	luavm.hpp
	luavm.cpp	
	hcell.hpp
	hcell.cpp)

ADD_EXECUTABLE(test-hcore test-hcore.cpp sample.lua)
TARGET_LINK_LIBRARIES(test-hcore h+core)
TARGET_LINK_YOCTO(test-hcore mk lua++ lua chem)

