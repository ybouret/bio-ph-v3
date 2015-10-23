ADD_EXECUTABLE(nbclac main.cpp
	cell.hpp cell.cpp
	nbclac.lua
)

TARGET_LINK_LIBRARIES(nbclac h+core)
TARGET_LINK_YOCTO(nbclac mk lua++ lua chem)

INSTALL( TARGETS nbclac RUNTIME DESTINATION bin )

