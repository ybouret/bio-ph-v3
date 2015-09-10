ADD_EXECUTABLE(modelplos main.cpp
	cell.hpp cell.cpp
	model.lua
	model-old.lua
)

TARGET_LINK_LIBRARIES(modelplos h+core)
TARGET_LINK_YOCTO(modelplos mk lua++ lua chem)

INSTALL( TARGETS modelplos  RUNTIME DESTINATION bin )

