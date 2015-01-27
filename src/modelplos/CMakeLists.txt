ADD_EXECUTABLE(modelplos main.cpp
	model.lua
)

TARGET_LINK_LIBRARIES(modelplos h+cell)
TARGET_LINK_YOCTO(modelplos mk lua++ lua chem)

INSTALL( TARGETS modelplos  RUNTIME DESTINATION bin )

