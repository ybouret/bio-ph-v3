ADD_EXECUTABLE(lactates main.cpp
lactates.lua)

TARGET_LINK_LIBRARIES(lactates h+core)
TARGET_LINK_YOCTO(lactates mk lua++ lua chem)

INSTALL( TARGETS lactates  RUNTIME DESTINATION bin )

