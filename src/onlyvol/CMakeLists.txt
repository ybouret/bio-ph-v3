ADD_EXECUTABLE(onlyvol main.cpp
volume.lua
)

TARGET_LINK_LIBRARIES(onlyvol h+cell)
TARGET_LINK_YOCTO(onlyvol mk lua++ lua chem)

INSTALL( TARGETS onlyvol  RUNTIME DESTINATION bin )

