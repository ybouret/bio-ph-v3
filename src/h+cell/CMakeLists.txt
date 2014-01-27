
ADD_LIBRARY(hcell cell.hpp cell.cpp cell-leak.cpp
	vm.hpp	vm.cpp
	../../sample.lua
)

ADD_EXECUTABLE(h+cell h+cell.cpp)
TARGET_LINK_LIBRARIES(h+cell hcell)
TARGET_LINK_YOCTO(h+cell chemical lua++ lua mk )

INSTALL( TARGETS h+cell RUNTIME DESTINATION bin)

