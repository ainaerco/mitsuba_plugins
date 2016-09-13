PROJECT(shapes)

SET( CXX_FILES 
	"subd.cpp" )

SET( H_FILES "threadsafeMap.h" )

ADD_LIBRARY( subd SHARED ${CXX_FILES} ${H_FILES} )
target_link_libraries(subd ${pkgs_LIBRARIES})

INSTALL(TARGETS subd
         RUNTIME DESTINATION plugins
)