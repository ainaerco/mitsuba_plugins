PROJECT(wardaniso)
SET( CXX_FILES 
	"wardaniso.cpp" )
ADD_LIBRARY( ${PROJECT_NAME} SHARED ${CXX_FILES} )
target_link_libraries(${PROJECT_NAME} ${pkgs_LIBRARIES})
INSTALL(TARGETS ${PROJECT_NAME}
         RUNTIME DESTINATION plugins
)

PROJECT(zmask)
SET( CXX_FILES 
	"zmask.cpp" )
ADD_LIBRARY( ${PROJECT_NAME} SHARED ${CXX_FILES} )
target_link_libraries(${PROJECT_NAME} ${pkgs_LIBRARIES})
INSTALL(TARGETS ${PROJECT_NAME}
         RUNTIME DESTINATION plugins
)
