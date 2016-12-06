PROJECT(shapes)

ADD_LIBRARY( subd SHARED "subd.cpp" "threadsafeMap.h" )
target_link_libraries(subd ${pkgs_LIBRARIES})

INSTALL(TARGETS subd
         RUNTIME DESTINATION plugins
)

ADD_LIBRARY( alembic SHARED "alembic.cpp" )
target_link_libraries(alembic ${pkgs_LIBRARIES})

INSTALL(TARGETS alembic
         RUNTIME DESTINATION plugins
)