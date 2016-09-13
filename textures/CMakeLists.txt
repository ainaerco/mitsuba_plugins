PROJECT(textures)

ADD_LIBRARY( oiio_texture SHARED "oiio_texture.cpp" )
target_link_libraries(oiio_texture ${pkgs_LIBRARIES})
INSTALL(TARGETS oiio_texture
         RUNTIME DESTINATION plugins
)

ADD_LIBRARY( adjust SHARED "adjust.cpp" )
target_link_libraries(adjust ${pkgs_LIBRARIES})
INSTALL(TARGETS adjust
         RUNTIME DESTINATION plugins
)

ADD_LIBRARY( checker SHARED "checker.cpp" )
target_link_libraries(checker ${pkgs_LIBRARIES})
INSTALL(TARGETS checker
         RUNTIME DESTINATION plugins
)

ADD_LIBRARY( bitmapudim SHARED "bitmapudim.cpp" )
target_link_libraries(bitmapudim ${pkgs_LIBRARIES})
INSTALL(TARGETS bitmapudim
         RUNTIME DESTINATION plugins
)

ADD_LIBRARY( combine SHARED "combine.cpp" )
target_link_libraries(combine ${pkgs_LIBRARIES})
INSTALL(TARGETS combine
         RUNTIME DESTINATION plugins
)

ADD_LIBRARY( blend SHARED "blend.cpp" )
target_link_libraries(blend ${pkgs_LIBRARIES})
INSTALL(TARGETS blend
         RUNTIME DESTINATION plugins
)

