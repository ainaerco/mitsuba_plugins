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

ADD_LIBRARY( variation SHARED "variation.cpp" )
target_link_libraries(variation ${pkgs_LIBRARIES})
INSTALL(TARGETS variation
         RUNTIME DESTINATION plugins
)

ADD_LIBRARY( colorramp SHARED "colorramp.cpp" )
target_link_libraries(colorramp ${pkgs_LIBRARIES})
INSTALL(TARGETS colorramp
         RUNTIME DESTINATION plugins
)

ADD_LIBRARY( osl SHARED "osl.cpp;simplerend.cpp" )
target_link_libraries(osl ${pkgs_LIBRARIES} ${LLVM_LIBRARY})
INSTALL(TARGETS osl
         RUNTIME DESTINATION plugins
)