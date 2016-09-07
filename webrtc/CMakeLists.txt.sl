set(webrtc_include
        common.h
        config.h
        transport.h
        typedefs.h
        common_types.h
        )

add_subdirectory(base)
add_subdirectory(common_audio)
add_subdirectory(modules)
add_subdirectory(system_wrappers)

install(FILES  ${webrtc_include}  DESTINATION ${PROJECT_BINARY_DIR}/include/webrtc)
