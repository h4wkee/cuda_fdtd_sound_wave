cmake_minimum_required(VERSION 3.16)
project(cuda_fdtd_sound_wave)

set(CMAKE_CXX_STANDARD 17)

if(WIN32)
    add_definitions(-D_WIN64)
endif()

find_package(OpenGL)

include_directories(
        include
        dependencies
        dependencies/glad/include
        dependencies/glfw/include
        dependencies/glm/include
)

add_subdirectory(dependencies/glad)
add_subdirectory(dependencies/GLFW)
add_subdirectory(dependencies/glm)

file(GLOB_RECURSE SRC_FILES src/*.cpp)

add_executable(cuda_fdtd_sound_wave ${SRC_FILES})

add_dependencies(cuda_fdtd_sound_wave glad glfw glm)
target_link_libraries(cuda_fdtd_sound_wave glad glfw glm ${OPENGL_gl_LIBRARY})