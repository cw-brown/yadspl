cmake_minimum_required(VERSION 3.20...4.0)

set(PROJECT_NAME yadspl)
project(${PROJECT_NAME} LANGUAGES CXX)


# Compiled settings
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(MSVC)
    add_compile_options(/W4 /permissive-)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Og -fconcepts-diagnostics-depth=2)
endif()

set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/include/external/imgui)
file(GLOB IMGUI_SRC ${IMGUI_DIR}/*.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp)

set(IMPLOT_DIR ${CMAKE_SOURCE_DIR}/include/external/implot)
file(GLOB IMPLOT_SRC ${IMPLOT_DIR}/*.cpp)

find_package(OpenGL REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GLFW REQUIRED glfw3)

# included libraries and headers
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/yadspl-filter)

# source files
file(GLOB SOURCES CONFIGURE_DEPENDS src/*.cpp)

add_executable(${PROJECT_NAME} ${SOURCES} ${IMGUI_SRC} ${IMPLOT_SRC})

target_include_directories(${PROJECT_NAME} PRIVATE ${IMGUI_DIR} ${IMGUI_DIR}/backends ${IMPLOT_DIR} ${GLFW_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME} ${GLFW_LIBRARIES} OpenGL::GL)