cmake_minimum_required(VERSION 3.20...4.0)

set(PROJECT_NAME yadspl)
project(${PROJECT_NAME} LANGUAGES CXX)

# Compiled settings
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(MSVC)
    add_compile_options(/W4 /permissive-)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -O3 -fconcepts-diagnostics-depth=2)
endif()

# included libraries and headers
include_directories(${CMAKE_SOURCE_DIR}/include)

# source files
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})
