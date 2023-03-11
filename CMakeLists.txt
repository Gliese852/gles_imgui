cmake_minimum_required(VERSION 3.23)

project(gles_imgui)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package (X11 REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules (EGL egl)
pkg_check_modules (GLESV2 glesv2)

add_subdirectory(imgui)
include_directories(imgui)

add_executable(gles_imgui main.cpp)
target_link_libraries(gles_imgui ${EGL_LIBRARIES} ${GLESV2_LIBRARIES} ${X11_X11_LIB} m imgui)
