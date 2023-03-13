project(imgui LANGUAGES CXX)
set(BUILD_SHARED_LIBS OFF)
add_definitions(-DIMGUI_IMPL_OPENGL_ES3)
add_library(${PROJECT_NAME}
	imgui.cpp
	imgui_draw.cpp
	imgui_demo.cpp
	imgui_tables.cpp
	imgui_widgets.cpp
	imgui_impl_opengl3.cpp
)
