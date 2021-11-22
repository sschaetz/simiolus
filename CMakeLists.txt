cmake_minimum_required(VERSION 3.10)

project(simiolus)

find_package(glfw3 3.3 REQUIRED)
find_package(OpenGL REQUIRED)

add_subdirectory(
  "vendor/portaudio"
  "build/vendor/portaudio"
  EXCLUDE_FROM_ALL
)

set(IMGUI_PATH "vendor/imgui")
file(GLOB IMGUI_SOURCES ${IMGUI_PATH}/*.cpp) 
add_library("ImGui" STATIC 
  ${IMGUI_SOURCES} 
  ${IMGUI_PATH}/backends/imgui_impl_glfw.cpp
  ${IMGUI_PATH}/backends/imgui_impl_opengl3.cpp
) 
target_include_directories("ImGui" PUBLIC ${IMGUI_PATH})
target_include_directories("ImGui" PUBLIC ${IMGUI_PATH}/backends)

find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3f IMPORTED_TARGET)
include_directories(PkgConfig::FFTW)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(simi main.cpp)
target_link_libraries(simi portaudio)
target_link_libraries(simi glfw)
target_link_libraries(simi ImGui)
target_link_libraries(simi OpenGL::GL)
target_link_libraries(simi ${CMAKE_DL_LIBS})
target_link_libraries(simi PkgConfig::FFTW)
