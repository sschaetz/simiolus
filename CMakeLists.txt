cmake_minimum_required(VERSION 3.10)

project(simiolus)

add_subdirectory(
  "vendor/portaudio"
  "build/vendor/portaudio"
  EXCLUDE_FROM_ALL
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(simi main.cpp)
target_link_libraries(simi portaudio)
