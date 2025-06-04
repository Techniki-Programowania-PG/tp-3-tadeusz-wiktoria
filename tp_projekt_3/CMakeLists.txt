cmake_minimum_required(VERSION 3.15)
project(signal_lib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ===== pybind11 =====
include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.11.1
)
FetchContent_MakeAvailable(pybind11)

# ===== matplot++ =====
add_subdirectory(matplotplusplus matplot_build)

# ===== Python-модуль =====
pybind11_add_module(signal_lib
    src/signal_lib.cpp
    src/signal_module.cpp
)

target_include_directories(signal_lib PRIVATE
    ${CMAKE_SOURCE_DIR}/matplotplusplus/source
)

target_link_libraries(signal_lib PRIVATE matplot)

# ===== Исполняемый файл main.cpp =====
add_executable(signal_plot main.cpp)

target_include_directories(signal_plot PRIVATE
    ${CMAKE_SOURCE_DIR}/matplotplusplus/source
)

target_link_libraries(signal_plot PRIVATE matplot)



# === Исполняемый файл main.cpp ===
add_executable(signal_plot main.cpp)

target_include_directories(signal_plot PRIVATE
    ${CMAKE_SOURCE_DIR}/matplotplusplus/source
)

target_link_libraries(signal_plot PRIVATE matplot)
