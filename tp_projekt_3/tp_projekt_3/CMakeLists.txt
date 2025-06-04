# CMakeList.txt : CMake project for p3, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (p3 "p3.cpp" "p3.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET p3 PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
