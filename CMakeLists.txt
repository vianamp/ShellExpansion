
cmake_minimum_required(VERSION 2.8)

PROJECT(ShellExpansion)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(ShellExpansion MACOSX_BUNDLE main _cc_finder _database _eblob)

if(VTK_LIBRARIES)
  target_link_libraries(ShellExpansion ${VTK_LIBRARIES})
else()
  target_link_libraries(ShellExpansion vtkHybrid vtkWidgets)
endif()

