# @author Munzir Zafar
# @date Dec 12, 2018

cmake_minimum_required(VERSION 2.6)

install(DIRECTORY "3DOF-WIP" "7DOFArm" "Krang" "KrangFixedWheels" "KrangWaist" "scenes"
  DESTINATION /usr/local/share/krang/urdf
  FILES_MATCHING
  PATTERN "*.urdf"
  PATTERN "*.dae"
  PATTERN "build" EXCLUDE
  PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

# uninstall target
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
