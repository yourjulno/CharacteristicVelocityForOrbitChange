cmake_minimum_required(VERSION 3.17)
project(CharacteristicVelocityForOrbitChange)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow -Wall -Werror")
file(GLOB_RECURSE SOURCES *.h *.cpp *.hpp)

#find_package(Eigen3  REQUIRED NO_MODULE)
add_library(CharacteristicVelocityForOrbitChange ${SOURCES})

#set(EIGEN_LIB Eigen3)
#
#target_link_libraries(Velocity Eigen3::Eigen)