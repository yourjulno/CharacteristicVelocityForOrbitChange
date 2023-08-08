cmake_minimum_required(VERSION 3.17)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow -Wall -Werror")
file(GLOB_RECURSE SOURCES *.h *.cpp *.hpp)
message(${SOURCES})
#find_package(Eigen3  REQUIRED NO_MODULE)
add_library(${PROJECT_NAME} ${SOURCES}

)

message(lalaa)
#set(EIGEN_LIB Eigen3)
#
#target_link_libraries(Velocity Eigen3::Eigen)