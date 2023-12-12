cmake_minimum_required(VERSION 3.26)

project(matrix)

set(CMAKE_CXX_STANDARD          20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address,alignment,bool,bounds,enum,float-cast-overflow,float-divide-by-zero,integer-divide-by-zero,nonnull-attribute,null,object-size,return,returns-nonnull-attribute,shift,signed-integer-overflow,undefined,unreachable,vla-bound,vptr")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(GTest REQUIRED)
enable_testing()

add_executable(matrix ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable(matrix_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/unit_test/test.cpp)

target_link_libraries(matrix_test PRIVATE ${GTEST_LIBRARIES})

target_include_directories(matrix PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(matrix_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
