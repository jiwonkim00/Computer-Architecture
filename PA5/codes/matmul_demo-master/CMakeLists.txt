cmake_minimum_required(VERSION 3.5)

project(matmul_demo C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "-O0  -Wall")

add_executable(matmul_demo main.c utils.h matmul.c utils.c)

set(LIBS pthread)

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(LIBS ${LIBS} rt)
endif ()

target_link_libraries(matmul_demo ${LIBS})
