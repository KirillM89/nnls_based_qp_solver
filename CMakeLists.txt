cmake_minimum_required(VERSION 3.14)
project(nnls_qp)
#set(ROOT "C:/Users/m00829527/nnls_based_qp_solver")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_CXX_FLAGS -g)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(PUB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/public)
add_library(${PROJECT_NAME} SHARED)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 20)
target_sources(${PROJECT_NAME} PRIVATE
    ${SRC_DIR}/log.cpp
    ${SRC_DIR}/timers.cpp
    ${SRC_DIR}/utils.cpp
    ${SRC_DIR}/decorators.cpp
    ${SRC_DIR}/core.cpp
    ${SRC_DIR}/linSolvers.cpp
    ${SRC_DIR}/callback.cpp
    ${SRC_DIR}/log.h
    ${SRC_DIR}/timers.h
    ${SRC_DIR}/utils.h
    ${SRC_DIR}/types.h
    ${SRC_DIR}/decorators.h
    ${SRC_DIR}/linSolvers.h
    ${SRC_DIR}/core.h
    ${SRC_DIR}/callback.h
)

