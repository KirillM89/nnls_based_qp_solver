cmake_minimum_required(VERSION 3.14)
project(nnls_qp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_CXX_FLAGS -g)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(PUB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
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
    ${SRC_DIR}/linSolvers.h
    ${SRC_DIR}/core.h
    ${SRC_DIR}/core_types.h
)
target_sources(${PROJECT_NAME} PUBLIC
    ${PUB_DIR}/decorators.h
    ${PUB_DIR}/types.h
    ${PUB_DIR}/callback.h
)
target_include_directories(${PROJECT_NAME} PRIVATE ${PUB_DIR})
add_executable(smoke)
set_property(TARGET smoke PROPERTY CXX_STANDARD 20)
target_sources(smoke PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/smoke/main.cpp
)
target_link_libraries(smoke PRIVATE ${PROJECT_NAME})
target_include_directories(smoke PRIVATE ${PUB_DIR})

add_executable(examples)
set_property(TARGET examples PROPERTY CXX_STANDARD 20)
target_sources(examples PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/examples/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/examples/tools.h
    ${CMAKE_CURRENT_SOURCE_DIR}/examples/tools.cpp)
target_link_libraries(examples PRIVATE ${PROJECT_NAME})
target_include_directories(examples PRIVATE ${PUB_DIR})

