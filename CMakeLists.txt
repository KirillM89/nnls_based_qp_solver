cmake_minimum_required(VERSION 3.14)
project(nnls_qp)
set(GTEST_CREATE_SHARED_LIBRARY 1)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(PUB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_library(${PROJECT_NAME} STATIC)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 20)
target_sources(${PROJECT_NAME} PRIVATE
    ${SRC_DIR}/log.cpp
    ${SRC_DIR}/timers.cpp
    ${SRC_DIR}/utils.cpp
    ${SRC_DIR}/decorator.cpp
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
    ${PUB_DIR}/decorator.h
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
target_include_directories(smoke PRIVATE
    ${PUB_DIR}
    #${CMAKE_CURRENT_SOURCE_DIR}/build/Debug
    #${CMAKE_CURRENT_SOURCE_DIR}/build/Release
)

set(GTEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests/gtest/install)
if(NOT EXISTS "${GTEST_DIR}")
    message(ERROR "Directory ${GTEST_DIR} does not exist!")
endif()
enable_testing()
add_executable(unit_tests)
target_link_directories(unit_tests PRIVATE
    ${GTEST_DIR}/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/build/Debug
    #${CMAKE_CURRENT_SOURCE_DIR}/build/Release
)
target_include_directories(unit_tests PRIVATE
    ${GTEST_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

#add_library(gt SHARED IMPORTED)
#set_target_properties(gt PROPERTIES
#    IMPORTED_LOCATION ${GTEST_DIR}/lib/gtest.lib
#    IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
#)
#add_library(gtm SHARED IMPORTED)
#set_target_properties(gtm PROPERTIES
#     IMPORTED_LOCATION ${GTEST_DIR}/lib/gtest_main.lib
#     IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
#)
target_link_libraries(unit_tests PRIVATE ${PROJECT_NAME} ${GTEST_DIR}/lib/gtest.lib)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)

#add_executable(examples)
#set_property(TARGET examples PROPERTY CXX_STANDARD 20)
#target_sources(examples PRIVATE
#    ${CMAKE_CURRENT_SOURCE_DIR}/examples/main.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/examples/tools.h
#    ${CMAKE_CURRENT_SOURCE_DIR}/examples/tools.cpp)
#target_link_libraries(examples PRIVATE ${PROJECT_NAME})
#target_include_directories(examples PRIVATE ${PUB_DIR})

