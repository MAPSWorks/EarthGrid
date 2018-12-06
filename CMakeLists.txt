cmake_minimum_required(VERSION 3.1)

project(EarthGrid)

set(PROJECT_INC_DIR "${PROJECT_SOURCE_DIR}/include")
set(PROJECT_SRC_DIR "${PROJECT_SOURCE_DIR}/src")
set(PROJECT_BIN_DIR "${PROJECT_SOURCE_DIR}/bin")
set(PROJECT_LIB_DIR "${PROJECT_SOURCE_DIR}/lib")
set(PROJECT_EXM_DIR "${PROJECT_SOURCE_DIR}/examples")
set(PROJECT_RC_DIR "${PROJECT_SOURCE_DIR}/resource")

include_directories(${PROJECT_INC_DIR})

add_subdirectory(${PROJECT_EXM_DIR})
add_subdirectory(${PROJECT_SRC_DIR})
