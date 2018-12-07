

include_directories("${PROJECT_SOURCE_DIR}/include")


aux_source_directory("${PROJECT_SRC_DIR}/third_party/absl/base" THIRDPARTY_SRCS)
aux_source_directory("${PROJECT_SRC_DIR}/third_party/absl/base/internal" THIRDPARTY_SRCS)
aux_source_directory("${PROJECT_SRC_DIR}/third_party/absl/numeric" THIRDPARTY_SRCS)
aux_source_directory("${PROJECT_SRC_DIR}/third_party/absl/strings" THIRDPARTY_SRCS)
aux_source_directory("${PROJECT_SRC_DIR}/third_party/absl/strings/internal" THIRDPARTY_SRCS)

aux_source_directory("${PROJECT_SRC_DIR}/util/bits" UTIL_SRCS)
aux_source_directory("${PROJECT_SRC_DIR}/util/coding" UTIL_SRCS)
aux_source_directory("${PROJECT_SRC_DIR}/util/math" UTIL_SRCS)
aux_source_directory("${PROJECT_SRC_DIR}/util/math/exactfloat" UTIL_SRCS)
aux_source_directory("${PROJECT_SRC_DIR}/util/unit" UTIL_SRCS)

aux_source_directory("${PROJECT_SRC_DIR}/base" BASE_SRCS)
aux_source_directory("${PROJECT_SRC_DIR}/core" CORE_SRCS)
aux_source_directory("${PROJECT_SRC_DIR}/cellid" CELLID_SRCS)
aux_source_directory("${PROJECT_SRC_DIR}/object" OBJECT_SRCS)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BIN_DIR})
set(LIBRARY_OUTPUT_PATH  ${PROJECT_LIB_DIR})

set(CMAKE_C_FLAGS "w")


add_library(third_party ${THIRDPARTY_SRCS})
add_library(util ${UTIL_SRCS})
add_library(base ${BASE_SRCS})
add_library(core  ${CORE_SRCS} )
add_library(cellid  ${CELLID_SRCS} )
add_library(object ${OBJECT_SRCS})