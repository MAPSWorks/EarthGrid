

include_directories("${PROJECT_SOURCE_DIR}/include")


aux_source_directory("${PROJECT_SRC_DIR}/base" BASE_SRCS)
aux_source_directory("${PROJECT_SRC_DIR}/core" CORE_SRCS)
aux_source_directory("${PROJECT_SRC_DIR}/object" OBJECT_SRCS)

add_library(core  ${CORE_SRCS} )
add_library(object ${OBJECT_SRCS})