set(SOURCE
  ${SOURCE}
  PARENT_SCOPE
)

set(HEADERS
  ${HEADERS}
  ${CMAKE_CURRENT_SOURCE_DIR}/bioio.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/catch.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cxxopts.h
  PARENT_SCOPE
)

add_subdirectory(spdlog)
