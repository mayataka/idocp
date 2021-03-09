# enable google test
include(cmake/gtest.cmake)

# build test helper 
add_library(
  test_helper
  SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/test_helper/test_helper.cpp
)
target_include_directories(
  test_helper
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/test_helper
  ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(
  test_helper
  PUBLIC 
  ${PROJECT_NAME} 
)

# macto for testing
macro(add_idocp_test TESTNAME)
  add_executable(
    ${TESTNAME} 
    ${TESTNAME}.cpp
  )
  target_include_directories(
    ${TESTNAME} 
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/test/test_helper
    ${PINOCCHIO_INCLUDE_DIRS}
  )
  target_link_libraries(
    ${TESTNAME} 
    PRIVATE
    GTest::GTest
    GTest::GMock
    ${PROJECT_NAME}
    test_helper
  )
  add_dependencies(
    ${TESTNAME} 
    googletest
  )
  add_test(
    NAME ${TESTNAME} 
    COMMAND $<TARGET_FILE:${TESTNAME}>
  )
endmacro()

# copy urdf files
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/urdf/iiwa14/iiwa14.urdf
  ${CMAKE_CURRENT_BINARY_DIR}/urdf/iiwa14/iiwa14.urdf
  COPYONLY
)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/urdf/anymal/anymal.urdf
  ${CMAKE_CURRENT_BINARY_DIR}/urdf/anymal/anymal.urdf
  COPYONLY
)

# add tests
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/robot)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cost)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/constraints)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/impulse)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/hybrid)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ocp)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unocp)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/line_search)