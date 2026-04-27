project(02_HMM_tests)

include_directories("${GTEST_SRC_DIR}/include")

function(add_gtest TEST_NAME)
    # Parse the additional arguments
    set(MULTI_ARGUMENTS LIBRARIES SOURCE_FILES)
    cmake_parse_arguments(MY_ARGS "" "" "${MULTI_ARGUMENTS}" ${ARGN})

    # Setup the executable
    add_executable(${TEST_NAME} "${TEST_NAME}.cpp" ${MY_ARGS_SOURCE_FILES})
    target_link_libraries(${TEST_NAME} gtest gtest_main test_driver ${MY_ARGS_LIBRARIES})

    # Add the test to CTest
    add_test(NAME ${TEST_NAME} WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests" COMMAND ${TEST_NAME})
endfunction()

add_library(test_driver main_test.cpp)

add_gtest(log_prob_test LIBRARIES HMM)
add_gtest(HMM_test LIBRARIES HMM)
add_gtest(distribution_test LIBRARIES HMM)