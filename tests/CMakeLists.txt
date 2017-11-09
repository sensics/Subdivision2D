
add_executable(BasicTests BasicTests.cpp)
target_link_libraries(BasicTests PRIVATE Subdivision2D sd2d-catch-vendored)
add_test(NAME BasicTests COMMAND BasicTests)
