include(ParseAndAddCatchTests)

add_executable(BasicTests BasicTests.cpp Container.cpp)
target_link_libraries(BasicTests PRIVATE Subdivision2D sd2d-catch-vendored)
set_property(TARGET BasicTests PROPERTY FOLDER Tests)
ParseAndAddCatchTests(BasicTests)
