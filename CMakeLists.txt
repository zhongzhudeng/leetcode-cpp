cmake_minimum_required(VERSION 3.14..)

project(
  leetcode 
  LANGUAGES CXX
)

find_package(Catch2 3 REQUIRED)

aux_source_directory(src EXECUTABLE)
add_executable(tests ${EXECUTABLE})
target_compile_features(tests PUBLIC cxx_std_20)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)
target_include_directories(tests PUBLIC include)
add_test(NAME tests COMMAND tests)