cmake_minimum_required(VERSION 3.14...3.21)

project(
  leetcode 
  LANGUAGES CXX)

include(FetchContent)
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.6)

FetchContent_MakeAvailable(catch)

aux_source_directory(src EXECUTABLE)
add_executable(tests main.cpp ${EXECUTABLE})
target_link_libraries(tests PRIVATE Catch2::Catch2)
target_include_directories(tests PUBLIC include)
target_compile_features(tests PUBLIC cxx_std_17)
# target_compile_definitions(tests PUBLIC CATCH_CONFIG_ENABLE_BENCHMARKING)
add_test(NAME tests COMMAND tests)