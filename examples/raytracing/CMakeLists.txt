add_executable(ArborX_Raytracing.exe example_raytracing.cpp)
target_link_libraries(ArborX_Raytracing.exe ArborX::ArborX)
add_test(NAME ArborX_Raytracing_Example COMMAND ./ArborX_raytracing.exe)
