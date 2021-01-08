add_executable(ArborX_Callback.exe example_callback.cpp)
target_link_libraries(ArborX_Callback.exe ArborX::ArborX)
add_test(NAME ArborX_Callback_Example COMMAND ./ArborX_Callback.exe)

add_executable(ArborX_Callback_Raytracing.exe example_callback_raytracing.cpp)
target_link_libraries(ArborX_Callback_Raytracing.exe ArborX::ArborX)
add_test(NAME ArborX_Callback_Raytracing_Example COMMAND ./ArborX_Callback_raytracing.exe)
