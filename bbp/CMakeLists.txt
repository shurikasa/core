add_executable(bbp_exec main.cc)
target_link_libraries(bbp_exec core)
add_executable(iterate_neuron_exec iterate_neuron_in_cube.cc)
target_link_libraries(iterate_neuron_exec core)
add_executable(ma_exec ma.cc)
target_link_libraries(ma_exec core)
bob_export_target(bbp_exec)
bob_end_subdir()

