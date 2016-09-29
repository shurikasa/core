add_executable(bbp_exec main.cc)
target_link_libraries(bbp_exec core)
bob_export_target(bbp_exec)
bob_end_subdir()

