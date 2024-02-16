add_executable(example example.cpp)
target_link_libraries(example ice_model)

configure_file(vertical_wall.inp vertical_wall.inp COPYONLY)

