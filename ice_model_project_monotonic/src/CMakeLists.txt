
# list of cpp files to compile
set(ICE_MODEL_SRC
    #matrix.cpp
    #newton_raphson_solver.cpp 
    FDM_solver.cpp
    FDM_ff.cpp
    FDM_icelens.cpp
    Ice_model.cpp
)

# list of header files to compile 
set(ICE_MODEL_HEADERS
    #common.hpp
    #matrix_inline.hpp
    #vector.hpp
    #base_matrix.hpp
    #matrix.hpp
    #nonlinear_solver.hpp
    #newton_raphson.hpp 
    Ice_model.hpp
    FDM_solver.hpp
    FDM_ff.hpp
    FDM_icelens.hpp
    input_file.hpp
)

set(ICE_MODEL_DEPEND_FILES ${ICE_MODEL_SRC} ${ICE_MODEL_HEADERS} PARENT_SCOPE)

add_library(ice_model ${ICE_MODEL_SRC})

# if linking to other libraries, e.g. MPI
target_link_libraries(ice_model ${ARMADILLO_LIBRARIES} ${EIGEN_LIBRARIES})
