add_subdirectory(action)
add_subdirectory(common)
add_subdirectory(distribution)
add_subdirectory(lattice)
add_subdirectory(montecarlo)
add_subdirectory(mpi)
add_subdirectory(qoi)
add_subdirectory(sampler)

if(MPI_CXX_FOUND)
  add_compile_definitions(USE_MPI)
endif()

add_executable(driver_qm driver_qm.cc ${SOURCES})
target_link_libraries (driver_qm Eigen3::Eigen GSL::gsl ${MPI_CXX_LIBRARIES})

add_executable(driver_qft driver_qft.cc ${SOURCES})
target_link_libraries (driver_qft Eigen3::Eigen GSL::gsl ${MPI_CXX_LIBRARIES})

add_executable(test_distribution test_distribution.cc ${SOURCES})
target_link_libraries (test_distribution Eigen3::Eigen GSL::gsl ${MPI_CXX_LIBRARIES})

add_executable(test_fastbessel test_fastbessel.cc ${SOURCES})
target_link_libraries (test_fastbessel Eigen3::Eigen GSL::gsl ${MPI_CXX_LIBRARIES})

add_executable(test_schwinger_fillin_distribution test_schwinger_fillin_distribution.cc ${SOURCES})
target_link_libraries (test_schwinger_fillin_distribution Eigen3::Eigen GSL::gsl ${MPI_CXX_LIBRARIES})
