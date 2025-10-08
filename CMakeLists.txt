cmake_minimum_required(VERSION 3.10...3.13)
project(cholesky)

find_package(MPI REQUIRED)

# === Test Jacobi ===
add_executable(test_jacobi source/test_jacobi.cpp)
target_include_directories(test_jacobi PRIVATE include)
target_compile_features(test_jacobi PRIVATE cxx_std_17)
target_compile_options(test_jacobi PRIVATE -Wall -fsanitize=address)
target_link_options(test_jacobi PRIVATE -fsanitize=address)

# === Test Gauss-Seidel ===
add_executable(test_gauss_seidel source/test_gauss_seidel.cpp)
target_include_directories(test_gauss_seidel PRIVATE include)
target_compile_features(test_gauss_seidel PRIVATE cxx_std_17)
target_compile_options(test_gauss_seidel PRIVATE -Wall -fsanitize=address)
target_link_options(test_gauss_seidel PRIVATE -fsanitize=address)

# === Exécutable 2 (version MPI) ===
add_executable(mpi source/test_jacobi_para.cpp)
target_include_directories(mpi PRIVATE include ${MPI_INCLUDE_PATH})
target_compile_features(mpi PRIVATE cxx_std_17)
target_compile_options(mpi PRIVATE -Wall -fsanitize=address ${MPI_CXX_COMPILE_FLAGS})
target_link_options(mpi PRIVATE -fsanitize=address ${MPI_CXX_LINK_FLAGS})
target_link_libraries(mpi PRIVATE ${MPI_CXX_LIBRARIES})