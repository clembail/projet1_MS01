cmake_minimum_required(VERSION 3.10...3.13)
project(cholesky)

find_package(MPI REQUIRED)

# === Test Jacobi ===
add_executable(jacobi source/jacobi.cpp)
target_include_directories(jacobi PRIVATE include)
target_compile_features(jacobi PRIVATE cxx_std_17)
target_compile_options(jacobi PRIVATE -Wall -fsanitize=address)
target_link_options(jacobi PRIVATE -fsanitize=address)

# === Test Gauss-Seidel ===
add_executable(gauss-seidel source/gauss-seidel.cpp)
target_include_directories(gauss-seidel PRIVATE include)
target_compile_features(gauss-seidel PRIVATE cxx_std_17)
target_compile_options(gauss-seidel PRIVATE -Wall -fsanitize=address)
target_link_options(gauss-seidel PRIVATE -fsanitize=address)

# === Test Jacobi parallélisé ===
add_executable(jacobi_para source/jacobi_para.cpp)
target_include_directories(jacobi_para PRIVATE include ${MPI_INCLUDE_PATH})
target_compile_features(jacobi_para PRIVATE cxx_std_17)
target_compile_options(jacobi_para PRIVATE -Wall -fsanitize=address ${MPI_CXX_COMPILE_FLAGS})
target_link_options(jacobi_para PRIVATE -fsanitize=address ${MPI_CXX_LINK_FLAGS})
target_link_libraries(jacobi_para PRIVATE ${MPI_CXX_LIBRARIES})

# === Test Gauss-Seidel parallélisé ===
add_executable(gauss-seidel_para source/gauss-seidel_para.cpp)
target_include_directories(gauss-seidel_para PRIVATE include ${MPI_INCLUDE_PATH})
target_compile_features(gauss-seidel_para PRIVATE cxx_std_17)
target_compile_options(gauss-seidel_para PRIVATE -Wall -fsanitize=address ${MPI_CXX_COMPILE_FLAGS})
target_link_options(gauss-seidel_para PRIVATE -fsanitize=address ${MPI_CXX_LINK_FLAGS})
target_link_libraries(gauss-seidel_para PRIVATE ${MPI_CXX_LIBRARIES})