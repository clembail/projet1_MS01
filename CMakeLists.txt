cmake_minimum_required(VERSION 3.10...3.13)
project(cholesky LANGUAGES CXX)

# Trouver MPI
find_package(MPI REQUIRED)

# Exécutable
add_executable(main
    source/main.cpp
    source/jacobi_seq.cpp
)

# Inclure les headers
target_include_directories(main PRIVATE
    include
    ${MPI_INCLUDE_PATH}
)

# Activer C++17
target_compile_features(main PRIVATE cxx_std_17)

# Options de compilation
target_compile_options(main PRIVATE
    -Wall
    -fsanitize=address
    ${MPI_CXX_COMPILE_FLAGS}
)

# Options de linkage
target_link_options(main PRIVATE
    -fsanitize=address
)

# Lier la bibliothèque MPI
target_link_libraries(main PRIVATE
    ${MPI_CXX_LIBRARIES}
)