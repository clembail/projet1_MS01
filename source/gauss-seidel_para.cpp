// pour exécuter : mpirun -np 4 ./buildPrenom/jacobi_para

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

const int Nx = 50;             // points intérieurs
const int Ny = 50;
const double TOLERANCE = 1e-6;
const int MAX_ITERATION = 10000;
const double alpha = 0.5;
const double a = 1.0;
const double b = 1.0;

double V(double y){
    return (1 - cos(2*M_1_PI*y/b));
}

double u0(double x, double y) {
    // Condition au bord, ici u0 = 1
    if (x==0){
        return 1.0*(1.0 + alpha*V(y));
    }

    return 1.0;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int nrows_global = Nx + 2;
    const int ncols_global = Ny + 2;

    const double dx = a / (Nx + 1);
    const double dy = b / (Ny + 1);

    // Découpe du domaine
    int nloc = nrows_global / size;
    int remainder = nrows_global % size;
    if (rank < remainder) {
        nloc++;
    }

    int start_row_global = rank * (nrows_global / size) + std::min(rank, remainder);

    // Taille locale incluant les lignes fantômes
    int nrows_loc = nloc;
    if (rank > 0) nrows_loc++;
    if (rank < size - 1) nrows_loc++;
    
    int start_idx_loc = (rank > 0) ? 1 : 0;

    std::vector<double> u(nrows_loc * ncols_global, 0.0);

    auto idx = [&](int i, int j) { return i * ncols_global + j; };

    // =============================================================
    // INITIALISATION DES CONDITIONS AUX BORDS (Logique correcte)
    // =============================================================
    for (int i_loc = 0; i_loc < nloc; ++i_loc) {
        int i_global = start_row_global + i_loc;
        double x = i_global * dx;
        
        // Bords BAS (y=0) et HAUT (y=b)
        u[idx(start_idx_loc + i_loc, 0)] = u0(x, 0.0);
        u[idx(start_idx_loc + i_loc, ncols_global - 1)] = u0(x, b);
    }
    if (rank == 0) {
        for (int j = 0; j < ncols_global; ++j) {
            u[idx(0, j)] = u0(0.0, j * dy);
        }
    }
    if (rank == size - 1) {
        int last_data_row = start_idx_loc + nloc - 1;
        for (int j = 0; j < ncols_global; ++j) {
            u[idx(last_data_row, j)] = u0(a, j * dy);
        }
    }

    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double denom = 2.0 * (dx2 + dy2);
    int iteration = 0;
    double global_diff;

    do {
        // 1. ÉCHANGE DES LIGNES FANTÔMES
        if (rank > 0) {
            MPI_Send(&u[idx(start_idx_loc, 0)], ncols_global, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u[idx(0, 0)], ncols_global, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            MPI_Send(&u[idx(start_idx_loc + nloc - 1, 0)], ncols_global, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u[idx(start_idx_loc + nloc, 0)], ncols_global, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // 2. MISE A JOUR DES BOULES ROUGES
        double local_diff = 0.0;
        int i_start = start_idx_loc;
        int i_end = start_idx_loc + nloc;
        if (rank == 0) i_start++;
        if (rank == size - 1) i_end--;

        for (int i = i_start; i < i_end; ++i) {
            for (int j = 1; j <= Ny; ++j) {
                if ((i+j)%2==0){
                    double u_old = u[idx(i,j)];
                    u[idx(i, j)] = (
                        (u[idx(i+1, j)] + u[idx(i-1, j)]) * dy2 +
                        (u[idx(i, j+1)] + u[idx(i, j-1)]) * dx2
                    ) / denom;
                    local_diff = std::max(local_diff, std::fabs(u_old - u[idx(i, j)]));}
            }
        }

        // 3. COMMUNICATION DE LA MISE A JOUR DES BOULES ROUGES
        if (rank > 0) {
            MPI_Send(&u[idx(start_idx_loc, 0)], ncols_global, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u[idx(0, 0)], ncols_global, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            MPI_Send(&u[idx(start_idx_loc + nloc - 1, 0)], ncols_global, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u[idx(start_idx_loc + nloc, 0)], ncols_global, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // 2. MISE A JOUR DES BOULES NOIRES
        // double local_diff = 0.0;
        // int i_start = start_idx_loc;
        // int i_end = start_idx_loc + nloc;
        // if (rank == 0) i_start++;
        // if (rank == size - 1) i_end--;

        for (int i = i_start; i < i_end; ++i) {
            for (int j = 1; j <= Ny; ++j) {
                if ((i+j)%2==1){
                    double u_old = u[idx(i,j)];
                    u[idx(i, j)] = (
                        (u[idx(i+1, j)] + u[idx(i-1, j)]) * dy2 +
                        (u[idx(i, j+1)] + u[idx(i, j-1)]) * dx2
                    ) / denom;
                    local_diff = std::max(local_diff, std::fabs(u_old - u[idx(i, j)]));}
            }
        }
        
        MPI_Allreduce(&local_diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        iteration++;

    } while (global_diff > TOLERANCE && iteration < MAX_ITERATION);

    if (rank == 0) {
        std::cout << "A convergé en " << iteration << " itérations avec une erreur de " << global_diff << " avec " << size << " processus." << std::endl;
    }

    // ===============================================
    //   RASSEMBLEMENT AVEC MPI_GATHERV
    // ===============================================

    // Chaque processus prépare ses données (sans lignes fantômes)
    std::vector<double> u_send(nloc * ncols_global);
    std::copy(u.begin() + idx(start_idx_loc, 0), 
              u.begin() + idx(start_idx_loc + nloc, 0), 
              u_send.begin());

    std::vector<int> recvcounts, displs;
    std::vector<double> u_global;

    if (rank == 0) {
        recvcounts.resize(size);
        displs.resize(size);
        u_global.resize(nrows_global * ncols_global);
    }

    // Le rang 0 collecte le nombre de lignes de chaque processus
    int send_count = nloc * ncols_global;
    MPI_Gather(&send_count, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Le rang 0 calcule les déplacements
    if (rank == 0) {
        displs[0] = 0;
        for (int i = 1; i < size; i++) {
            displs[i] = displs[i-1] + recvcounts[i-1];
        }
    }

    // Rassemblement des données avec Gatherv
    MPI_Gatherv(u_send.data(), send_count, MPI_DOUBLE,
                u_global.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // =========================
    //   ÉCRITURE CSV
    // =========================

    if (rank == 0) {
        std::string pathData {"data_gauss-seidel_para.csv"};
        std::ofstream file(pathData);
        for (int i = 0; i < nrows_global; ++i) {
            for (int j = 0; j < ncols_global; ++j) {
                file << u_global[i * ncols_global + j];
                if (j < ncols_global - 1) file << ",";
            }
            file << "\n";
        }
        file.close();
        std::cout << "Résultat écrit dans " + pathData << std::endl;
    }

    MPI_Finalize();
    return 0;
}