// pour exécuter : mpirun -np 4 ./buildPrenom/jacobi_para

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

const int Nx = 100;             // points intérieurs
const int Ny = 100;
const double TOLERANCE = 1e-6;
const int MAX_ITERATION = 10000;
const double alpha = 0.5;

double V(double y){
    return (1 - cos(2*M_1_PI*y/1.0));
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

    const int nrows_global = Nx;
    const int ncols = Ny;

    const double dx = 1.0 / (Nx + 1);
    const double dy = 1.0 / (Ny + 1);

    // Découpe du domaine : chaque processus gère nloc lignes
    int nloc = nrows_global / size;
    int start_row = rank * nloc + 1;
    int end_row = start_row + nloc - 1;

    // Taille locale incluant les lignes fantômes
    int nrows_loc = nloc + 2;
    int ncols_loc = ncols + 2;

    // Tableau 2D contigu
    std::vector<double> u(nrows_loc * ncols_loc, 0.0);
    std::vector<double> u_new = u;

    // Fonction d’indexation 2D → 1D
    auto idx = [&](int i, int j) { return i * ncols_loc + j; };

    // Conditions aux bords gauche/droite : u = 1
    for (int i = 0; i < nrows_loc; ++i) {
        u[idx(i, 0)] = u[idx(i, ncols_loc - 1)] = 1.0;
        u_new[idx(i, 0)] = u_new[idx(i, ncols_loc - 1)] = 1.0;
    }

    // Conditions sur le bord global haut/bas
    if (rank == 0) {
        for (int j = 0; j < ncols_loc; ++j) {
            u[idx(1, j)] = 1.0;
            u_new[idx(1, j)] = 1.0;
        }
    }

    if (rank == size - 1) {
        for (int j = 0; j < ncols_loc; ++j) {
            u[idx(nloc, j)] = 1.0;
            u_new[idx(nloc, j)] = 1.0;
        }
    }

    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double denom = 2.0 * (dx2 + dy2);

    int iteration = 0;
    double global_diff;

    do {
        // 1. Échange des lignes fantômes
        if (rank > 0) {
            MPI_Send(&u[idx(1, 0)], ncols_loc, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u[idx(0, 0)], ncols_loc, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            MPI_Send(&u[idx(nloc, 0)], ncols_loc, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u[idx(nloc + 1, 0)], ncols_loc, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // 2. Calcul de Jacobi local
        double local_diff = 0.0;
        for (int i = 1; i <= nloc; ++i) {
            for (int j = 1; j <= ncols; ++j) {
                u_new[idx(i, j)] = (
                    (u[idx(i+1, j)] +
                    u[idx(i-1, j)])*dx2 +
                    (u[idx(i, j+1)] +
                    u[idx(i, j-1)])*dy2
                )/denom;
                local_diff = std::max(local_diff, std::fabs(u_new[idx(i, j)] - u[idx(i, j)]));
            }
        }

        // 3. Calcul de l’erreur globale
        MPI_Allreduce(&local_diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // 4. Mise à jour
        std::swap(u, u_new);
        ++iteration;

    } while (global_diff > TOLERANCE && iteration < MAX_ITERATION);

    if (rank == 0) {
        std::cout << "À convergé en " << iteration << " itérations pour = " << global_diff << " avec " << size << " processus." << std::endl;
    }

    // =========================
    //   RASSEMBLEMENT MPI
    // =========================

    int local_size = nloc * ncols_loc;
    std::vector<double> u_send(nloc * ncols_loc);

    // On enlève les lignes fantômes
    for (int i = 0; i < nloc; ++i) {
        std::copy(u.begin() + idx(i + 1, 0),
                  u.begin() + idx(i + 2, 0),
                  u_send.begin() + i * ncols_loc);
    }

    std::vector<double> u_global;
    if (rank == 0) {
        u_global.resize(nrows_global * ncols_loc);
    }

    MPI_Gather(u_send.data(), local_size, MPI_DOUBLE,
               u_global.data(), local_size, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // =========================
    //   ÉCRITURE CSV
    // =========================
    if (rank == 0) {
        std::string pathData {"data_jacobi_para.csv"};
        std::ofstream file(pathData);
        for (int i = 0; i < nrows_global; ++i) {
            for (int j = 0; j < ncols_loc; ++j) {
                file << u_global[i * ncols_loc + j];
                if (j < ncols_loc - 1) file << ",";
            }
            file << "\n";
        }
        file.close();
        std::cout << "Résultat écrit dans " + pathData << std::endl;
    }

    MPI_Finalize();
    return 0;
}