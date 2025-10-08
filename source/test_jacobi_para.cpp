#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

const int N = 100;             // points intérieurs
const double TOLERANCE = 1e-6;
const int MAX_ITERATION = 10000;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int nrows_global = N;
    const int ncols = N;

    const double h = 1.0 / (N + 1);

    // Découpe du domaine : chaque processus gère nloc lignes
    int nloc = nrows_global / size;
    int start_row = rank * nloc + 1;
    int end_row = start_row + nloc - 1;

    // Ajouter 2 lignes fantômes (top + bottom)
    std::vector<std::vector<double>> u(nloc + 2, std::vector<double>(ncols + 2, 0.0));
    std::vector<std::vector<double>> u_new = u;

    // Conditions aux bords gauche/droite : u = 1
    for (int i = 0; i < nloc + 2; ++i) {
        u[i][0] = u[i][ncols + 1] = 1.0;
        u_new[i][0] = u_new[i][ncols + 1] = 1.0;
    }

    // Conditions sur le bord global haut/bas
    if (rank == 0) {
        for (int j = 0; j <= ncols + 1; ++j) {
            u[1][j] = 1.0;
            u_new[1][j] = 1.0;
        }
    }

    if (rank == size - 1) {
        for (int j = 0; j <= ncols + 1; ++j) {
            u[nloc][j] = 1.0;
            u_new[nloc][j] = 1.0;
        }
    }

    int iteration = 0;
    double global_diff;

    do {
        // 1. Échange des lignes fantômes
        if (rank > 0) {
            MPI_Send(u[1].data(), ncols + 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(u[0].data(), ncols + 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            MPI_Send(u[nloc].data(), ncols + 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(u[nloc + 1].data(), ncols + 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // 2. Calcul de Jacobi local
        double local_diff = 0.0;
        for (int i = 1; i <= nloc; ++i) {
            for (int j = 1; j <= ncols; ++j) {
                u_new[i][j] = 0.25 * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]);
                local_diff = std::max(local_diff, std::fabs(u_new[i][j] - u[i][j]));
            }
        }

        // 3. Calcul de l’erreur globale
        MPI_Allreduce(&local_diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // 4. Mise à jour
        std::swap(u, u_new);
        ++iteration;

    } while (global_diff > TOLERANCE && iteration < MAX_ITERATION);

    if (rank == 0) {
        std::cout << "Converged in " << iteration << " iterations with diff = " << global_diff << std::endl;
    }

    MPI_Finalize();
    return 0;
}
