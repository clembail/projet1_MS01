// pour exécuter : mpirun -np 4 ./buildPrenom/jacobi_para

#include <mpi.h>
#include <iostream>
// #include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

// const int Nx = 100;             // points intérieurs
// const int Ny = 100;
// const double TOLERANCE = 1e-5;
// const int MAX_ITERATION = 25000;
const double alpha = 0.5;
const double a = 1.0;
const double b = 1.0;

double V(double y){
    return (1 - cos(2*M_PI*y/b));
}

double u0(double x, double y) {
    double U0 = 1.0;
    if (x==0){
        return U0*(1.0 + alpha*V(y));
    }

    return U0;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

        // --- DEBUT CHANGEMENT ---
    if (argc != 5) {
        // Pour MPI, n'imprimez que sur le rang 0
        // if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <Nx> <Ny> <maxIter> <tolerance>" << std::endl;
        // }
        // MPI_Finalize();
        return 1;
    }

    const int Nx = std::atoi(argv[1]);
    const int Ny = std::atoi(argv[2]);
    const int maxIter = std::atoi(argv[3]);
    const double tol = std::atof(argv[4]);
    // --- FIN CHANGEMENT ---

    const int nrows_global = Nx + 2;
    const int ncols_global = Ny + 2;

    const double dx = a / (Nx + 1);
    const double dy = b / (Ny + 1);

    MPI_Barrier(MPI_COMM_WORLD);
    double time1 = MPI_Wtime();

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
    std::vector<double> u_new = u;

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
    u_new = u;

    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double denom = 2.0 * (dx2 + dy2);
    int iteration = 0;
    double global_diff;

    do {
// 1. ÉCHANGE DES LIGNES FANTÔMES (Corrigé avec MPI_Sendrecv)
        
        // Définir les voisins, MPI_PROC_NULL gère les bords (rank 0 et size-1)
        int rank_up = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
        int rank_down = (rank < size - 1) ? rank + 1 : MPI_PROC_NULL;

        // 1. Envoyer vers le HAUT, Recevoir d'en HAUT
        MPI_Sendrecv(
            &u[idx(start_idx_loc, 0)], ncols_global, MPI_DOUBLE, rank_up, 0, // Données envoyées
            &u[idx(0, 0)], ncols_global, MPI_DOUBLE, rank_up, 0,             // Tampon de réception
            MPI_COMM_WORLD, MPI_STATUS_IGNORE
        );

        // 2. Envoyer vers le BAS, Recevoir d'en BAS
        MPI_Sendrecv(
            &u[idx(start_idx_loc + nloc - 1, 0)], ncols_global, MPI_DOUBLE, rank_down, 0, // Données envoyées
            &u[idx(start_idx_loc + nloc, 0)], ncols_global, MPI_DOUBLE, rank_down, 0,     // Tampon de réception
            MPI_COMM_WORLD, MPI_STATUS_IGNORE
        );

        // 2. CALCUL JACOBI LOCAL
        double local_diff = 0.0;
        int i_start = start_idx_loc;
        int i_end = start_idx_loc + nloc;
        if (rank == 0) i_start++;
        if (rank == size - 1) i_end--;

        for (int i = i_start; i < i_end; ++i) {
            for (int j = 1; j <= Ny; ++j) {
                u_new[idx(i, j)] = (
                    (u[idx(i+1, j)] + u[idx(i-1, j)]) * dy2 +
                    (u[idx(i, j+1)] + u[idx(i, j-1)]) * dx2
                ) / denom;
                local_diff = std::max(local_diff, std::fabs(u_new[idx(i, j)] - u[idx(i, j)]));
            }
        }
        
        MPI_Allreduce(&local_diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        std::swap(u, u_new);
        iteration++;

    } while (global_diff > tol && iteration < maxIter);

    // if (rank == 0) {
    //     std::cout << "A convergé en " << iteration << " itérations avec une erreur de " << global_diff << " avec " << size << " processus." << std::endl;
    // }

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

    // // =========================
    // //   ÉCRITURE CSV
    // // =========================

    // if (rank == 0) {
    //     std::string pathData {"data_jacobi_para.csv"};
    //     std::ofstream file(pathData);
    //     for (int i = 0; i < nrows_global; ++i) {
    //         for (int j = 0; j < ncols_global; ++j) {
    //             file << u_global[i * ncols_global + j];
    //             if (j < ncols_global - 1) file << ",";
    //         }
    //         file << "\n";
    //     }
    //     file.close();
    //     std::cout << "Résultat écrit dans " + pathData << std::endl;
    // }

    MPI_Barrier(MPI_COMM_WORLD);
    double time2 = MPI_Wtime();

    if (rank == 0) {
        // --- DEBUT CHANGEMENT ---
        // Ancien cout :
        // std::cout << "gauss-seidel_para, " << Nx << ", " << size << ", " 
        //           << time2-time1 << ", " << iteration << ", " 
        //           << global_diff << std::endl;

        // Nouveau cout (plus facile à parser) :
        std::cout << "Temps: " << (time2 - time1) << std::endl;
        std::cout << "Iterations: " << iteration << std::endl;
        std::cout << "Error: " << global_diff << std::endl;
        // --- FIN CHANGEMENT ---
    }


    MPI_Finalize();
    return 0;
}