// pour exécuter : ./buildPrenom/jacobi

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>

using namespace std;

// const int Nx = 800;       // Nombre de points intérieurs en x
// const int Ny = 800;       // Nombre de points intérieurs en y
// const int maxIter = 50000;
// const double tol = 1e-4;
const double alpha = 0.5;

double V(double y){
    return (1 - cos(2*M_PI*y/1.0));
}

double u0(double x, double y) {
    double U0 = 1.0;
    if (x==0){
        return U0*(1.0 + alpha*V(y));
    }

    return U0;
}

int main(int argc, char** argv) {

    if (argc != 5) {
            std::cerr << "Usage: " << argv[0] << " <Nx> <Ny> <maxIter> <tolerance>" << std::endl;
        return 1;
    }

    const int Nx = std::atoi(argv[1]);
    const int Ny = std::atoi(argv[2]);
    const int maxIter = std::atoi(argv[3]);
    const double tol = std::atof(argv[4]);

    double dx = 1.0 / (Nx + 1);
    double dy = 1.0 / (Ny + 1);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Matrices pour les itérations de Jacobi
    vector<vector<double>> u(Nx + 2, vector<double>(Ny + 2, 0.0));
    vector<vector<double>> u_new(Nx + 2, vector<double>(Ny + 2, 0.0));

    // Initialiser les conditions au bord
    for (int i = 0; i <= Nx + 1; ++i) {
        double x = i * dx;
        u[i][0] = u0(x, 0);             // bas
        u[i][Ny + 1] = u0(x, 1);        // haut
        u_new[i][0] = u[i][0];
        u_new[i][Ny + 1] = u[i][Ny + 1];
    }

    for (int j = 0; j <= Ny + 1; ++j) {
        double y = j * dy;
        u[0][j] = u0(0, y);             // gauche
        u[Nx + 1][j] = u0(1, y);        // droite
        u_new[0][j] = u[0][j];
        u_new[Nx + 1][j] = u[Nx + 1][j];
    }

    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double denom = 2.0 * (dx2 + dy2);

    int iteration = 0;
    double error;

    do {
        error = 0.0;

        // Jacobi iteration: update interior points
        for (int i = 1; i <= Nx; ++i) {
            for (int j = 1; j <= Ny; ++j) {
                u_new[i][j] = ((u[i+1][j] + u[i-1][j]) * dy2 +
                               (u[i][j+1] + u[i][j-1]) * dx2) / denom;

                error = max(error, fabs(u_new[i][j] - u[i][j]));
            }
        }

        // Swap u and u_new
        u.swap(u_new);

        iteration++;

        // if (iteration % 1000 == 0)
        //     cout << "Iteration " << iteration << ", error = " << error << endl;

    } while (error > tol && iteration < maxIter);


    // // ÉCRITURE CSV
    // std::ofstream file("data/data_jacobi.csv");
    // for(int i = 0; i<=Nx+1 ; i++){
    //     for(int j = 0; j<=Ny+1 ; j++){
    //         file << u[i][j];
    //         if (j < Ny + 1){ file << ",";}
    //     }
    //     file << "\n";
    // }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    std::cout << "Temps: " << (duration * 0.000001) << std::endl;
    std::cout << "Iterations: " << iteration << std::endl;
    std::cout << "Error: " << error << std::endl;

    return 0;
}
