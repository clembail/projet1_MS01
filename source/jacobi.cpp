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

    const int largeur = Ny + 2; 
    auto idx = [largeur](int i, int j) {
        return i * largeur + j;
    };

    const int size = (Nx + 2) * (Ny + 2);
    vector<double> u(size, 0.0);
    vector<double> u_new(size, 0.0);

    // Initialiser les conditions au bord avec la fonction idx
    for (int i = 0; i <= Nx + 1; ++i) {
        double x = i * dx;
        u[idx(i, 0)] = u0(x, 0);             // bas
        u[idx(i, Ny + 1)] = u0(x, 1);        // haut
        u_new[idx(i, 0)] = u[idx(i, 0)];
        u_new[idx(i, Ny + 1)] = u[idx(i, Ny + 1)];
    }

    for (int j = 0; j <= Ny + 1; ++j) {
        double y = j * dy;
        u[idx(0, j)] = u0(0, y);             // gauche
        u[idx(Nx + 1, j)] = u0(1, y);        // droite
        u_new[idx(0, j)] = u[idx(0, j)];
        u_new[idx(Nx + 1, j)] = u[idx(Nx + 1, j)];
    }

    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double denom = 2.0 * (dx2 + dy2);

    int iteration = 0;
    double error;

    do {
        error = 0.0;

        for (int i = 1; i <= Nx; ++i) {
            for (int j = 1; j <= Ny; ++j) {
                int k = idx(i, j); 

                u_new[k] = ((u[k + largeur] + u[k - largeur]) * dy2 +
                            (u[k + 1]    + u[k - 1])    * dx2) / denom;

                error = max(error, fabs(u_new[k] - u[k]));
            }
        }

        u.swap(u_new);

        iteration++;

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
