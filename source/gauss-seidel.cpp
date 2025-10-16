// pour exécuter : ./buildPrenom/gauss-seidel

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

const int Nx = 50;       // Nombre de points intérieurs en x
const int Ny = 50;       // Nombre de points intérieurs en y
const int maxIter = 10000;
const double tol = 1e-6;
const double alpha = 0.5;

double V(double y){
    return (1 - cos(2*M_1_PI*y/1.0));
}

double u0(double x, double y) {
    double U0 = 1.0;
    if (x==0){
        return U0*(1.0 + alpha*V(y));
    }

    return U0;
}

int main() {
    double dx = 1.0 / (Nx + 1);
    double dy = 1.0 / (Ny + 1);

    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double denom = 2.0 * (dx2 + dy2);

    // Création d'une grille Nx * Ny avec les bords
    std::vector<std::vector<double>> u(Nx + 2, std::vector<double>(Ny + 2, 0.0));

    // Condition au bord :
    for (int i = 0; i <= Nx + 1; ++i) {
        double x = i * dx;
        u[i][0] = u0(x, 0);             // bas
        u[i][Ny + 1] = u0(x, 1);        // haut
    }

    for (int j = 0; j <= Ny + 1; ++j) {
        double y = j * dy;
        u[0][j] = u0(0, y);             // gauche
        u[Nx + 1][j] = u0(1, y);        // droite
    }

    // f(i,j) = 0 pour tout le domaine (Laplace homogène)
    std::vector<std::vector<double>> f(Nx + 2, std::vector<double>(Ny + 2, 0.0));

    int iter = 0;
    double diff;

    do {
        diff = 0.0;

        // Mise à jour selon Gauss-Seidel
        for (int i = 1; i <= Nx; ++i) {
            for (int j = 1; j <= Ny; ++j) {
                double old_u = u[i][j];
                u[i][j] = ((u[i+1][j] + u[i-1][j]) * dy2 +
                               (u[i][j+1] + u[i][j-1]) * dx2 - f[i][j]*dx2*dy2) / denom;
                diff = std::max(diff, std::fabs(u[i][j] - old_u));
            }
        }

        ++iter;
    } while (diff > tol && iter < maxIter);

    std::cout << "Convergence en " << iter << " iterations, avec diff = " << diff << "\n";

    // Affichage d'une partie de la solution (optionnel)
    for (int i = 0; i <= Nx + 1; i += Nx/10) {
        for (int j = 0; j <= Ny + 1; j += Ny/10) {
            std::cout << u[i][j] << " ";
        }
        std::cout << "\n";
    }

    std::ofstream file("data_gauss_seidel.csv");
    for(int i = 0; i<=Nx+1 ; i++){
        for(int j = 0; j<=Ny+1 ; j++){
            file << u[i][j];
            if (j < Ny + 1){ file << ",";}
        }
        file << "\n";
    }

    return 0;
}
