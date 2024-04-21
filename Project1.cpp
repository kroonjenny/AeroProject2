#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

// Define grid parameters
const int fact = 1; // Factor 
const int num_points_xi = 10 * fact; // number of grid points in xi dir
const int num_points_eta = 10 * fact;// number of grid points in eta dir
const int num_points_x = 20 * fact;
const int num_points_total = 50 * fact;

const int num_iterations = 1000;
const int num_iterationss = 1000;
const double tolerance = 1e-5;

std::vector<double> linspace(double start, double end, int num)
{
    std::vector<double> result;
    if (num == 0) {
        return result;
    }
    
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result.push_back(start + i * step);
    }
    
    return result;
}

void write_to_csv(const std::string& filename, const std::vector<double>& data) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    for (size_t i = 0; i < data.size(); ++i) {
        file << data[i];

        // Add comma if it's not the last element
        if (i != data.size() - 1) {
            file << ",";
        }
    }

    file.close();
}

int main(){
// Define xi and eta values
    double ad = 0.1 / fact; // adjustment for bump 
   // std::vector<double> xi(num_points_xi);// define space for vectors 
    std::vector<double> eta1(num_points_x);
    std::vector<double> eta2(num_points_x);
    std::vector<double> xi = linspace((2+ad),(3-ad), num_points_xi);

    for (int i = 0; i < num_points_xi; ++i) {
        eta1[i] = 1.0 - 0.1 * std::sin((xi[i] - 2) * M_PI);
        eta2[i] = 0.1 * std::sin((xi[i] - 2) * M_PI);
    }

    std::vector<double> x = linspace(0,2, num_points_x);
    std::vector<double> x2 = linspace(3,5, num_points_x);


    // Define boundary conditions
    std::vector<double> etax(x); // set up vector 
    std::vector<double> etax1(x);
    std::vector<double> etax2((x2));
    std::vector<double> etax3((x2));

    for (int i = 0; i < num_points_xi; ++i) {
        etax[i] = 0; // fill with ones and zeros 
        etax1[i] = 1;
    }

    std::vector<double> etaT(num_points_total);
    std::vector<double> etaB(num_points_total);
    // adds together two vectors 
    for (int i = 0; i < num_points_total; ++i) {
        etaT[i] = (i < num_points_x) ? etax1[i] : eta1[2*num_points_x];
        etaB[i] = (i < num_points_x) ? etax[i] : eta2[i - num_points_x];
    }

    std::vector<double> xcombine;
    xcombine.insert(xcombine.end(), x.begin(), x.end());
    xcombine.insert(xcombine.end(), xi.begin(), xi.end());
    xcombine.insert(xcombine.end(), x2.begin(), x2.end());

    //std::vector<std::vector<double> > XI(num_points_eta, std::vector<double>(num_points_total + num_points_x));
    std::vector<double> XI[num_points_xi] = {xcombine, xcombine, xcombine,xcombine, xcombine, xcombine,xcombine, xcombine, xcombine,xcombine};
    std::vector<double> xi_old[num_points_xi];

    copy(begin(XI), end(XI), begin(xi_old));

    std::vector<double> yi(num_points_eta);
  

    // Initialize potential grid
    //std::vector<std::vector<double> > V(num_points_eta, std::vector<double>(num_points_total + num_points_x));
    std::vector<double> V[num_points_xi] = {etaB, xcombine, xcombine,xcombine, xcombine, xcombine,xcombine, xcombine, xcombine,etaT};
    //std::vector<double> 
    // Apply boundary conditions
    /*for (int i = 0; i < num_points_eta; ++i) {
        V[i][0] = etaB[i];
        V[i][num_points_eta] = etaT[i];
        //V[i][0] = 0;
        //V[i][num_points_total + num_points_x - 1] = 1;
        //V[i][0] = i * 1.0 / (num_points_eta - 1);
        V[i][num_points_total + num_points_x - 1] = i * 1.0 / (num_points_eta - 1);
    }*/

    // Perform iterative solution of Laplace's equation
    /*for (int iter = 0; iter < num_iterations; ++iter) {
        std::vector<std::vector<double> > V_old = V;

        // Iterate through interior points
        for (int i = 1; i < num_points_xi - 1; ++i) {
            for (int j = 1; j < num_points_total + num_points_x - 1; ++j) {
                // Update potential using Laplace's equation and finite difference 
                V[i][j] = 0.25 * (V_old[i + 1][j] + V_old[i - 1][j] + V_old[i][j + 1] + V_old[i][j - 1]);
            }
        }

        // Check convergence
        double max_diff = 0.0;
        for (int i = 0; i < num_points_eta; ++i) {
            for (int j = 0; j < num_points_total + num_points_x; ++j) {
                max_diff = std::max(max_diff, std::abs(V[i][j] - V_old[i][j]));
            }
        }
        if (max_diff < tolerance) {
            std::cout << "Convergence achieved after " << iter << " iterations.\n";
            break;
        }*/
    //}

 // Iterate through interior points
   for (int iterr = 0; iterr < num_iterationss; ++iterr) {
        for (int i = 1; i < num_points_xi - 1; ++i) {
            for (int j = 1; j < num_points_total + num_points_x - 1; ++j) {
                // Update potential using Laplace's equation and finite difference 
                XI[i][j] = 0.25 * (xi_old[i + 1][j] + xi_old[i - 1][j] + xi_old[i][j + 1] + xi_old[i][j - 1]);
            }
        }

        // Check convergence
        double max_diff = 0.0;
        for (int i = 0; i < num_points_eta; ++i) {
            for (int j = 0; j < num_points_total + num_points_x; ++j) {
                max_diff = std::max(max_diff, std::abs(XI[i][j] - xi_old[i][j]));
            }
        }
        if (max_diff < tolerance) {
            std::cout << "Convergence achieved after " << iterr << " iterations.\n";
            break;
        }
    }

        // Visualize the potential distribution
        // Open a file for writing
        std::ofstream outFile("output.txt");
 
        // Check if the file was opened successfully
       // if (!outFile)
       // {
        //    std::cerr << "Unable to open file for writing." << std::endl;
       //    return 1;
       // }
    std::vector<double> YI = linspace(0,5, num_points_total);

    return 0;
}


/* //Project 2 
// Euler equations. The code should use a finite volume method with the JST scheme.

const double machNum = 0.8; // Mach number 
const double c = 343; // m/s speed of sound 

double f[num_points_xi][num_points_eta]; // inital flux x dir 
double g[num_points_xi][num_points_eta]; // intial flux y dir 
double dt = 1.0; // time step 
double result[num_points_xi][num_points_eta];
std::vector<std::vector<double> > A(num_points_xi, std::vector<double>(num_points_eta));
std::vector<std::vector<double> > R(num_points_xi, std::vector<double>(num_points_eta));

void fvd(){
// Finite Volume Discretization 

 // Perform iterative solution of Aij 
  
    for (int iterrr = 0; iterrr < num_iterations; ++iterrr) {
       // std::vector<std::vector<double> > A_old = A;

        // Iterate through interior points
        for (int i = 1; i < num_points_xi - 1; ++i) {
            for (int j = 1; j < num_points_total + num_points_x - 1; ++j) {
                // Equation 4.3.7
                A[i][j] = 0.5 * (((XI[i+.5][j+.5] - XI[i - .5][j-.5])* (YI[i-.5][j + .5] - YI[i+.5][j-.5]))- ((YI[i+.5][j-.5] - YI[i-.5][j-.5])*(XI[i-.5][j+.5] - XI[i+.5][j-.5]))) ;
                R[i][j] = (0.5 (f[i][j]+ f[i+1][j])*deltaY[i+.5][j]) - (.5(g[i][j]- g[i+1][j])*deltaX[i+.5][j])+(0.5 (f[i][j]+ f[i][j+1])*deltaY[i][j+.5]) - (.5(g[i][j]- g[i][j+1])*deltaX[i][j+.5])+(0.5 (f[i][j]+ f[i-1][j])*deltaY[i-.5][j]) - (.5(g[i][j]- g[i-1][j])*deltaX[i+.5][j])+(0.5 (f[i][j]+ f[i][j-1])*deltaY[i][j-.5]) - (.5(g[i][j]- g[i][j-1])*deltaX[i][j+.5]);
            }
        }
    
    }
} 
//int main() {
 //return 0 ;
//} */



