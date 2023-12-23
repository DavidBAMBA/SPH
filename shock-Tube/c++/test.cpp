#include <iostream>
#include <Eigen/Dense>
#include <vector>

const int N = 10;       // Mock value for the sake of the test
const int NParams = 5;  // Mock value for the sake of the test

Eigen::MatrixXd integrate(double t, const Eigen::MatrixXd& S, const Eigen::ArrayXd& mass) {
    return Eigen::MatrixXd::Constant(S.rows(), S.cols(), 1.0);
}

template <typename Derived1, typename Derived2>
std::vector<Eigen::MatrixXd> Integration(const Eigen::ArrayBase<Derived1>& In_S_flat, const Eigen::ArrayBase<Derived2>& mass) {
    const double tstep = 0.005;
    const double tmax = tstep * 2000;
    const int    NSteps = static_cast<int>((tmax - tstep) / tstep);

    Eigen::ArrayXd S_flat = In_S_flat;

    // Vector of MatrixXd for results
    std::vector<Eigen::MatrixXd> S_i(NSteps, Eigen::MatrixXd( N, NParams));
    
    double t = 0.005;

    // Loop for the integration
    for (int ii = 0; ii < NSteps; ++ii) {
        Eigen::ArrayXd S_temp;
        Eigen::ArrayXd k1(N * NParams);
        Eigen::ArrayXd k2(N * NParams);
        Eigen::ArrayXd k3(N * NParams);
        Eigen::ArrayXd k4(N * NParams);

        // k1
        k1 = tstep * integrate(t, S_flat.matrix(), mass);
        
        // k2
        S_temp = S_flat + 0.5 * k1;
        k2 = tstep * integrate(t + 0.5 * tstep, S_temp.matrix(), mass);
        
        // k3
        S_temp = S_flat + 0.5 * k2;
        k3 = tstep * integrate(t + 0.5 * tstep, S_temp.matrix(), mass);
        
        // k4
        S_temp = S_flat + k3;
        k4 = tstep * integrate(t + tstep, S_temp.matrix(), mass);
        
        // Update S_flat using RK4
        S_flat += (1.0/6.0) * (k1 + 2*k2 + 2*k3 + k4).array();
        
        // Reshape and store result in vector of matrices
        S_i[ii] = Eigen::Map<Eigen::MatrixXd>(S_flat.data(), N, NParams);

        t += tstep;
        //std::cout << "\t i: " << ii << "\t S: " << S_i[ii] <<"\t"<< std::endl;  // Print just one element of the matrix for simplicity. You can modify this to print more or the whole matrix.
    }

    return S_i;
}

Eigen::MatrixXd analyticalSolution(double t, const Eigen::ArrayXd& initial) {
    Eigen::ArrayXd resultArray = initial + t;
    return Eigen::Map<Eigen::MatrixXd>(resultArray.data(), N, NParams);
}


void test_Integration() {
    Eigen::ArrayXd In_S_flat(N * NParams);
    Eigen::ArrayXd mass(N);

    // Mock initial conditions
    In_S_flat.setConstant(2.0);
    mass.setConstant(1.0);

    // Call the Integration function
        //std::cout << "First matrix from the results:\n" << In_S_flat[0] << std::endl;

    std::vector<Eigen::MatrixXd> results = Integration(In_S_flat, mass);
    std::cout << "First matrix from the results:\n" << results[100] << std::endl;

    double t = 0.005;
    bool allClose = true;
    for (const auto& result : results) {
        Eigen::MatrixXd analytical = analyticalSolution(t, In_S_flat);
        
        // Check if the RK4 result is close to the analytical solution
        if (!result.isApprox(analytical)) {
            allClose = false;
            break;
        }
        
        t += 0.005;
    }
    //std::cout << "First matrix from the results:\n" << analytica[90] << std::endl;

    if (allClose) {
        std::cout << "RK4 results are close to the analytical solution!" << std::endl;
    } else {
        std::cout << "There's a discrepancy between RK4 and the analytical solution." << std::endl;
    }
}
#include <iostream>
#include <Eigen/Dense>

Eigen::VectorXd linspace(double start, double end, int num, bool endpoint = true) {
    Eigen::VectorXd linspaced = Eigen::VectorXd::Zero(endpoint ? num : num - 1);
    double delta = (end - start) / (endpoint ? num - 1 : num);

    for (int i = 0; i < (endpoint ? num : num - 1); ++i) {
        linspaced(i) = start + delta * i;
    }

    return linspaced;
}

int main() {
    // Generar 320 puntos de -0.6 a 0, sin incluir 0
    Eigen::VectorXd xs1 = linspace(-0.6, 0.001, 320+1, false);
    
    // Generar 80 puntos de 0 a 0.6, sin incluir 0.6
    Eigen::VectorXd xs2 = linspace(0, 0.6059, 80+1, false);

    // Imprimir los vectores
    std::cout << std::endl << xs1 << std::endl;
    std::cout << std::endl << xs2 << std::endl;

    return 0;
}


//int main() {
  //  test_Integration();
    //return 0;
//}


