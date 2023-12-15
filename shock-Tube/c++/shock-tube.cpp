#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>

// Constant and Parameters
const int    N             = 400;
const double mass_constant = 0.001875;
const double kappa         = 2;
const double nu            = 1.4;
const double gamma1        = 1.4;
const int    NParams       = 5;
const double Gamma         = 1.4;


template <typename Derived1, typename Derived2>
Eigen::VectorXd h_len(const Eigen::ArrayBase<Derived1>& mass, const Eigen::ArrayBase<Derived2>& density, double nu) {

    Eigen::ArrayXd in = Eigen::VectorXd::Zero(mass.size());
    return in + nu * (mass.array() / density.array());
}

double Kernel(double x_ij, double h)
{
    double ad = 1.0 / h;
    double r = std::abs(x_ij) / h;  // Normalizing distance by smoothing length
    
    // Initial values
    double q = 0.0; 

    if ( (r >= 0.0)  && (r < 1) ) {
        q = ad * (2.0/3.0 - r*r + 0.5 * r*r*r);
    }
    else if ( (r >= 1) && (r <= 2)){
        q = ad * (2 - r)*(2 - r)*(2 - r) / 6.0;
    }

    return q;
}

double D_Kernel(double x_ij, double h)
{
    double ad = 1.0 / h;
    double r  = std::abs(x_ij) / h;  // Normalizing distance by smoothing length
    
    // Initial values
    double dq = 0.0;

    // Masks for different conditions
    if ( (r >= 0) && (r < 1) ){
        dq = ad * (-2.0 + 1.5 * r)*(x_ij / ( h*h ));
    }
    else if ( (r >= 1) && (r < 2) ) {
        dq = -ad * (0.5 * ((2 - r)*(2 - r))) * (x_ij / (h * std::abs(x_ij)));
    }

    return dq;
}

double  A_Viscosity(double  x_i,double  x_j, double  v_i, double  v_j,
                    double rho_i, double  rho_j,double e_i, double  e_j, double  h_i, double h_j)
{   
    double alpha  = 1.0, beta = 1.0; 
    double c_i    = std::sqrt((Gamma - 1.0) * e_i);
    double c_j    = std::sqrt((Gamma - 1.0) * e_j);
    double x_ij   = x_i - x_j;
    double v_ij   = v_i - v_j;
    double c_ij   = (c_i + c_j) * 0.5;
    double rho_ij = (rho_i + rho_j) * 0.5;
    double h_ij   = (h_i + h_j) * 0.5;
    double theta  = 0.1 * h_ij;
    double phi_ij = ( h_ij * v_ij * x_ij) / (x_ij*x_ij + theta*theta); 
    
    double visc = 0.0;
    if( x_ij * v_ij < 0){
        visc = (-alpha * c_ij * phi_ij + beta * phi_ij*phi_ij ) / rho_ij;
    }

    return visc;

}

template <typename Derived>
Eigen::VectorXd Pressure (const Eigen::ArrayBase<Derived>& rho, const Eigen::ArrayBase<Derived>& e) 
{
    return (Gamma-1) * rho.array() * e.array();
}

double Acceleration(double mass,double rho_i, double rho_j, double P_i, double P_j,double dq, double Visc)
                        
{
    return  -mass * ( P_i/(rho_i*rho_i) + P_j/(rho_j*rho_j) + Visc) * dq;
}

double Energy(double mass, double rho_i, double rho_j, double P_i,double  P_j, double v_i,double  v_j,double dq, double Visc) 
{
        return 0.5 * mass * ( P_i/(rho_i*rho_i) + P_j/(rho_j*rho_j) + Visc) * (v_i - v_j) * dq;
}

template <typename Derived1, typename Derived2>
std::pair<Eigen::ArrayXd, Eigen::ArrayXd> integrate(double t, const Eigen::ArrayBase<Derived1>& S_flat, const Eigen::ArrayBase<Derived2>& mass) {


    Eigen::MatrixXd S  = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Map(static_cast<const Derived1&>(S_flat).data(), N, NParams);

    std::cout << "Initial S matrix: \n" << S << "\n";

    Eigen::VectorXd h = h_len(mass.array(), S.col(2).array(), nu);

    //std::cout << "Calculated h: \n" << h.transpose() << "\n";

    std::vector<int>    pair_i;
    std::vector<int>    pair_j;
    std::vector<double> q;
    std::vector<double> dq;

    // Nearest Neighbors
    for (int ii = 0; ii < N - 1; ++ii) {
        for (int jj = ii + 1; jj < N; ++jj) {

            double x_ij = S(ii, 0) - S(jj, 0);
            double h_ij = 0.5 * (h(ii) + h(jj));

            if (std::abs(x_ij) <= kappa * h_ij) {
                //std::cout << "Pair " << ii << " and " << jj << " with q: " << Kernel(x_ij, h_ij) << " and dq: " << D_Kernel(x_ij, h_ij) << "\n";

                pair_i.push_back(ii);
                pair_j.push_back(jj);

                q.push_back( Kernel(x_ij, h_ij) );
                dq.push_back( D_Kernel(x_ij, h_ij) );
            }
        }
    }

    int NPairs = pair_i.size();

    //std::cout << "Total pairs: " << pair_i.size() << "\n\n";

    // Update the self Density
    S.col(2) = mass.array() * (2.0 / (3.0 * h.array()));
    //std::cout << "density: " << S.col(2).transpose() << "\n";

    Eigen::MatrixXd dS = Eigen::MatrixXd::Zero(S.rows(), S.cols());

    // Update Density
    for (int kk = 0; kk < NPairs; ++kk) {
        int pi = pair_i[kk];
        int pj = pair_j[kk];

        S(pi, 2) += mass(pj) * q[kk];
        S(pj, 2) += mass(pi) * q[kk];

        //std::cout << "Density update for pair " << pair_i[kk] << " and " << pair_j[kk] << ": " << S(pair_i[kk], 2) << " and " << S(pair_j[kk], 2) << "\n";

    }

    //Update pressure
    S.col(4) = Pressure(S.col(2).array(), S.col(3).array());
    //std::cout << "pressure: " << S.col(4).transpose() << "\n";
    // Calculate the System Equations 
    for (int kk = 0; kk < NPairs; ++kk) {
        int pi = pair_i[kk];
        int pj = pair_j[kk];

        double a_visc = A_Viscosity(S(pi, 0), S(pj, 0), S(pi, 1), S(pj, 1),
                                    S(pi, 2), S(pj, 2), S(pi, 3), S(pj, 3), h(pi), h(pj));

        dS(pi, 1) += Acceleration( mass(pj), S(pi, 2), S(pj, 2), S(pi, 4), S(pj, 4), dq[kk], a_visc);
        dS(pj, 1) -= Acceleration( mass(pi), S(pj, 2), S(pi, 2), S(pj, 4), S(pi, 4), dq[kk], a_visc);

        dS(pi, 3) += Energy(mass(pj), S(pi, 2), S(pj, 2), S(pi, 4), S(pj, 4), S(pi, 1), S(pj, 1), dq[kk], a_visc);
        dS(pj, 3) -= Energy(mass(pi), S(pj, 2), S(pi, 2), S(pj, 4), S(pi, 4), S(pj, 1), S(pi, 1), dq[kk], a_visc);
    
        //std::cout << "Acceleration for pair " << pair_i[kk] << " and " << pair_j[kk] << ": " << dS(pair_i[kk], 1) << " and " << dS(pair_j[kk], 1) << "\n";
        //std::cout << "Energy for pair " << pair_i[kk] << " and " << pair_j[kk] << ": " << dS(pair_i[kk], 3) << " and " << dS(pair_j[kk], 3) << "\n";
    }


    dS.col(2).setZero();
    dS.col(4).setZero();
    dS.col(0) = S.col(1);

    //std::cout << "Final S matrix: \n" << S << "\n";
    Eigen::VectorXd S_flat2 = Eigen::Map<Eigen::VectorXd>(S.data(), N * NParams);

    Eigen::VectorXd dS_flat = Eigen::Map<Eigen::VectorXd>(dS.data(), N * NParams);

    return std::make_pair(dS_flat.array(), S_flat2.array());
}

template <typename Derived1, typename Derived2>
Eigen::MatrixXd Integration(const Eigen::ArrayBase<Derived1>& In_S_flat, const Eigen::ArrayBase<Derived2>& mass) {
    const double tstep = 0.00005;
    const double tmax = tstep * 500;
    const int    NSteps = static_cast<int>((tmax - tstep) / tstep);

    Eigen::ArrayXd S_flat  = In_S_flat;
    

    // MatrixXd for results
    Eigen::MatrixXd S_i(N, NParams); // No need to multiply by NSteps
    Eigen::MatrixXd S_item(N, NParams); 
    Eigen::ArrayXd S_temp(N*NParams), dS_flat(N*NParams);

    double t = 0.00005;

    // Loop for the integration
    for (int ii = 0; ii <= NSteps; ++ii) {
        Eigen::ArrayXd k1(N * NParams);
        Eigen::ArrayXd k2(N * NParams);
        Eigen::ArrayXd k3(N * NParams);
        Eigen::ArrayXd k4(N * NParams);

        std::tie(dS_flat, S_temp) = integrate(t, S_flat, mass);
        // k1
        k1 = tstep * dS_flat;
        
        // k2
        S_flat = S_temp + 0.5 * k1;
        //S_item = Eigen::Map<Eigen::MatrixXd>(S_flat.data(), N, NParams);
        //std::cout <<"k2 : "<< ii <<"\n\n\n S_in: \n\n\n"<< S_item << std::endl; 


        std::tie(dS_flat, S_temp) = integrate(t, S_flat, mass);
        k2 = tstep * dS_flat;
        // k3
        S_flat = S_temp + 0.5 * k2;

        std::tie(dS_flat, S_temp) = integrate(t, S_flat, mass);
        k3 = tstep * dS_flat;
        // k4
        S_flat = S_temp + k3;

        std::tie(dS_flat, S_temp) = integrate(t, S_flat, mass);
        k4 = tstep * dS_flat;
        
        // Update S_flat
        S_flat = S_temp + (1.0/6.0) * (k1 + 2*k2 + 2*k3 + k4);
        
        t += tstep;
        //std::cout << "\t i: " << ii << " W[i]: \n" << S_flat << std::endl; 
        S_i = Eigen::Map<Eigen::MatrixXd>(S_flat.data(), N, NParams);
        std::cout <<"i: "<< ii <<"\n\n\n S_in: \n\n\n"<< S_i << std::endl; 

        }

    S_i = Eigen::Map<Eigen::MatrixXd>(S_flat.data(), N, NParams);
    //std::cout << "\n\n\n S_int: \n\n\n"<< S_i.col(0) << std::endl; 


    return S_i;
}

void save_to_csv(const Eigen::MatrixXd& matrix, const std::string& filename) {
    std::ofstream out(filename);
    
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            out << matrix(i, j);
            if (j < matrix.cols() - 1) 
                out << ",";  // separate values with commas
        }
        out << "\n";  // new line for each row
    }
    out.close();
}

Eigen::VectorXd linspace(double start, double end, int num, bool endpoint = true) {
    Eigen::VectorXd linspaced = Eigen::VectorXd::Zero(endpoint ? num : num - 1);
    double delta = (end - start) / (endpoint ? num - 1 : num);

    for (int i = 0; i < (endpoint ? num : num - 1); ++i) {
        linspaced(i) = start + delta * i;
    }

    return linspaced;
}

int main(){
    // NUmber of particles 
    int p_left  = 200;
    int p_right = 200;
    double t1   = 1.0;

    //Initialize Variables
    Eigen::ArrayXd mass(N);
    mass.setConstant(mass_constant);

    Eigen::VectorXd x(N);
    Eigen::VectorXd x_left  = linspace(-0.6, 0.0, p_left+1, false);
    Eigen::VectorXd x_right = linspace(0.0, 0.6, p_right+1, false);

    x.head(p_left)  = x_left;
    x.tail(p_right) = x_right;
    // Density
    Eigen::VectorXd rho(N);
    rho.head(p_left).setConstant(1.0);
    rho.tail(p_right).setConstant(1.0);
    // Energy
    Eigen::VectorXd e(N);
    e.head(p_left).setConstant(2.5);  
    e.tail(p_right).setConstant(1.75);
    // Velocities
    Eigen::VectorXd v(N);
    v.setConstant(0.0);  
    // Pressures
    Eigen::VectorXd P(N);
    P.head(p_left).setConstant(1.0);
    P.tail(p_right).setConstant(0.0001);

    Eigen::MatrixXd S(N, 5);
    S.col(0) = x;
    S.col(1) = v;
    S.col(2) = rho;
    S.col(3) = e;
    S.col(4) = P;

    Eigen::VectorXd S_flat = Eigen::Map<Eigen::VectorXd>(S.data(), N * 5);

    Eigen::MatrixXd S_int = Integration(S_flat.array(), mass.array());
    
    save_to_csv(S_int, "output.csv");

    //std::cout<<"shape: "<<S_int.size()<<std::endl;
    
    return 0;
}
