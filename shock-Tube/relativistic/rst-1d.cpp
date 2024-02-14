#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>

// Constant and Parameters

double x1 = -0.6, x2 = 0.6;
const int Nx_l  = 820; const int Nx_r  = 280;

const int Nx = Nx_l + Nx_r;
const int N = Nx;

const double kappa   = 2.0;
const double nu      = 1.4;
const double gamma1  = 4.0/3.0;
const double G   = 1.4;

 
struct Particle {

    // Propierties
    double r = 0.0;
    double v = 0.0;
    double rho = 0.0, e = 0.0, P = 0.0;
    double mass = 0.001875;

    // Derivatives
    double d_r = 0.0;
    double d_v = 0.0;
    double d_e = 0.0;
    double d_rho = 0.0;
    double d_P = 0.0;

};


std::vector<Particle> Mesh(double x1, double x2) {
    // Calculate steps
    double xdim = x2 - x1;
    double xstep_l = (xdim/2.0) / Nx_l; 
    double xstep_r = (xdim/2.0) / Nx_r;

    // Create a mesh of particles
    std::vector<Particle> mesh;
    mesh.reserve(N);

    for(int ii = 1; ii <= Nx ; ++ii){
    
        Particle p;

        // Left side
        if (ii < Nx_l ){
            p.r = x1 + ii * xstep_l;
            p.v = 0.0;
            p.rho = 1.0;
            p.e = 2.5;
            p.P = 1.0;
        }
        // Right side
        else {
            p.r =  x2 - (Nx - ii) * xstep_r;
            p.v = 0.0;
            p.rho = 0.125;
            p.e = 2.0;
            p.P = 0.125;
        }

    mesh.push_back(p);
    
    }
    return mesh;
}


void Boundary_Periodic(std::vector<Particle>& mesh, double x2){
    for (auto& p : mesh){
        double pos = std::fmod(p.r, x2);
        
        if (pos < 0) {
            std::cout<<"pos_ori: "<<pos<<std::endl;
            p.r = pos + x2;
            std::cout<<"pos-corect: "<< p.r <<std::endl;
        }
        else {
            p.r = pos;
        }
    }
}

void Boundary_Reflective(std::vector<Particle>& mesh, double x1, double x2) {
    for (auto& p : mesh) {
        // Reflect in x
        if (p.r < x1) {
            p.r = x1 + (x1 - p.r);
            p.v *= -1; // Inviert velocity
        } else if (p.r > x2) {
            p.r = x2 - (p.r - x2);
            p.v *= -1;
        }
    }
}  


void Rigid_boundary(std::vector<Particle>& mesh, double x1, double x2) {
    for (auto& p: mesh) {
        if (p.r <= x1 || p.r >= x2) {
            p.v = 0; // Poner la velocidad a 0 en las fronteras
        }
    }
}

double h_len(double mass, double rho) {
    double in = 0.0;
    return  nu * (mass / rho);
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
    double c_i    = std::sqrt((G - 1.0) * e_i);
    double c_j    = std::sqrt((G - 1.0) * e_j);
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

double energy_disipation(double  x_i,double  x_j, double  v_i, double  v_j,
                    double rho_i, double  rho_j,double e_i, double  e_j, double  h_i, double h_j){

    double alpha  = 1.0; double beta = 2.0;
    double c_i = std::sqrt((G - 1.0) * e_i);
    double c_j = std::sqrt((G - 1.0) * e_j);
    double x_ij = x_i - x_j;
    double j = x_ij / std::sqrt(x_ij*x_ij);
    double v_ij = v_i - v_j;
    double c_ij   = (c_i + c_j) * 0.5;
    double rho_ij = (rho_i + rho_j) * 0.5;
    double h_ij   = (h_i + h_j) * 0.5;
    double v_sig = 0.0;

    if (v_ij * x_ij <= 0.0 ){
        v_sig = 0.5 * ( c_i + c_j - beta * v_ij*j );
    } 

    double ei = 0.5 * alpha * v_sig * j*v_i*j*v_i + e_i;
    double ej = 0.5 * alpha * v_sig * j*v_j*j*v_j + e_j;

    double u_diss = (ei - ej) / rho_ij * j ;

    return u_diss;

}


double Pressure (double rho, double e) 
{
    return (gamma-1) * rho * e;
}


double Acceleration(double mass,double rho_i, double rho_j, double P_i, double P_j,double dq, double Visc){
    return  - mass * ( P_i/(rho_i*rho_i) + P_j/(rho_j*rho_j) + Visc) * dq;
}


double Energy(double mass, double rho_i, double rho_j, double P_i,double  P_j, double v_i,double  v_j,double dq, double Visc) 
{
        return   mass * ( P_i/(rho_i*rho_i) + P_j/(rho_j*rho_j) ) * (v_i - v_j) * dq;
}


void nearest_neight(const std::vector<Particle>& mesh,
                    std::vector<int>& pair_i, std::vector<int>& pair_j,
                    std::vector<double>& q, std::vector<double>& dq) {

    double h_ij;

    for (int ii = 0; ii < N - 1; ++ii) {
        for (int jj = ii + 1 ; jj < N; ++jj){

            // relative distance
            double r_ij = mesh[ii].r - mesh[jj].r;

            h_ij = 0.5 * ( h_len(mesh[ii].mass, mesh[ii].rho) + h_len(mesh[jj].mass, mesh[jj].rho) );
        
            // Condition of nearest
            if (std::abs(r_ij) <= kappa * h_ij){
                //std::cout << "Pair " << ii << " and " << jj << " with q: " << Kernel(r_ij, h_ij) << " and dq: " << D_Kernel(r_ij, h_ij) << "\n";

                //std::cout<<"ok"<< std::endl;
                pair_i.push_back(ii);
                pair_j.push_back(jj);
                q.push_back( Kernel(r_ij, h_ij) );
                dq.push_back( D_Kernel(r_ij, h_ij) );

            }
        }       
    }
}


std::vector<Particle> System(std::vector<Particle>& mesh1) {

    std::vector<Particle> mesh = mesh1;
    std::vector<int> pair_i;
    std::vector<int> pair_j;
    std::vector<double> q;
    std::vector<double> dq;

    // Nearest Neighbors
    nearest_neight(mesh, pair_i, pair_j, q, dq);

    int NPairs = pair_i.size();

    // Update Self Density
    for (Particle& p : mesh) {
        p.rho = p.mass * (2.0 / (3.0 * h_len(p.mass, p.rho))); 
        //std::cout << "Density: " << p.rho << "\n";
    }

    // Update Density
    for (int kk = 0; kk < NPairs; ++kk) {
        int pi = pair_i[kk];
        int pj = pair_j[kk];
        
        mesh[pi].rho += mesh[pj].mass * q[kk];
        mesh[pj].rho += mesh[pi].mass * q[kk];
        //std::cout << "Density update for pair " << pair_i[kk] << " and " << pair_j[kk] << ": " <<  mesh[pi].rho << " and " << mesh[pj].rho << "\t q: "<< q[kk] <<"\n";
    }
    
    // Update Pressure
    for (Particle& p : mesh) {
        p.P = Pressure(p.rho, p.e);
    }

    // Calculate the System Equations 
    for (int kk = 0; kk < NPairs; ++kk) {
        int pi = pair_i[kk];
        int pj = pair_j[kk];

        double a_visc = A_Viscosity(mesh[pi].r, mesh[pj].r, mesh[pi].v, mesh[pj].v,
                                    mesh[pi].rho, mesh[pj].rho, mesh[pi].e, mesh[pj].e, h_len(mesh[pi].mass, mesh[pi].rho), h_len(mesh[pj].mass, mesh[pj].rho));
        double u_visc = energy_disipation(mesh[pi].r, mesh[pj].r, mesh[pi].v, mesh[pj].v,
                                    mesh[pi].rho, mesh[pj].rho, mesh[pi].e, mesh[pj].e, h_len(mesh[pi].mass, mesh[pi].rho), h_len(mesh[pj].mass, mesh[pj].rho));
        
        mesh[pi].d_v += Acceleration( mesh[pj].mass, mesh[pi].rho, mesh[pj].rho, mesh[pi].P, mesh[pj].P, dq[kk], a_visc);
        mesh[pj].d_v -= Acceleration( mesh[pi].mass, mesh[pj].rho, mesh[pi].rho, mesh[pj].P, mesh[pi].P, dq[kk], a_visc);

        mesh[pi].d_e += Energy(mesh[pj].mass, mesh[pi].rho, mesh[pj].rho, mesh[pi].P, mesh[pj].P, mesh[pi].v, mesh[pj].v, dq[kk], u_visc);
        mesh[pj].d_e -= Energy(mesh[pi].mass, mesh[pj].rho, mesh[pi].rho, mesh[pj].P, mesh[pi].P, mesh[pj].v, mesh[pi].v, dq[kk], u_visc);

    }

    for (Particle& p : mesh) {
        p.d_r = p.v; 
    }

    return mesh;
}


std::vector<Particle> Integration(std::vector<Particle>& mesh, double x1, double x2) {

    double tstep = 0.0005;
    const double tmax = tstep * 200;
    const int    NSteps = static_cast<int>((tmax - tstep) / tstep);
    double t = 0.0005;

    // Loop for the integration
    for (int ii = 0; ii <= NSteps; ++ii) {

        std::vector<Particle> int_mesh = mesh;
        std::vector<Particle> k1;
        std::vector<Particle> k2;
        std::vector<Particle> k3;
        std::vector<Particle> k4;
        
        k1 = System(int_mesh);

        // k2
        for (size_t ii = 0; ii < int_mesh.size(); ++ii) {
            int_mesh[ii].r = k1[ii].r + 0.5 * tstep * k1[ii].d_r;
            int_mesh[ii].v = k1[ii].v + 0.5 * tstep * k1[ii].d_v;
            int_mesh[ii].e = k1[ii].e + 0.5 * tstep * k1[ii].d_e;
        }
        //Boundary_Periodic(int_mesh, x2);
        //Boundary_Reflective(int_mesh, x1, x2);
        k2 = System(int_mesh);

        // k3
        for (size_t ii = 0; ii < int_mesh.size(); ++ii) {
            int_mesh[ii].r = k2[ii].r + 0.5 * tstep * k2[ii].d_r;
            int_mesh[ii].v = k2[ii].v + 0.5 * tstep * k2[ii].d_v;
            int_mesh[ii].e = k2[ii].e + 0.5 * tstep * k2[ii].d_e;
        }
        //Boundary_Periodic(int_mesh, x2);
        //Boundary_Reflective(int_mesh, x1, x2);
        k3 = System(int_mesh);

        // k4
        for (size_t ii = 0; ii < int_mesh.size(); ++ii) {
            int_mesh[ii].r = k3[ii].r + tstep * k3[ii].d_r;
            int_mesh[ii].v = k3[ii].v + tstep * k3[ii].d_v;
            int_mesh[ii].e = k3[ii].e + tstep * k3[ii].d_e;
        }
        //Boundary_Periodic(int_mesh, x2);
        //Boundary_Reflective(int_mesh, x1, x2);
        k4 = System(int_mesh);

        // Update Mesh
        for (size_t ii = 0; ii < mesh.size(); ++ii) {
            mesh[ii].r   += (tstep / 6.0) * (k1[ii].d_r + 2 * k2[ii].d_r + 2 * k3[ii].d_r + k4[ii].d_r);
            mesh[ii].v   += (tstep / 6.0) * (k1[ii].d_v + 2 * k2[ii].d_v + 2 * k3[ii].d_v + k4[ii].d_v);
            mesh[ii].e   += (tstep / 6.0) * (k1[ii].d_e + 2 * k2[ii].d_e + 2 * k3[ii].d_e + k4[ii].d_e);
            mesh[ii].P   = k4[ii].P; 
            mesh[ii].rho = k4[ii].rho; 

        }
        //Boundary_Periodic(mesh, x2);
        //Boundary_Reflective(mesh, x1, x2);
        //Boundary_Periodic(mesh, x2);

        t += tstep;

        std::cout <<"  bi: "<< ii <<std::endl; 
        
        }

    return mesh;

}


std::vector<Particle> EulerIntegration(std::vector<Particle>& mesh) {
    double tstep = 0.000005;
    int NSteps = 30000;
    
    std::vector<Particle> int_mesh = mesh;

    for (int step = 0; step < NSteps; ++step) {

        int_mesh = System(mesh);

        for (size_t ii = 0; ii < int_mesh.size(); ++ii) {
            mesh[ii].r   = int_mesh[ii].r + tstep * int_mesh[ii].d_r;
            mesh[ii].v   = int_mesh[ii].v + tstep * int_mesh[ii].d_v;
            mesh[ii].e   = int_mesh[ii].e + tstep * int_mesh[ii].d_e;
            mesh[ii].P   = int_mesh[ii].P + tstep * int_mesh[ii].d_P;
            mesh[ii].rho = int_mesh[ii].rho + tstep * int_mesh[ii].d_rho;
        }
    }
    

    return mesh;
}


void csv(const std::vector<Particle>& mesh){
    std::ofstream file("mesh-1d.csv");

    // Escribe los encabezados en el archivo CSV
    file << "x,vx,rho,e,P\n";
    
    for (const auto& p : mesh) {
        file << p.r << "," << p.v << "," << p.rho << "," << p.e << "," << p.P << "\n";
    }

    file.close();
}

int main(){

    // Define the dimentions of tube

    std::vector<Particle> mesh = Mesh(x1,x2);
    std::cout<<" Mesh created" << std::endl;
    std::vector<Particle> int_mesh = Integration(mesh, x1, x2);

    csv(int_mesh);

    
    return 0;
}
