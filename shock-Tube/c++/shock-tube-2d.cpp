#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>

// Constant and Parameters



const int Nx_l  = 50; const int Nx_r  = 50;
const int Ny_l  = 10; const int Ny_r  = 10;
const int Nx = Nx_l + Nx_r;
const int Ny = Ny_l + Ny_r;

const int N = Nx * Ny;

const double kappa  = 2.0;
const double nu     = 1.4;
const double gamma1 = 1.4;
const double Gamma  = 1.4;

struct Particle {

    // Propierties
    Eigen::Vector2d r = Eigen::Vector2d(0.0,0.0);
    Eigen::Vector2d v = Eigen::Vector2d(0.0,0.0);
    double rho = 0.0, e = 0.0 , P = 0.0;
    double mass = 0.001875;
    
    // Derivatives
    Eigen::Vector2d d_r  = Eigen::Vector2d(0.0,0.0);
    Eigen::Vector2d d_v  = Eigen::Vector2d(0.0,0.0);
    double d_e = 0;
    double d_rho = 0;
    double d_P = 0;
};

std::vector<Particle> Mesh(double x1, double x2, double y1, double y2) {
    // Calculate steps

    double xdim = x2 - x1;
    double ydim = y2 - y1;
    double xstep_l = (xdim/2.0) / Nx_l; double xstep_r = (xdim/2.0) / Nx_r;
    double ystep = ydim / (Ny-1); 

    // Create a mesh of particles
    std::vector<Particle> mesh;
    mesh.reserve(N);

    for(int jj = 1; jj <= Ny ; ++jj){
        for(int ii = 1; ii <= Nx ; ++ii){
    
            Particle p;

            // Left side
            if (ii < Nx_l ){
                p.r = Eigen::Vector2d(x1 + ii*xstep_l, y1 + jj*ystep);
                p.v = Eigen::Vector2d(0.0,0.0);
                p.rho = 1;
                p.e = 2.5;
                p.P = 1;
            }
            // Right side
            else {
                p.r = Eigen::Vector2d( x2 - (Nx - ii) * xstep_r, y1 + jj*ystep);
                p.v = Eigen::Vector2d(0.0,0.0);
                p.rho = 0.2;
                p.e = 1.795;
                p.P = 0.1795;
            }

            mesh.push_back(p);
        }
    }

    return mesh;
}


void Boundary_Reflective(std::vector<Particle>& mesh, double x1, double x2, double y1, double y2) {
    for (auto& p : mesh) {
        // Reflect in x
        if (p.r.x() < x1) {
            p.r.x() = x1 + (x1 - p.r.x());
            p.v.x() *= -1; // Inviert velocity
        } else if (p.r.x() > x2) {
            p.r.x() = x2 - (p.r.x() - x2);
            p.v.x() *= -1;
        }

        // Reflect in y
        if (p.r.y() < y1) {
            p.r.y() = y1 + (y1 - p.r.y());
            p.v.y() *= -1;
        } else if (p.r.y() > y2) {
            p.r.y() = y2 - (p.r.y() - y2);
            p.v.y() *= -1;
        }
    }
}



void Boundary_Periodic(std::vector<Particle>& mesh, double x2, double y2){
    for (auto& p : mesh){

        double posx = std::fmod(p.r.x(), x2);
        
        if (posx < 0) {
            //std::cout<<"pos_ori: "<<posx<<std::endl;
            p.r.x() = posx + x2;
            //std::cout<<"pos-corect: "<< p.r <<std::endl;
        }
        else {
            p.r.x() = posx;
        }

        double posy = std::fmod(p.r.y(), y2);
        
        if (posy < 0) {
            //std::cout<<"pos_ori: "<<pos<<std::endl;
            p.r.y() = posy + y2;
            //std::cout<<"pos-corect: "<< p.r.y() <<std::endl;
        }
        else {
            p.r.y() = posy;
        }
    }
}


double h_len(double mass, double rho) {
    return nu * mass / rho;
}


double Kernel(const Eigen::Vector2d& r_ij, double h)
{
    double ad = 1.0 / (M_PI * h*h*h);

    double r  = r_ij.norm();

    double q = 0.0; 

    if ( (r >= 0.0)  && (r < 1) ) {
        q = ad * (1.0 - 1.5*r*r + 0.75 * r*r*r);
    }
    else if ( (r >= 1) && (r <= 2)){
        q = ad * (2 - r)*(2 - r)*(2 - r) / 4.0;
    }

    return q;
}


Eigen::Vector2d GradKernel2D(const Eigen::Vector2d& r_ij, double h) 
{
   double ad = 10.0 / (7 * M_PI * h * h);
   double r  = r_ij.norm() / h;
   Eigen::Vector2d gradW;
    
    if ( (r >= 0) && (r < 1) ){
        gradW.x() = ad * (-3.0 + 2.25 * r) / (h * h) * r_ij.x();
        gradW.y() = ad * (-3.0 + 2.25 * r) / (h * h) * r_ij.y();
    }
    else if ( (r >= 1) && (r < 2) ) {
        gradW.x() = -ad * (0.75 * ((2 - r)*(2 - r))) / (h * r) *  r_ij.x() ;
        gradW.y() = -ad * (0.75 * ((2 - r)*(2 - r))) / (h * r) *  r_ij.y() ;
    } 
    return gradW;
}


double A_Viscosity(const Eigen::Vector2d& r_i, const Eigen::Vector2d& r_j, const Eigen::Vector2d& v_i, const Eigen::Vector2d& v_j,
                            double rho_i, double  rho_j, double e_i, double e_j, double  h_i, double h_j)
{   
    double alpha  = 1.0, beta = 1.0; 
    double c_i    = std::sqrt((Gamma - 1.0) * e_i);
    double c_j    = std::sqrt((Gamma - 1.0) * e_j);
    Eigen::Vector2d r_ij = r_i - r_j;
    Eigen::Vector2d v_ij = v_i - v_j;
    double c_ij   = (c_i + c_j) * 0.5;
    double rho_ij = (rho_i + rho_j) * 0.5;
    double h_ij   = (h_i + h_j) * 0.5;
    double eta    = 0.1 * h_ij;
    double mu_ij = ( h_ij * v_ij.dot(r_ij)) / (r_ij.dot(r_ij) + eta*eta); 
    
    double visc = 0.0;
    if( r_ij.dot(v_ij) < 0){
        visc = (-alpha * c_ij * mu_ij + beta * mu_ij*mu_ij ) / rho_ij;
    }

    return visc;

}


double Pressure (double rho, double e) 
{
    return (Gamma-1) * rho * e;
}


Eigen::Vector2d Acceleration(double mass,double rho_i, double rho_j, double P_i, double P_j,
                             const Eigen::Vector2d& gradW, double Visc)
{
    Eigen::Vector2d Acc;
    Acc.x() = -mass * ( P_i/(rho_i*rho_i) + P_j/(rho_j*rho_j) + Visc) * gradW.x();
    Acc.y() = -mass * ( P_i/(rho_i*rho_i) + P_j/(rho_j*rho_j) + Visc) * gradW.y();

    return  Acc;
}


double Energy(double mass, double rho_i, double rho_j, double P_i,double  P_j, const Eigen::Vector2d& v_i,
                        const Eigen::Vector2d& v_j, const Eigen::Vector2d& gradW, double Visc) 
{
        return 0.5 * mass * ( P_i/(rho_i*rho_i) + P_j/(rho_j*rho_j) + Visc) * (v_i - v_j).dot(gradW);
}


void nearest_neight(const std::vector<Particle>& mesh,
                    std::vector<int>& pair_i, std::vector<int>& pair_j,
                    std::vector<double>& q, std::vector<Eigen::Vector2d>& dq) {

    double h_ij;

    for (int ii = 0; ii < N; ++ii) {
        for (int jj =ii + 1 ; jj < N; ++jj){

            // relative distance
            Eigen::Vector2d r_ij = mesh[ii].r - mesh[jj].r;

            h_ij = 0.5 * ( h_len(mesh[ii].mass, mesh[ii].rho) + h_len(mesh[jj].mass, mesh[jj].rho) );
        
            // Condition of nearest
            if (r_ij.norm() <= kappa * h_ij){
                pair_i.push_back(ii);
                pair_j.push_back(jj);
                q.push_back( Kernel(r_ij, h_ij) );
                dq.push_back( GradKernel2D(r_ij, h_ij) );

            }
        }       
    }
}


std::vector<Particle> System(std::vector<Particle>& mesh1) {


    std::vector<Particle> mesh = mesh1;
    std::vector<int> pair_i;
    std::vector<int> pair_j;
    std::vector<double> q;
    std::vector<Eigen::Vector2d> dq;

    // Nearest Neighbors
    nearest_neight(mesh, pair_i, pair_j, q, dq);

    int NPairs = pair_i.size();

    // Update self Density
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
        //std::cout << "Pressure " << p.P << "\n";
    }

    // Calculate the System Equations 
    for (int kk = 0; kk < NPairs; ++kk) {
        int pi = pair_i[kk];
        int pj = pair_j[kk];

        double a_visc = A_Viscosity(mesh[pi].r, mesh[pj].r, mesh[pi].v, mesh[pj].v,
                                    mesh[pi].rho, mesh[pj].rho, mesh[pi].e, mesh[pj].e, h_len(mesh[pi].mass, mesh[pi].rho), h_len(mesh[pj].mass, mesh[pj].rho));
        
        mesh[pi].d_v += Acceleration( mesh[pj].mass, mesh[pi].rho, mesh[pj].rho, mesh[pi].P, mesh[pj].P, dq[kk], a_visc);
        mesh[pj].d_v -= Acceleration( mesh[pi].mass, mesh[pj].rho, mesh[pi].rho, mesh[pj].P, mesh[pi].P, dq[kk], a_visc);

        mesh[pi].d_e += Energy(mesh[pj].mass, mesh[pi].rho, mesh[pj].rho, mesh[pi].P, mesh[pj].P, mesh[pi].v, mesh[pj].v, dq[kk], a_visc);
        mesh[pj].d_e -= Energy(mesh[pi].mass, mesh[pj].rho, mesh[pi].rho, mesh[pj].P, mesh[pi].P, mesh[pj].v, mesh[pi].v, dq[kk], a_visc);

    }

    for (Particle& p : mesh) {
        p.d_r = p.v; 
    }

    return mesh;
}


std::vector<Particle> Integration(std::vector<Particle>& mesh, double x1, double x2) {

    double tstep = 0.0005;
    const double tmax = tstep * 400;
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
        std::cout<<"i: "<< step <<std::endl;
    }
    

    return mesh;
}


void csv(const std::vector<Particle>& mesh){
    std::ofstream file("mesh-2d.csv");

    // Escribe los encabezados en el archivo CSV
    file << "x,y,vx,vy,rho,e,P\n";
    
    for (const auto& p : mesh) {
        file << p.r(0) << "," << p.r(1) << "," << p.v(0) << "," << p.v(1) << "," 
        << p.rho << "," << p.e << "," << p.P << "\n";
    }

    file.close();
}


int main(){

    // Define the dimentions of tube
    double t1 = 1.0;
    double x1 = 0.0; double x2 = 1.2;
    double y1 = 0.0; double y2 = 0.3;
    std::vector<Particle> mesh = Mesh(x1,x2,y1,y2);
    std::cout<<" Mesh created" << std::endl;
    std::vector<Particle> int_mesh = Integration(mesh, x1, x2);

    csv(int_mesh);

    return 0;
}
