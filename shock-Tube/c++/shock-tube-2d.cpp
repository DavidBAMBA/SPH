#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>

// Constant and Parameters

const int Nx_l = 10; const int Nx_r = 10;
const int Ny_l = 10; const int Ny_r = 10;
const int Nx = Nx_l + Nx_r;
const int Ny = Ny_l + Ny_r;

const int N = Nx * Ny;

const double kappa  = 2000000000.0;
const double nu     = 1.4;
const double gamma1 = 5.0 / 3.0;
const double Gamma  = 1.4;

struct Particle {

    // Propierties
    Eigen::Vector2d r = Eigen::Vector2d(0.0,0.0);
    Eigen::Vector2d v = Eigen::Vector2d(0.0,0.0);
    double rho = 0.0, e = 0.0 , P = 0.0;
    double mass = 20.0/N;
    
    // Derivatives
    Eigen::Vector2d d_r  = Eigen::Vector2d(0.0,0.0);
    Eigen::Vector2d d_v  = Eigen::Vector2d(0.0,0.0);
    double d_e = 0.0;
    double d_rho = 0.0;
    double d_P = 0.0;
};

std::vector<Particle> Mesh(double x1, double x2, double y1, double y2) {
    // Calculate dimensions
    double xdim = x2 - x1;
    double ydim = y2 - y1;

    // Calculate steps for left and right sides
    double xstep_l = (xdim / 2.0) / Nx_l;
    double xstep_r = (xdim / 2.0) / Nx_r;
    double ystep = ydim / ((Ny - 1) * sqrt(3) / 2);

    // Create a mesh of particles
    std::vector<Particle> mesh;
    mesh.reserve(Nx_l * Ny + Nx_r * Ny);

    // Offsets for hexagonal close packing
    double xOffset_l = xstep_l / 2;
    double xOffset_r = xstep_r / 2;
    double yOffset = ystep * sqrt(3) / 2;

    for (int jj = 1; jj <= Ny; ++jj) {
        for (int ii = 1; ii <= (Nx_l + Nx_r); ++ii) {
            Particle p;

            // Calculate position with offsets for hexagonal packing
            double x, y;
            y = y1 + (jj - 1) * yOffset;

            if (ii <= Nx_l) { // Left side
                x = x1 + (ii - 1) * xstep_l + (jj % 2) * xOffset_l;
                p.r = Eigen::Vector2d(x, y);
                p.v = Eigen::Vector2d(0.0, 0.0);
                p.rho = 1.0;
                p.e = 2.5;
                p.P = 1.0;
            } else { // Right side
                int ii_r = ii - Nx_l; // Adjusted index for right side
                x = x1 + xdim / 2 + (ii_r - 1) * xstep_r + (jj % 2) * xOffset_r;
                p.r = Eigen::Vector2d(x, y);
                p.v = Eigen::Vector2d(0.0, 0.0);
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



void Boundary_Periodic(std::vector<Particle>& mesh, double x1, double x2, double y1, double y2) {
    for (auto& p : mesh) {
        // Ajustar la posición en el eje x
        if (p.r.x() < x1) {
            p.r.x() += (x2 - x1);
        } else if (p.r.x() >= x2) {
            p.r.x() -= (x2 - x1);
        }

        // Ajustar la posición en el eje y
        if (p.r.y() < y1) {
            p.r.y() += (y2 - y1);
        } else if (p.r.y() >= y2) {
            p.r.y() -= (y2 - y1);
        }
    }
}


double h_len(double mass, double rho) {
    return 0.99; //nu * mass / rho;
}



double M6Kernel(const Eigen::Vector2d& r_ij, double h) {
    const double sigma = 7.0/ (478.0 * M_PI * h * h);
    double r = r_ij.norm();
    double q = r / h;
    double w = 0.0;

    if (r > 0.0){

        if (q >= 0.0 && q < 1.0) {
        w = sigma * ((3 - q) * (3 - q) * (3 - q) * (3 - q) * (3 - q) - 
                     6 * (2 - q) * (2 - q) * (2 - q) * (2 - q) * (2 - q) +
                     15 * (1 - q) * (1 - q) * (1 - q) * (1 - q) * (1 - q));
        } else if (q >= 1.0 && q < 2.0) {
            w = sigma * ((3 - q) * (3 - q) * (3 - q) * (3 - q) * (3 - q) - 
                     6 * (2 - q) * (2 - q) * (2 - q) * (2 - q) * (2 - q));
        } else if (q >= 2.0 && q < 3.0) {
            w = sigma * ((3 - q) * (3 - q) * (3 - q) * (3 - q) * (3 - q));
        }

    }

    
    return w;
}



Eigen::Vector2d GradM6Kernel(const Eigen::Vector2d& r_ij, double h) {
    const double sigma = 7.0/ (478.0 * M_PI * h * h);
    double r = r_ij.norm();
    //std::cout<<"r: "<<r<<std::endl;
    double q = r / h;
    //std::cout<<"q: "<<q<<std::endl;
    Eigen::Vector2d gradW = Eigen::Vector2d::Zero();
    Eigen::Vector2d r_hat = r_ij.normalized() / h;

    if (r > 0.0){

        if (q >= 0.0 && q < 1.0) {
            gradW = (sigma * 5 * ((3 - q) * (3 - q) * (3 - q) * (3 - q) - 
                6  * (2 - q) * (2 - q) * (2 - q) * (2 - q) +
                15 * (1 - q) * (1 - q) * (1 - q) * (1 - q))) * r_hat;
        } else if (q >= 1.0 && q < 2.0) {
            gradW = (sigma * ((3 - q) * (3 - q) * (3 - q) * (3 - q) * (3 - q) - 
                     6 * (2 - q) * (2 - q) * (2 - q) * (2 - q) * (2 - q))) * r_hat;
        } else if (q >= 2.0 && q < 3.0) {
            gradW = (sigma * ((3 - q) * (3 - q) * (3 - q) * (3 - q) * (3 - q))) * r_hat;
        }
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
    return (gamma1-1) * rho * e;
}


Eigen::Vector2d Acceleration(double mass,double rho_i, double rho_j, double P_i, double P_j,
                             const Eigen::Vector2d& gradW, double Visc)
{
    Eigen::Vector2d Acc;
    Acc = - mass * ( P_i/(rho_i*rho_i) + P_j/(rho_j*rho_j)) * gradW;
    //std::cout<<"\n acc: "<< Acc.transpose() << std::endl;

    return  Acc;
}


double Energy(double mass, double rho_i, double rho_j, double P_i,double  P_j, const Eigen::Vector2d& v_i,
                        const Eigen::Vector2d& v_j, const Eigen::Vector2d& gradW, double Visc) 
{
        return  mass * ( P_i*v_j/(rho_i*rho_i) + P_j*v_i/(rho_j*rho_j)).dot(gradW);
}


void nearest_neight(const std::vector<Particle>& mesh,
                    std::vector<int>& pair_i, std::vector<int>& pair_j,
                    std::vector<double>& q, std::vector<Eigen::Vector2d>& dq) {

    double h_ij;

    for (int ii = 0; ii < N; ++ii) {
        for (int jj = 0 ; jj < N; ++jj){

            // relative distance
            Eigen::Vector2d r_ij = mesh[ii].r - mesh[jj].r;

            h_ij = 0.5 * ( h_len(mesh[ii].mass, mesh[ii].rho) + h_len(mesh[jj].mass, mesh[jj].rho) );
        
            // Condition of nearest
            if (r_ij.norm() <= kappa * h_ij){
                pair_i.push_back(ii);
                pair_j.push_back(jj);
                q.push_back( M6Kernel(r_ij, h_ij) );
                dq.push_back( GradM6Kernel(r_ij, h_ij) );

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
    //for (Particle& p : mesh) {
        //p.rho = p.mass * (2.0 / (3.0 * h_len(p.mass, p.rho))); 
        //std::cout << "Density: " << p.rho << "\n";
    //}

    // Update Density
    for (int kk = 0; kk < NPairs; ++kk) {

        int pi = pair_i[kk];
        int pj = pair_j[kk];
        
        mesh[pi].rho += mesh[pj].mass * q[kk];
        //mesh[pj].rho += mesh[pi].mass * q[kk];
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
        //mesh[pj].d_v -= Acceleration( mesh[pi].mass, mesh[pj].rho, mesh[pi].rho, mesh[pj].P, mesh[pi].P, dq[kk], a_visc);

        mesh[pi].d_e += Energy(mesh[pj].mass, mesh[pi].rho, mesh[pj].rho, mesh[pi].P, mesh[pj].P, mesh[pi].v, mesh[pj].v, dq[kk], a_visc);
        //mesh[pj].d_e -= Energy(mesh[pi].mass, mesh[pj].rho, mesh[pi].rho, mesh[pj].P, mesh[pi].P, mesh[pj].v, mesh[pi].v, dq[kk], a_visc);

    }

    for (Particle& p : mesh) {
        p.d_r = p.v; 
    }

    return mesh;
}


std::vector<Particle> Integration(std::vector<Particle>& mesh, double x1, double x2, double y1, double y2) {

    double tstep = 0.0005;
    const double tmax = tstep * 400;
    const int    NSteps = static_cast<int>((tmax - tstep) / tstep);

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
        Boundary_Periodic(int_mesh, x1, x2, y1, y2);
        //Boundary_Reflective(int_mesh, x1, x2);
        k2 = System(int_mesh);

        // k3
        for (size_t ii = 0; ii < int_mesh.size(); ++ii) {
            int_mesh[ii].r = k2[ii].r + 0.5 * tstep * k2[ii].d_r;
            int_mesh[ii].v = k2[ii].v + 0.5 * tstep * k2[ii].d_v;
            int_mesh[ii].e = k2[ii].e + 0.5 * tstep * k2[ii].d_e;
        }
        Boundary_Periodic(int_mesh, x1, x2, y1, y2);
        //Boundary_Reflective(int_mesh, x1, x2);
        k3 = System(int_mesh);

        // k4
        for (size_t ii = 0; ii < int_mesh.size(); ++ii) {
            int_mesh[ii].r = k3[ii].r + tstep * k3[ii].d_r;
            int_mesh[ii].v = k3[ii].v + tstep * k3[ii].d_v;
            int_mesh[ii].e = k3[ii].e + tstep * k3[ii].d_e;
        }
        Boundary_Periodic(int_mesh, x1, x2, y1, y2);
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
        Boundary_Periodic(mesh, x1, x2, y1, y2);
        //Boundary_Reflective(mesh, x1, x2);
        //Boundary_Periodic(mesh, x2);

        std::cout <<"  bi: "<< ii <<std::endl; 
        
        }

    return mesh;

}


std::vector<Particle> EulerIntegration(std::vector<Particle>& mesh) {
    double tstep = 0.005;
    int NSteps = 100;
    
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



void InitializeDensity(std::vector<Particle>& mesh) {
    for (int ii = 0; ii < mesh.size(); ++ii) {
        auto& pi = mesh[ii];
        for (const auto& pj : mesh) {
            Eigen::Vector2d r_ij = pi.r - pj.r; 
            double h = h_len(pj.mass, pj.rho); 
            pi.rho += pj.mass * M6Kernel(r_ij, h); 
        }
    }
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
    double x1 = 0.0; double x2 = 0.5;
    double y_1 = 0.0; double y2 = 0.3;
    std::vector<Particle> mesh = Mesh(x1,x2,y_1,y2);
    std::cout<<" Mesh created" << std::endl;
    //InitializeDensity(mesh);
    std::vector<Particle> int_mesh = Integration(mesh, x1, x2, y_1, y2);

    csv(int_mesh);

    return 0;
}
