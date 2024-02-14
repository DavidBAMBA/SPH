#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>


// Constant and Parameters

const int Nx_l = 10; const int Nx_r = 5;
const int Ny_l = 5; const int Ny_r = 5;
const int Nx = Nx_l + Nx_r;
const int Ny = Ny_l + Ny_r;

const int N = Nx * Ny;

const double nu     = 1.4;
const double gamma1 = 5.0 / 3.0;
const double Gamma  = 1.4;

struct Particle {
    Eigen::Vector2d r = Eigen::Vector2d(0.0, 0.0); // Position in 2D
    Eigen::Vector2d v = Eigen::Vector2d(0.0, 0.0); // Velocity in 2D
    double rho = 0.0, e = 0.0, P = 0.0; // Density, internal energy, pressure
    double mass = 2.0/N;
    // Derivatives
    Eigen::Vector2d d_r = Eigen::Vector2d(0.0, 0.0);
    Eigen::Vector2d d_v = Eigen::Vector2d(0.0, 0.0);
    double d_e = 0.0;
    double d_rho = 0.0;
    double d_P = 0.0;
    bool ghost = false; // Indicates if the particle is a ghost particle
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
    mesh.reserve(N);

    // Offsets for hexagonal close packing
    double xOffset_l = xstep_l / 2.0;
    double xOffset_r = xstep_r / 2.0;
    double yOffset = ystep * sqrt(3) / 2.0;

    for (int jj = 1; jj <= Ny; ++jj) {
        for (int ii = 0; ii <= (Nx_l + Nx_r) + 1 ; ++ii) {
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
                p.rho = 0.125;
                p.e = 1.795;
                p.P = 0.1;
            }

            mesh.push_back(p);
        }
    }



    return mesh;
}

void UpdateGhostParticles(std::vector<Particle>& particles) {
    for (auto& p : particles) {
        if (p.ghost) {
            // Mantener la posición estática (opcional, depende de si quieres ajustar la posición en algún momento)
            // p.r = Eigen::Vector2d(valor_x_fijo, valor_y_fijo);

            // Restablecer la velocidad a cero para simular una pared sólida
            p.v = Eigen::Vector2d(0.0, 0.0);

            // Restablecer derivadas si es necesario
            p.d_v = Eigen::Vector2d(0.0, 0.0);
        }
    }
}

void AddGhostParticles(std::vector<Particle>& mesh, double x1, double x2, double y1, double y2) {
    // Calculate the step size based on the ghost_layer_width
    // Assuming a uniform distribution of ghost particles along the layer
    double ghost_layer_width_x = 0.05;
    double ghost_layer_width_y = 0.01;
    double mass = 2.0;;
    double x_step = (x2 - x1) / 20; // Example division - adjust based on your simulation needs
    double y_step = (y2 - y1) / 20; // Example division - adjust based on your simulation needs
    
     // Top and bottom ghost layers
    for (double x = x1 - ghost_layer_width_x; x <= x2 + ghost_layer_width_x; x += x_step) {
        // Bottom layer
        Particle bottom_ghost;
        bottom_ghost.r = Eigen::Vector2d(x, y1 - ghost_layer_width_y);
        bottom_ghost.mass = mass;
        bottom_ghost.ghost = true;
        mesh.push_back(bottom_ghost);
        
        // Top layer
        Particle top_ghost;
        top_ghost.r = Eigen::Vector2d(x, y2 + ghost_layer_width_y);
        top_ghost.mass = mass;
        top_ghost.ghost = true;
        mesh.push_back(top_ghost);
    }
    
    // Left and right ghost layers
    for (double y = y1; y <= y2 + 0.01; y += y_step) {
        // Left layer
        Particle left_ghost;
        left_ghost.r = Eigen::Vector2d(x1 - ghost_layer_width_x, y);
        left_ghost.mass = mass;
        left_ghost.ghost = true;
        mesh.push_back(left_ghost);
        
        // Right layer
        Particle right_ghost;
        right_ghost.r = Eigen::Vector2d(x2 + ghost_layer_width_x, y);
        right_ghost.mass = mass;
        right_ghost.ghost = true;
        mesh.push_back(right_ghost);
    }
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
        //if (p.r.x() < x1) {
          //  p.r.x() += (x2 - x1);
        //} else if (p.r.x() >= x2) {
         //   p.r.x() -= (x2 - x1);
        //}

        // Ajustar la posición en el eje y
        if (p.r.y() < y1) {
            p.r.y() += (y2 - y1);
        } else if (p.r.y() >= y2) {
            p.r.y() -= (y2 - y1);
        }
    }
}


double h_len(double mass, double rho) {
    return std::sqrt(mass / rho);
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


double momentum_viscosity(const Eigen::Vector2d& r_i, const Eigen::Vector2d& r_j, const Eigen::Vector2d& v_i, const Eigen::Vector2d& v_j,
                            double rho_i, double  rho_j, double e_i, double e_j, double  h_i, double h_j)
{   
    double alpha  = 1.0; double beta = 2.0;
    double c_i = std::sqrt((Gamma - 1.0) * e_i);
    double c_j = std::sqrt((Gamma - 1.0) * e_j);
    Eigen::Vector2d r_ij = r_i - r_j;
    Eigen::Vector2d j = r_ij.normalized();
    Eigen::Vector2d v_ij = v_i - v_j;
    double c_ij   = (c_i + c_j) * 0.5;
    double rho_ij = (rho_i + rho_j) * 0.5;
    double h_ij   = (h_i + h_j) * 0.5;
    double v_sig = 0.0; 

    if (v_ij.dot(j) <= 0.0 ){
        v_sig = 0.5 * ( c_i + c_j - beta * v_ij.dot(j) );

    }

    double visc = -alpha * v_sig * v_ij.dot(j) / rho_ij;

    return visc;

}

Eigen::Vector2d energy_disipation(const Eigen::Vector2d& r_i, const Eigen::Vector2d& r_j, const Eigen::Vector2d& v_i, const Eigen::Vector2d& v_j,
                                double rho_i, double  rho_j, double e_i, double e_j, double  h_i, double h_j){

    double alpha  = 1.0; double beta = 2.0;
    double c_i = std::sqrt((Gamma - 1.0) * e_i);
    double c_j = std::sqrt((Gamma - 1.0) * e_j);
    Eigen::Vector2d r_ij = r_i - r_j;
    Eigen::Vector2d j = r_ij.normalized();
    Eigen::Vector2d v_ij = v_i - v_j;
    double c_ij   = (c_i + c_j) * 0.5;
    double rho_ij = (rho_i + rho_j) * 0.5;
    double h_ij   = (h_i + h_j) * 0.5;
    double v_sig = 0.0; 
    if (v_ij.dot(j) <= 0.0 ){
        v_sig = 0.5 * ( c_i + c_j - beta * v_ij.dot(j) );
    } 

    double ei = 0.5 * alpha * v_sig * j.dot(v_i)*j.dot(v_i) + e_i;
    double ej = 0.5 * alpha * v_sig * j.dot(v_j)*j.dot(v_j) + e_j;

    Eigen::Vector2d u_diss = (ei -ej) / rho_ij * j ;

    return - u_diss;

}

double Pressure (double rho, double e) 
{
    return (gamma1-1) * rho * e;
}


Eigen::Vector2d Acceleration(double mass,double rho_i, double rho_j, double P_i, double P_j,
                             const Eigen::Vector2d& gradW, double Visc)
{
    Eigen::Vector2d Acc;
    
    Acc.x() = -mass * ( P_i/(rho_i*rho_i) + P_j/(rho_j*rho_j) ) * gradW.x();
    Acc.y() = -mass * ( P_i/(rho_i*rho_i) + P_j/(rho_j*rho_j) ) * gradW.y();
    //std::cout<<"\n acc: "<< Acc.transpose() << std::endl;

    return  Acc;
}


double Energy(double mass, double rho_i, double rho_j, double P_i,double  P_j, const Eigen::Vector2d& v_i,
                        const Eigen::Vector2d& v_j, const Eigen::Vector2d& gradW, Eigen::Vector2d diss) 
{
        return -mass * ( P_i*v_j/(rho_i*rho_i) + P_j*v_i/(rho_j*rho_j) ).dot(gradW); //0.5??
}



std::vector<Particle> System(std::vector<Particle>& mesh) {

    //std::vector<Particle> mesh = mesh1;
    std::vector<double> q;
    std::vector<Eigen::Vector2d> dq;
    
    for (int ii = 0; ii < N; ++ii) {
        for (int jj = ii + 1; jj < N; ++jj) { 

            Eigen::Vector2d r_ij = mesh[ii].r - mesh[jj].r;
            double h_ij = h_len(mesh[ii].mass, mesh[ii].rho); 
            double k = M6Kernel(r_ij, h_ij); 

            // Density sum
            mesh[ii].rho += mesh[jj].mass * k;
            mesh[jj].rho += mesh[ii].mass * k;

            Eigen::Vector2d dk = GradM6Kernel(r_ij, h_ij); 
            dq.push_back(dk); 
        }
    }

    // Update Pressure
    for (Particle& p : mesh) {
        p.P = Pressure(p.rho, p.e);
        //std::cout << "Pressure " << p.P << "\n";
    }

    double kk = 0;
    for (int ii = 0; ii < N; ++ii) {
        for (int jj = ii + 1; jj < N; ++jj) { 

            double p_visc = momentum_viscosity(mesh[ii].r, mesh[jj].r, mesh[ii].v, mesh[jj].v,
                                    mesh[ii].rho, mesh[jj].rho, mesh[ii].e, mesh[jj].e, h_len(mesh[ii].mass, mesh[ii].rho), h_len(mesh[jj].mass, mesh[jj].rho));

            Eigen::Vector2d e_visc = energy_disipation(mesh[ii].r, mesh[jj].r, mesh[ii].v, mesh[jj].v,
                                    mesh[ii].rho, mesh[jj].rho, mesh[ii].e, mesh[jj].e, h_len(mesh[ii].mass, mesh[ii].rho), h_len(mesh[jj].mass, mesh[jj].rho));

            mesh[ii].d_v += Acceleration( mesh[jj].mass, mesh[ii].rho, mesh[jj].rho, mesh[ii].P, mesh[jj].P, dq[kk], p_visc);
            mesh[jj].d_v -= Acceleration( mesh[ii].mass, mesh[jj].rho, mesh[ii].rho, mesh[jj].P, mesh[ii].P, dq[kk], p_visc);

            mesh[ii].d_e += Energy(mesh[jj].mass, mesh[ii].rho, mesh[jj].rho, mesh[ii].P, mesh[jj].P, mesh[ii].v, mesh[jj].v, dq[kk], e_visc);
            mesh[jj].d_e -= Energy(mesh[ii].mass, mesh[jj].rho, mesh[ii].rho, mesh[jj].P, mesh[ii].P, mesh[jj].v, mesh[ii].v, dq[kk], e_visc);
            kk++;
        }
    }


    for (Particle& p : mesh) {
        p.d_r = p.v; 
    }

    return mesh;
}


std::vector<Particle> Integration(std::vector<Particle>& mesh, double x1, double x2, double y1, double y2, int n) {

    double tstep = 0.0005;
    const double tmax = tstep * n;
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
        UpdateGhostParticles(mesh);
        //Boundary_Reflective(mesh, x1, x2);
        //Boundary_Periodic(mesh, x2);

        std::cout <<" i: "<< ii <<std::endl; 
        
        }

    return mesh;

}


std::vector<Particle> EulerIntegration(std::vector<Particle>& mesh, double x1, double x2, double y1, double y2) {
    double tstep = 0.0005;
    int NSteps = 400;
    
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
        Boundary_Periodic(mesh, x1, x2, y1, y2);
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
    double x1 = 0.0; double x2 = 1.0;
    double y_1 = 0.0; double y2 = 0.1;
    int steps = 1;
    std::vector<Particle> mesh = Mesh(x1,x2,y_1,y2);
    std::cout<<" Mesh created" << std::endl;
    AddGhostParticles(mesh, x1, x2, y_1, y2);
    //InitializeDensity(mesh);
    std::vector<Particle> int_mesh = Integration(mesh, x1, x2, y_1, y2, steps);
    //std::vector<Particle> int_mesh = EulerIntegration(mesh, x1, x2, y_1, y2);



    csv(int_mesh);

    return 0;
}
