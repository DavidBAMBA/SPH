#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>


constexpr double PI = 3.14159265358979323846;
const double par = 1.0;
const double k = 0.1;
const double nu = 1.0;
const double lmbda  = 2.01;
const int Nx = 3;
const int Ny = 3;
const int Nz = 3;
const int N = Nx * Ny * Nz;
const double kappa = 400000000.0;

struct Particle {

    // Propierties
    Eigen::Vector3d r;
    Eigen::Vector3d v;
    double rho, P;
    double mass = 2.0/N;

    // Derivatives
    Eigen::Vector3d d_r = Eigen::Vector3d::Zero();
    Eigen::Vector3d d_v = Eigen::Vector3d::Zero();
    double d_rho = 0;
    double d_P = 0;

};

std::vector<Particle> Mesh(double x1, double x2, double y1, double y2, double z1, double z2) {
    // Calculate steps
    double xdim = x2 - x1;
    double ydim = y2 - y1;
    double zdim = z2 - z1;
    // The xstep is modified to ensure that hexagonal close packing is achieved
    double xstep = xdim / Nx;
    double ystep = ydim / ((Ny - 1) * sqrt(3) / 2); 
    double zstep = zdim / (Nz - 1); 

    // Create a mesh of particles
    std::vector<Particle> mesh;
    mesh.reserve(N);

    // Offsets for hexagonal close packing
    double xOffset = xstep / 2;
    double yOffset = ystep * sqrt(3) / 2;
    
    for (int kk = 1; kk <= Nz; ++kk) {
        for (int jj = 1; jj <= Ny; ++jj) {
            for (int ii = 1; ii <= Nx; ++ii) {
                Particle p;
                // Calculate position with offsets for hexagonal packing
                double x = x1 + (ii - 1) * xstep + (jj % 2) * xOffset;
                double y = y1 + (jj - 1) * yOffset;
                double z = z1 + (kk - 1) * zstep + ((kk % 2) * yOffset - yOffset / 2) * (ii % 2);
                
                p.r = Eigen::Vector3d(x, y, z);
                p.v = Eigen::Vector3d(0.0, 0.0, 0.0);
                p.rho = 0.0;
                p.P = 0.0;

                mesh.push_back(p);
            }
        }
    }
    return mesh;
}

double h_len(double mass, double rho){
    
    return 0.1;
}


double Kernel(const Eigen::Vector3d& r_ij, double h) {

    double r = r_ij.norm();

    double w = std::pow((1.0 / (h * std::sqrt(M_PI))),3) * std::exp((-r * r) / (h * h));

    return w;
}

Eigen::Vector3d GradKernel3D(const Eigen::Vector3d& r_ij, double h) 
{
    double r  = r_ij.norm();
    Eigen::Vector3d gradW;
    double n = -2 * std::exp(- (r*r) / (h*h) ) / (std::pow(h,5)*std::pow(M_PI,1.5));
    
    gradW.x() = n * r_ij.x();
    gradW.y() = n * r_ij.y();
    gradW.z() = n * r_ij.z();
    
    return gradW;
}


void nearest_neight(const std::vector<Particle>& mesh,
                    std::vector<int>& pair_i, std::vector<int>& pair_j,
                    std::vector<double>& q, std::vector<Eigen::Vector3d>& dq) {

    double h_ij;
    
    for (int ii = 0; ii < N; ++ii) {
        int count = 0;
        for (int jj = 0; jj < N; ++jj){
            // relative distance
            Eigen::Vector3d r_ij = mesh[ii].r - mesh[jj].r;

            h_ij = 0.5 * ( h_len(mesh[ii].mass, mesh[ii].rho) + h_len(mesh[jj].mass, mesh[jj].rho) );

            // Condition of nearest
            if (r_ij.norm() <= kappa * h_ij){
                pair_i.push_back(ii);
                pair_j.push_back(jj);
                q.push_back( Kernel(r_ij, h_ij) );
                dq.push_back( GradKernel3D(r_ij, h_ij) );
                count += 1;
            }
        }
        //std::cout<<"vecinos : "<<count<<std::endl;       
    }


}

void InitializeDensity(std::vector<Particle>& mesh) {
    for (auto& pi : mesh) {
        pi.rho = 0.0; 
        for (const auto& pj : mesh) {
            Eigen::Vector3d r_ij = pi.r - pj.r; 
            double h = h_len(pj.mass, pj.rho); 
            pi.rho += pj.mass * Kernel(r_ij, h); 
        }
    }
}


double Pressure(double rho) {
    return k * std::pow(rho,(1.0 + 1.0 / par));
}

Eigen::Vector3d gravForce(const Eigen::Vector3d& r) {
    return -lmbda * r;
}

Eigen::Vector3d viscosForce(const Eigen::Vector3d& v) {

    return -nu * v;
}

Eigen::Vector3d Acceleration(double mass,double rho_i, double rho_j, double P_i, double P_j,
                             const Eigen::Vector3d& gradW)
{
    Eigen::Vector3d Acc;
    Acc.x() = -mass * ( P_i/(rho_i*rho_i) + P_j/(rho_j*rho_j)) * gradW.x();
    Acc.y() = -mass * ( P_i/(rho_i*rho_i) + P_j/(rho_j*rho_j)) * gradW.y();
    Acc.z() = -mass * ( P_i/(rho_i*rho_i) + P_j/(rho_j*rho_j)) * gradW.z();

    return  Acc;
}


std::vector<Particle> System(std::vector<Particle>& mesh) {

    std::vector<int> pair_i;
    std::vector<int> pair_j;
    std::vector<double> q;
    std::vector<Eigen::Vector3d> dq;

    // Nearest Neighbors
    nearest_neight(mesh, pair_i, pair_j, q, dq);
    //std::cout<<"aqui"<<std::endl;
    int NPairs = pair_i.size();
    //std::cout<<"Npairs"<<pair_i.size()<<std::endl;

   
    //std::cout<<"aqui"<<std::endl;

    // Update Density
    for (int kk = 0; kk < NPairs; ++kk) {

        int pi = pair_i[kk];
        int pj = pair_j[kk];
        
        mesh[pi].rho += mesh[pj].mass * q[kk];
        mesh[pj].rho += mesh[pi].mass * q[kk];
    }
    //std::cout<<"aqui"<<std::endl;
    // Update Pressure
    for (Particle& p : mesh) {
        p.P = Pressure(p.rho);
    }

    // Calculate the System Equations 
    for (int kk = 0; kk < NPairs; ++kk) {
        int pi = pair_i[kk];
        int pj = pair_j[kk];

        mesh[pi].d_v += Acceleration( mesh[pj].mass, mesh[pi].rho, mesh[pj].rho, mesh[pi].P, mesh[pj].P, dq[kk]) + viscosForce(mesh[pi].v) + gravForce(mesh[pi].r);
        mesh[pj].d_v -= Acceleration( mesh[pi].mass, mesh[pj].rho, mesh[pi].rho, mesh[pj].P, mesh[pi].P, dq[kk]) - viscosForce(mesh[pj].v) - gravForce(mesh[pj].r);
        
        }

    for (Particle& p : mesh) {
        p.d_r = p.v; 
    }

    return mesh;
}


std::vector<Particle> Integration(std::vector<Particle>& mesh, double x1, double x2, int n) {

    double tstep = 0.005;
    const double tmax = tstep * n;
    const int    NSteps = static_cast<int>((tmax - tstep) / tstep);
    double t = 0.005;

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
        }

        k2 = System(int_mesh);

        // k3
        for (size_t ii = 0; ii < int_mesh.size(); ++ii) {
            int_mesh[ii].r = k2[ii].r + 0.5 * tstep * k2[ii].d_r;
            int_mesh[ii].v = k2[ii].v + 0.5 * tstep * k2[ii].d_v;
        }

        k3 = System(int_mesh);

        // k4
        for (size_t ii = 0; ii < int_mesh.size(); ++ii) {
            int_mesh[ii].r = k3[ii].r + tstep * k3[ii].d_r;
            int_mesh[ii].v = k3[ii].v + tstep * k3[ii].d_v;
        }

        k4 = System(int_mesh);

        // Update Mesh
        for (size_t ii = 0; ii < mesh.size(); ++ii) {
            mesh[ii].r += (tstep / 6.0) * (k1[ii].d_r + 2 * k2[ii].d_r + 2 * k3[ii].d_r + k4[ii].d_r);
            mesh[ii].v += (tstep / 6.0) * (k1[ii].d_v + 2 * k2[ii].d_v + 2 * k3[ii].d_v + k4[ii].d_v);
            mesh[ii].P = k4[ii].P; 
            mesh[ii].rho = k4[ii].rho;
        }

        t += tstep;

        std::cout <<"  bi: "<< ii <<std::endl; 
        
        }

    return mesh;

}

void csv(const std::vector<Particle>& mesh){
    std::ofstream file("mesh-3d.csv");

    // Escribe los encabezados en el archivo CSV
    file << "x,y,z,vx,vy,vz,rho,P\n";
    
    for (const auto& p : mesh) {
        file << p.r(0) << "," << p.r(1) << "," << p.r(2) << "," 
        << p.v(0) << "," << p.v(1) << "," << p.v(2) << "," 
        << p.rho << "," << p.P << "\n";
    }

    file.close();
}

int main(int argc, char* argv[]) {

     // Define the dimentions of tube
    int n = std::atoi(argv[1]);
    double x1 = -3.0, x2 = 3.0;
    double y1 = -3.0, y2 = 3.0;
    double z1 = -3.0, z2 = 3.0;

    std::vector<Particle> mesh = Mesh(x1,x2,y1,y2,z1,z2);
    std::cout<<" Mesh created" << std::endl;
    InitializeDensity(mesh);
    std::vector<Particle> int_mesh = Integration(mesh, x1, x2,n);

    csv(int_mesh);
    

    return 0;
}
