#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <Eigen/Dense>
#include <Eigen/Core>

const double pi = 3.14159265359;

struct Particle {
    double mass;
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;

    // Constructor
    Particle(double m, const Eigen::Vector3d& pos, const Eigen::Vector3d& vel)
        : mass(m), position(pos), velocity(vel) {}
};

// 3D Gaussian Smoothing kernel
double W(double r, double h) {
    return (1.0 / (h * sqrt(pi))) * (1.0 / (h * sqrt(pi))) * (1.0 / (h * sqrt(pi))) *
           exp(-r * r / (h * h));
}

// Gradient of the 3D Gaussian Smoothing kernel
Eigen::Vector3d gradW(const Eigen::Vector3d& r, double h) {
    double norm_r = r.norm();
    double dw = -2.0 * exp(-norm_r * norm_r / (h * h)) / (pow(h, 5) * pow(pi, 1.5));
    return dw * r / norm_r;
}

// Pairwise separations
Eigen::Vector3d pairwiseSeparations(const Particle& p1, const Particle& p2) {
    return p1.position - p2.position;
}

// Density
double density(const Particle& particle, const std::vector<Particle>& particles, double h) {
    double sum_w = 0.0;
    for (const auto& neighbor : particles) {
        Eigen::Vector3d dx = pairwiseSeparations(particle, neighbor);
        double r = dx.norm();
        double w = W(r, h);
        sum_w += neighbor.mass * w;
    }
    return sum_w;
}

// Equation of State
double pressure(double rho, double k, double n) {
    return k * pow(rho, (1.0 + 1.0/n));
}

// Acceleration
Eigen::Vector3d acceleration(const Particle& particle, const std::vector<Particle>& particles,
                             double m, double h, double k, double n, double lmbda, double nu) {
    Eigen::Vector3d acc_pressure(0.0, 0.0, 0.0);
    Eigen::Vector3d acc_gravity = -lmbda * particle.position;
    Eigen::Vector3d acc_viscosity = -nu * particle.velocity;

    double rho = density(particle, particles, h);
    double P = pressure(rho, k, n);

    for (const auto& neighbor : particles) {
        Eigen::Vector3d dx = pairwiseSeparations(particle, neighbor);
        double r = dx.norm();
        Eigen::Vector3d dw = gradW(dx, h);

        double neighbor_rho = density(neighbor, particles, h);
        double neighbor_P = pressure(neighbor_rho, k, n);

        acc_pressure -= neighbor.mass * (P/rho/rho + neighbor_P/neighbor_rho/neighbor_rho) * dw;
    }

    return acc_pressure + acc_gravity + acc_viscosity;
}

int main() {
    // Simulation parameters
    int    N    = 1000000;       // Number of p
    double t    = 0.0;     // current time of the simulation
    double tEnd = 12.0; // time at which simulation ends
    double dt   = 0.04;   // timestep
    double M    = 2.0;     // star mass
    double R    = 0.75;    // star radius
    double h    = 0.1;     // smoothing length
    double k    = 0.1;     // equation of state constant
    double n    = 1.0;     // polytropic index
    double nu   = 2.0;    // damping


    std::vector<Particle> particles;
    srand(static_cast<unsigned>(time(nullptr)));
    
    // Create some particles based on your parameters
    for (int i = 0; i < N; ++i) {
        Eigen::Vector3d pos(rand() % static_cast<int>(R * 1000) / 1000.0, 
                            rand() % static_cast<int>(R * 1000) / 1000.0, 
                            rand() % static_cast<int>(R * 1000) / 1000.0);
        Eigen::Vector3d vel(0, 0, 0);  // Assuming particles are stationary initially
        particles.push_back(Particle(1.0, pos, vel));
    }

    // Simulate the system over some timesteps
    int numTimesteps = 100;
    //double dt = 0.01;

    for (int step = 0; step < numTimesteps; ++step) {
        // Here you would call a function like timeIntegration() or any other to evolve your system
        // timeIntegration(particles, dt, ...); 

        // Output for visualization or analysis
        if (step % 10 == 0) {  // e.g., output every 10 steps
            std::cout << "Timestep: " << step << std::endl;
            for (const auto& p : particles) {
                std::cout << "Position: " << p.position.transpose() 
                          << " Velocity: " << p.velocity.transpose() << std::endl;
            }
        }
    }

    return 0;
}
