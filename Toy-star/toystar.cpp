#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

constexpr double EPSILON = 1e-10; // Small value to prevent division by zero
constexpr double PI = 3.14159265358979323846;

struct Particle {
    double mass;
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;

    Particle(double m, const Eigen::Vector3d& pos, const Eigen::Vector3d& vel)
        : mass(m), position(pos), velocity(vel) {}
};

Eigen::Vector3d W(const Eigen::Vector3d& r, double h) {
    double r_norm = r.norm();
    double factor = 1.0 / (h * sqrt(PI)); // Use the constant PI value
    double exp_term = exp(-(r_norm * r_norm) / (h * h));
    return (factor * factor * factor) * exp_term * r;
}

Eigen::Vector3d gradW(const Eigen::Vector3d& r, double h) {
    double r_norm = r.norm();
    double dw = -2.0 * exp(-r_norm * r_norm / (h * h)) / (pow(h, 5) * pow(PI, 1.5)); // Use the constant PI value
    return dw * r / r_norm;
}

Eigen::VectorXd Density(const std::vector<Particle>& particles, double h) {
    int N = particles.size();
    Eigen::VectorXd rho(N);
    rho.setZero();

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            Eigen::Vector3d dx = particles[i].position - particles[j].position;
            rho(i) += particles[j].mass * W(dx, h).norm();
        }
    }
    return rho;
}

double individualPressure(double rho_value, double k, double n) {
    return k * pow(rho_value + EPSILON, (1.0 + 1.0 / n));
}

Eigen::Vector3d gravForce(double lmbda, const Eigen::Vector3d& r) {
    return -lmbda * r;
}

Eigen::Vector3d viscosForce(double nu, const Eigen::Vector3d& v) {
    return -nu * v;
}

Eigen::VectorXd Acceleration(const std::vector<Particle>& particles, double h, double k, double n, double lmbda, double nu) {
    int N = particles.size();
    Eigen::VectorXd acc(3 * N);
    acc.setZero();

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            Eigen::Vector3d dx = particles[i].position - particles[j].position;
            Eigen::Vector3d dw = gradW(dx, h);
            
            double P_i = individualPressure(rho(i), k, n);
            double P_j = individualPressure(rho(j), k, n);

            Eigen::Vector3d acc_pressure = -particles[j].mass *
                ((P_j / pow(rho(j) + EPSILON, 2)) + (P_i / pow(rho(i) + EPSILON, 2))) * dw;

            acc.segment(3 * i, 3) += acc_pressure;
        }

        acc.segment(3 * i, 3) += gravForce(lmbda, particles[i].position) + viscosForce(nu, particles[i].velocity);
    }
    return acc;
}

int main() {
    int N = 1000;
    double tEnd = 12.0;
    double dt = 0.04;
    double M = 2.0;
    double R = 0.75;
    double h = 0.1;
    double k = 0.1;
    double n = 1.0;
    double nu = 2.0;
    double lmbda = 2.01;

    std::vector<Particle> particles;

    // Initialize particles with random positions within the star's radius and zero velocities.
    for (int i = 0; i < N; ++i) {
        Eigen::Vector3d rand_pos;
        do {
            rand_pos = R * Eigen::Vector3d::Random();
        } while (rand_pos.norm() > R);

        particles.push_back(Particle(M / N, rand_pos, Eigen::Vector3d::Zero()));
    }



    // Simulation loop with leapfrog integration
    int numTimesteps = static_cast<int>(tEnd / dt);
    std::vector<Eigen::VectorXd> densityData(numTimesteps, Eigen::VectorXd(N));

    for (int step = 0; step < numTimesteps; ++step) {
        Eigen::VectorXd acc = Acceleration(particles, h, k, n, lmbda, nu);

        for (int i = 0; i < N; ++i) {
            particles[i].velocity += acc.segment(3 * i, 3) * (dt / 2);
            particles[i].position += particles[i].velocity * dt;
        }

        acc = Acceleration(particles, h, k, n, lmbda, nu);

        for (int i = 0; i < N; ++i) {
            particles[i].velocity += acc.segment(3 * i, 3) * (dt / 2);
        }



        densityData[step] = Density(particles, h);

        std::cout << "\rStep # " << step << std::flush;
    }

    return 0;
}
