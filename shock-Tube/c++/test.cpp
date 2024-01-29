#include <iostream>
#include <cmath>
#include <Eigen/Dense>


Eigen::Vector2d GradM6QuinticKernel2D(const Eigen::Vector2d& r_ij, double h) {
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



int main() {
    
    Eigen::Vector2d r  = Eigen::Vector2d(0.25,0.0);
    double h = 0.1;
    Eigen::Vector2d q = GradM6QuinticKernel2D(r, h);

    std::cout<< "value: " << q.transpose() << std::endl;


    
    return 0;
}
