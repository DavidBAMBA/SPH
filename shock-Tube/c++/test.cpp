#include <iostream>
#include <cmath>

// Function to apply periodic boundary conditions in 1D
double applyPeriodicBoundary1D(double position, double box_length) {
    // The fmod function is used to apply the wrap-around effect
    position = std::fmod(position, box_length);
    std::cout<<position<<std::endl;
    // If position is negative, wrap it to the other end of the box
    if (position < 0) {
        position += box_length;
    }
    return position;
}

int main() {
    // Define the size of the 1D box
    double box_length = 1.2;
    
    // Example position of a particle in 1D
    double position = 1.5;
    
    // Apply periodic boundary conditions
    double wrapped_position = applyPeriodicBoundary1D(position, box_length);
    
    // Output the result
    std::cout << "The wrapped position is: " << wrapped_position << std::endl;
    
    return 0;
}
