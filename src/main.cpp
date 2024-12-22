#include <iostream>
#include "Vector.hpp"
#include "Matrix.hpp"


int main(void){

    Vector<float, 3> v1 = {1, 2, 3};
    Vector<float, 3> v2 = {4, 5, 6};

    std::cout << "v1: " << v1 << std::endl;
    std::cout << "v2: " << v2 << std::endl;

    return 0;
}