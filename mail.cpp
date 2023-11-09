#include <iostream>
#include "sim202.hpp"
#include "sim202.cpp"


int main()
{
    
    std::cout << "DEMARAGE DU TEST\n";
    
    std::cout << "Test de géométrie.h\n";
    
    std::cout << "Test de la classe vecteur\n";
    vecteur A(1,2);
    vecteur B(A);
    
    std::cout << vecteur() << vecteur(1) << vecteur(1,2) << vecteur(1,2,3) << B;
    
    std::cout << rotation_d2(vecteur(1 ,0), PI/2);
    
    std::cout << "\n" << segment(sommet(0,0),sommet(1,0));
    
    std::cout << "Test de scene\n";
    
    polygone poly(4,vecteur(2,0),2);
    for (int i =0; i<4; i++) {
        std::cout << poly.segments[i];
    };
    
    std::cout << EPSILON;
    
    // Test de Translation et dilatation
//    
//    sommet S(1,1);
//    vecteur v(0,1);
//    
//    std::cout << translation(v, S);
    
    poly.transformation_poly(vecteur(-2,0), 2, PI/4);
    for (int i =0; i<4; i++) {
        std::cout << poly.segments[i];
    };
    
    std::cout << poly;
    
    vector<sommet> vectsom;
    vectsom.resize(4);
    vectsom[0]=sommet(-5,0);
    vectsom[1]=sommet(-6,-1);
    vectsom[2]=sommet(-4,0);
    vectsom[3]=sommet(-6,1);
    
    std::cout << polygone(vectsom);
    return 0;
}
