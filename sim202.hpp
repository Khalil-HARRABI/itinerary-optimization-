#ifndef Planification_de_trajectoire
#define Planification_de_trajectoire

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>


using namespace std;


const double PI= 3.141592653589793238462643383279;
const double EPSILON = 1.e-12;

//=================================================================================
//                        class sommet
//=================================================================================

class sommet
{ public:
    double x,y;  //2D sommet coordinates
    sommet(double xi, double yi)  : x(xi),y(yi){}
	sommet& operator+=(const sommet& p);
    sommet& operator-=(const sommet& p);
    sommet& operator*=(double a);
};
sommet operator+(const sommet& p1, const sommet& p2);
sommet operator-(const sommet& p1, const sommet& p2);
sommet operator*(double a, const sommet& p);
double operator*(const sommet& p1, const sommet& p2); // produit vectoriel 2D
double operator|(const sommet& p1, const sommet& p2); // produit scalaire  2D
double norme(const sommet& p);                       // norme euclidienne
ostream& operator<<(ostream& out, const sommet& p);

//=================================================================================
//                        class segment
//=================================================================================
class segment;
typedef segment sommet;

class segment {
    public :
    double* sommet;
// CONSTRUCTEURS
    segment(); // Point de dimension 2 par défaut nul
    segment(double x,double y);
    segment(const segment& s);
    segment(const sommet& v1,const sommet& v2);
    segment& operator= (const segment& s);
// AFFICHAGE
    friend std::ostream& operator <<(std::ostream &,const segment &);
    void print();
// ECRITURE
    void print_fichier(std::ostream&) const;
// COMBINAISONS LINEAIRES
    segment & operator +=(const segment & v);
    segment & operator +=(const double);
    segment & operator -=(const segment & v);
    segment & operator -=(const double);
    segment & operator *=(const double);
    segment & operator /=(const double);
// ACCES
    double& operator ()(int i) {return sommet[i-1];}
    double& operator [](int i) {return sommet[i];}
    double operator ()(int i) const{return sommet[i-1];}
    double operator [](int i) const{return sommet[i];}
};


// FONCTIONS DE COMBINAISONS DE segmentS
bool operator==(const segment & a,const segment & b);
bool operator!=(const segment & a,const segment & b);
segment operator +(const segment &,const segment &);
segment operator +(const segment &,const double & lambda);
segment operator +(const double & lambda,const segment &);
segment operator -(const segment &,const segment &);
segment operator -(const segment &,const double & lambda);
segment operator -(const double & lambda,const segment &);
segment operator *(const segment &,const double & lambda);
segment operator *(const double & lambda,const segment &);
segment operator /(const segment &,const double & lambda);

// FONCTIONS EUCLIDIENNES
double ps ( const segment & ,const segment &);
double norm (const segment & v);
double det_d2(const segment & v1, const segment & v2);
double maxim(double,double);
segment rotation_d2 (const segment& v, double theta);
sommet translation (const segment& v, const sommet& s);
sommet dillatation (const sommet& centre, const sommet& s,double dill);

//=================================================================================
//                        class obstacle
//=================================================================================
class obstacle {
    // On tourne dans le sens trigo
public:
    int nb_sommet;
    vector<sommet> sommets;
    vector<segment> segments;

    // CONSTRUCTEUR
// obstacle régulier, centre ve, n >=3cotés,
    obstacle(int n,const segment& v ,double d);
    // A partir d'une liste de points
    obstacle(const vector<sommet>&);
    // Translation, dilatation, rotation par rapport au centre de gravité du obstacle
    void transformation_poly(const segment&,double dill ,double theta);
    void remplissage_segm();

    // AFFICHAGE
    friend std::ostream& operator <<(std::ostream &,const obstacle &);
    // ECRITURE
    void print_fichier(std::ostream&) const;
};

//=================================================================================
//                        class scene
//=================================================================================
class scene {
public:
    int nb_obstacle;
    sommet depart;
    sommet objectif;
    vector<polygone> obstacles;

    // ECRITURE DANS UN FICHIER
    void exporte(string titre);

    // LECTURE D'UN FICHIER
    void importe(string titre);


};

    

//=================================================================================
//                        class graphe
//=================================================================================
class graphe {
public:
    // La matrice de graphe est ordonnée comme suit :
    int dim;
    double** dist;
    graphe(const graphe&);
    graphe(const scene& );
    double& distance(int i,int j){return dist[i][j];}
    friend std::ostream & operator <<(std::ostream &,const graphe &);
    
};

/* FONCTIONS INTERMEDIAIRES pour construire un graphe à partir d'une scène */

// Savoir si un segment est sur le chemin d'une trajectoire (point plus vecteur)
bool intersection (const sommet&,const sommet&,const segment&);

// Savoir si il y a un obstacle entre deux points d'une scène.
bool intersection_totale(const scene&, const sommet& source,int,int, const sommet& arrivee,int,int);
// Calcule la dimention de la matrice pour une scene
int calcule_dimension (const scene&);

// Initialise la matrice dist a la dimention dim
void initialise(double** dist, int dim);

// Vérifie si le sommet i est accessible par le sommet s d'un même polygone
bool accessible_sur_soimeme(const scene& scn, const polygone& p,int numpol,int i,int s);

// Vérifie si le sommet i du polygone p est accessible depuis le sommet s (du polygone [int], numéro [int])
bool accessible_sur_autre(const scene& scn,const polygone& p,int numpol,int i,const sommet s,int s_p,int s_n);

/* CALCUL DU CHEMIN */

// Calcule le chemin : Donne le vecteur des predecesseur
void calcule_chemin (const graphe&, int*);

// minimum partiel
int minimum (const double*,const bool*,int);

// Retraite une liste des prédécesseurs en un tableau de sommets consécutifs, le premier étant le nombre de sommets de la chaine, le deuxième 0, et le dernier la sortie
int* liste_sommet (int* brut,int n);

void cherche_coord(int nb_som, int dim, int& poly,int& som,int num,const scene& scn);

// Ajoute a un fichier Texte la liste des sommets
void ajoute_au_fichier(int* sol,int dim,const scene& scn, string titre);
#endif /* defined(__Projet_Xcode___Planification_de_trajectoire__graph__) */

/* FONCTION FINALE, PREND UNE SCENE, ET LA TRANSFORME EN .txt */

void calcule_le_plus_court_chemin (scene& scn, string titre);

void calcule_le_plus_court_chemin_padding_cercle ( scene scn, string titre,double rayon);