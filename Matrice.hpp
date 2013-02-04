#ifndef DEF_MATRICE
#define DEF_MATRICE

#include <iostream>
#include <iomanip>
#include <cmath>

#include"Vecteur.hpp"


class Matrice{
	public:
/*Constructors and destructor*/
		Matrice();
		Matrice(unsigned int N);
		Matrice(unsigned int N, double val);
		Matrice(Matrice const& mat);
		~Matrice();

/*operators*/
		Matrice& operator=(Matrice const& mat);
		double const& operator()(unsigned int const& i, unsigned int const& j) const;
		double& operator()(unsigned int const& i, unsigned int const& j);
		Matrice& operator*=(Matrice const& mat); // m1 *= m2 : m1 = m1*m2
		Matrice& operator+=(Matrice const& mat); 
		Matrice& operator-=(Matrice const& mat); 
		
/*methods that modify the class*/
		void chop();

/*methods that return something related to the class*/
		Matrice transpose() const;
		
/*other methods*/
		inline double* ptr() const { return m; };
		inline unsigned int size() const { return N; };

		void print() const;

	private:
		double *m; // m = [[ligne0],[ligne1],...]
		unsigned int N;
		
/*methods that modify the class*/
		void fill_matrice(double val);

};

std::ostream& operator<<(std::ostream& flux, Matrice const& mat);
Matrice operator*(Matrice const& mat1, Matrice const& mat2);
Matrice operator+(Matrice const& mat1, Matrice const& mat2);
Matrice operator-(Matrice const& mat1, Matrice const& mat2);
Matrice operator^(Vecteur const& vec1, Vecteur const& vec2);
#endif
