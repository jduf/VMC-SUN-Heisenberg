#ifndef DEF_MATRICE
#define DEF_MATRICE


#include <iostream>
#include <cmath>

#include "Vecteur.hpp"

class Matrice{
	public:
		Matrice();
		Matrice(unsigned int N);
		Matrice(unsigned int N, double val);
		Matrice(Matrice const& mat);
		~Matrice();

		void Test() const;
		inline unsigned int size() const{ return N; };
		void Chop();
		Matrice Transpose() const;

		Matrice& operator=(Matrice const& mat);
		Matrice& operator*=(Matrice const& mat);
		double const& operator()(unsigned int i, unsigned int j) const;
		double& operator()(unsigned int i, unsigned int j);

		double* ptr(){ return m; };

		void eigenvalue(Vecteur& EVal);
		void print() const;

	private:
		double *m;
		unsigned int N;
		
		void fill_matrice(double val);
};

std::ostream& operator<<(std::ostream& flux, Matrice const& mat);
Matrice operator*(Matrice const& mat1, Matrice const& mat2);

extern "C" void dsyev_(char *jobz,
		char *uplo,
		unsigned int *n,
		double *a,
		unsigned int *lda, 
		double *w, 
		double *work, 
		unsigned int *lwork, 
		int *info);
#endif
