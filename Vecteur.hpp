#ifndef DEF_VECTEUR
#define DEF_VECTEUR


#include <iostream>
#include <cmath>

class Vecteur{
	public:
		Vecteur();
		Vecteur(unsigned int N);
		Vecteur(unsigned int N, double val);
		Vecteur(Vecteur const& mat);
		~Vecteur();

		void Test() const;
		inline unsigned int size() const{ return N; };
		void Chop();

		Vecteur& operator=(Vecteur const& mat);
		//Vecteur& operator*=(double const& d);
		double const& operator()(unsigned int i) const;
		double& operator()(unsigned int i);

		double* ptr(){ return v; };
		void print() const;

	private:
		double *v;
		unsigned int N;
		
		void fill_vecteur(double val);
};

//double operator*(Vecteur const& vec1, Vecteur const& vec2);
#endif
