#ifndef DEF_VECTEUR
#define DEF_VECTEUR


#include <iostream>
#include <cmath>

class Vecteur{
	public:
/*Constructors and destructor*/
		Vecteur();
		Vecteur(unsigned int N);
		Vecteur(unsigned int N, double val);
		Vecteur(Vecteur const& mat);
		~Vecteur();

/*operators*/
		Vecteur& operator=(Vecteur const& mat);
		inline double const& operator()(unsigned int const& i) const { return v[i]; };
		inline double& operator()(unsigned int const& i) { return v[i]; };
		//Vecteur& operator*=(double const& d);

/*methods that modify the class*/
		void chop();

/*other methods*/
		inline double* ptr() const { return v; };
		inline unsigned int size() const { return N; };

		void test() const;
		void print() const;

	private:
		double *v;
		unsigned int N;
		
/*methods that modify the class*/
		void fill_vecteur(double val);
};

//double operator*(Vecteur const& vec1, Vecteur const& vec2);
#endif
