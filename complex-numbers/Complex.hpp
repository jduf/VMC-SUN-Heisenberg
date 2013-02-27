#ifndef DEF_COMPLEX
#define DEF_COMPLEX

#include<iostream>
#include<cmath>
class Complex{
	private:
		double a,b;
		double const disgard_imaginary_part, disgard_real_part;

	public:
// définitions des constructeurs
		Complex();
		Complex(Complex const& c);
		Complex(double const& a, double const& b);

// accesseurs à la partie réelle et imaginaire
		double const& Re() const;
		double const& Im() const;

// méthodes diverses
		Complex conj() const;
		double module() const;
		double const& ProjRe() const;
		double const& ProjIm() const;

// surcharge d'opérateurs internes
		Complex& operator+=(Complex const& c);
		Complex& operator-=(Complex const& c);
		Complex& operator*=(Complex const& c);
		Complex& operator/=(Complex const& c);
		Complex& operator*=(double const& r);
		Complex& operator/=(double const& r);
};

// surcharge d'opérateurs externes
std::ostream& operator<<(std::ostream &flux, Complex const& c);
Complex operator+(Complex const& c1, Complex const& c2);
Complex operator-(Complex const& c1, Complex const& c2);
Complex operator*(Complex const& c1, Complex const& c2);
Complex operator/(Complex const& c1, Complex const& c2);
Complex operator*(Complex const& c, double const& r);
Complex operator/(Complex const& c, double const& r);
Complex operator*(double const& r, Complex const& c);
Complex operator/(double const& r, Complex const& c);

// fonctions liées à la classe
Complex exp(Complex const& c);
Complex ii(double const& r);
Complex ii(Complex const& r);
#endif
