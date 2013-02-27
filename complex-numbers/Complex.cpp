#include"Complex.hpp"

// définitions des constructeurs
Complex::Complex()
	:a(0.0),b(0.0),disgard_imaginary_part(1e-10),disgard_real_part(1e-10)
{ }

Complex::Complex(Complex const& c)
	:a(c.Re()),b(c.Im()),disgard_imaginary_part(1e-10),disgard_real_part(1e-10)
{ }

Complex::Complex(double const& a, double const& b)
	:a(a),b(b),disgard_imaginary_part(1e-10),disgard_real_part(1e-10)
{ }

// accesseurs à la partie réelle et imaginaire
double const& Complex::Re() const {return a;}

double const& Complex::Im() const {return b;}

// méthodes diverses
Complex Complex::conj() const{
	return Complex(a,-b);
}
double Complex::module() const{
	return sqrt(a*a+b*b);
}
double const& Complex::ProjRe() const {
	if(b > disgard_imaginary_part){
		std::cerr<<"imaginary part too large to be disgarded safly"<<std::endl;
	}
	return a;	
}
double const& Complex::ProjIm() const {
	if(a > disgard_real_part){
		std::cerr<<"real part too large to be disgarded safly"<<std::endl;
	}
	return b;	
}

// surcharge d'opérateurs internes
Complex& Complex::operator+=(Complex const& c){
	a += c.a;
	b += c.b;
	return *this;
}
Complex& Complex::operator-=(Complex const& c){
	a -= c.a;
	b -= c.b;
	return *this;
}
Complex& Complex::operator*=(Complex const& c){
	double a_old(a);
	a = a*c.a - b*c.b;
	b = b*c.a + a_old*c.b;
	return *this;
}
Complex& Complex::operator/=(Complex const& c){
	double tmp(c.module());
	double a_old(a);
	tmp *=tmp;
	a = (a*c.a + b*c.b)/tmp;
	b = (b*c.a - a_old*c.b)/tmp;
	return *this;
}
Complex& Complex::operator*=(double const& r){
	a *= r;
	b *= r;
	return *this;
}
Complex& Complex::operator/=(double const& r){
	a /= r;
	b /= r;
	return *this;
}

// surcharge d'opérateurs externes
std::ostream& operator<<(std::ostream& flux, Complex const& c){
	if (c.Im()<0.0){
		flux << c.Re()<<" - "<<std::abs(c.Im())<<"i";
	} else{
		flux << c.Re()<<" + "<<c.Im()<<"i";
	}
	return flux;
}
Complex operator+(Complex const& c1, Complex const& c2){
	Complex c_out(c1);
	c_out += c2;
	return c_out;
}
Complex operator-(Complex const& c1, Complex const& c2){
	Complex c_out(c2); // attention ici on veut que c = c1-c2...
	c_out -= c1;
	return c_out;
}
Complex operator*(Complex const& c1, Complex const& c2){
	Complex c_out(c1);
	c_out *= c2;
	return c_out;
}
Complex operator/(Complex const& c1,Complex const& c2){
	Complex c_out(c1);
	c_out /= c2;
	return c_out;
}
Complex operator*(Complex const& c, double const& r){
	return r*c;
}
Complex operator/(Complex const& c, double const& r){ 
	Complex c_out(c);
	c_out /= r;
	return c_out;
}
Complex operator*(double const& r, Complex const& c){
	Complex c_out(c);
	c_out *= r;
	return c_out;
}
Complex operator/(double const& r, Complex const& c){
	double tmp(c.module());
	tmp *= tmp;
	return r * c.conj() / tmp;
}

// fonctions liées à la classe
Complex exp(Complex const& c){
	Complex c_cout(cos(c.Im()),sin(c.Im()));
	return exp(c.Re())*c_cout ;
}
Complex ii(double const& r){
	return Complex(0,r);
}
Complex ii(Complex const& c){
	return Complex(-c.Im(),c.Re());
}


