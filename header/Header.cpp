#include "Header.hpp"

Header::Header(){ }

Header::~Header(){}

void Header::add(std::string const& s){
	rst.text(s);
}

void Header::add(std::string const& s, double const& d){
	rst.def(s + " = " + tostring(d), "double");
}

void Header::add(std::string const& s, bool const& d){
	if(d){rst.def(s + " = true", "bool");}
	else{rst.def(s + " = false", "bool");}
}

void Header::add(std::string const& s, std::string const& d){
	rst.def(s + " = " + d, "string");
}

void Header::add(std::string const& s, unsigned int const& d){
	rst.def(s + " = " + tostring(d), "unsigned int");
}

void Header::add(std::string const& s, int const& d){
	rst.def(s + " = " + tostring(d), "int");
}

void Header::add(std::string const& s, std::complex<double> const& d){
	rst.def(s + " = " + tostring(d), "complex");
}

void Header::add(std::string const& s, Vector<unsigned int> const& vec){
	rst.def(s + "(" + tostring(vec.size()) + ")","Vector<unsigned int>");
}

void Header::add(std::string const& s, Vector<double> const& vec){
	rst.def(s + "(" + tostring(vec.size()) + ")","Vector<double>");
}

void Header::add(std::string const& s, Matrix<int> const& mat){
	rst.def(s + "(" + tostring(mat.row()) + "x" + tostring(mat.col())+ ")","Matrix<int>");
}

void Header::add(std::string const& s, Matrix<unsigned int> const& mat){
	rst.def(s + "(" + tostring(mat.row()) + "x" + tostring(mat.col())+ ")","Matrix<unsigned int>");
}

void Header::add(std::string const& s, Matrix<double> const& mat){
	rst.def(s + "(" + tostring(mat.row()) + "x" + tostring(mat.col())+ ")","Matrix<double>");
}

void Header::add(std::string const& s, Matrix<std::complex<double> > const& mat){
	rst.def(s + "(" + tostring(mat.row()) + "x" + tostring(mat.col())+ ")","Matrix<complex>");
}

std::string Header::when(){
	Time t;
	std::string h;
	h = "Created on";
	h += " " + tostring(t.day());
	h += "." + tostring(t.month());
	h += "." + tostring(t.year());
	h += " at " + tostring(t.hour());
	h += ":" + tostring(t.min());
	h += ":" + tostring(t.sec());
	return h;
}

void Header::init(std::string title){
	rst.title(title,"=");
	rst.textit(when());
	rst.np();
}

std::ostream& operator<<(std::ostream& flux, Header const& h){
	flux<<h.get();
	return flux;
}
