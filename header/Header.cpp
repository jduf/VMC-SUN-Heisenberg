#include "Header.hpp"

Header::Header(){ }

Header::~Header(){}

std::string Header::when(){
	Time t;
	std::string s("Created on");
	s += " " + tostring(t.day());
	s += "." + tostring(t.month());
	s += "." + tostring(t.year());
	s += " at " + tostring(t.hour());
	s += ":" + tostring(t.min());
	s += ":" + tostring(t.sec());
	return s;
}

void Header::init(std::string t){
	title(t,"=");
	textit(when());
	nl();
}

void Header::add(std::string const& s){
	text(s);
}

void Header::add(std::string const& s, double const& d){
	def(s + " = " + tostring(d), "double");
}

void Header::add(std::string const& s, bool const& d){
	if(d){def(s + " = true", "bool");}
	else{def(s + " = false", "bool");}
}

void Header::add(std::string const& s, std::string const& d){
	def(s + " = " + d, "string");
}

void Header::add(std::string const& s, unsigned int const& d){
	def(s + " = " + tostring(d), "unsigned int");
}

void Header::add(std::string const& s, int const& d){
	def(s + " = " + tostring(d), "int");
}

void Header::add(std::string const& s, std::complex<double> const& d){
	def(s + " = " + tostring(d), "complex");
}

//void Header::add(std::string const& s, Vector<unsigned int> const& vec){
	//def(s + "(" + tostring(vec.size()) + ")","Vector<unsigned int>");
//}
//
//void Header::add(std::string const& s, Vector<double> const& vec){
	//def(s + "(" + tostring(vec.size()) + ")","Vector<double>");
//}
//
//void Header::add(std::string const& s, Matrix<int> const& mat){
	//def(s + "(" + tostring(mat.row()) + "x" + tostring(mat.col())+ ")","Matrix<int>");
//}
//
//void Header::add(std::string const& s, Matrix<unsigned int> const& mat){
	//def(s + "(" + tostring(mat.row()) + "x" + tostring(mat.col())+ ")","Matrix<unsigned int>");
//}
//
//void Header::add(std::string const& s, Matrix<double> const& mat){
	//def(s + "(" + tostring(mat.row()) + "x" + tostring(mat.col())+ ")","Matrix<double>");
//}
//
//void Header::add(std::string const& s, Matrix<std::complex<double> > const& mat){
	//def(s + "(" + tostring(mat.row()) + "x" + tostring(mat.col())+ ")","Matrix<complex>");
//}
//

std::ostream& operator<<(std::ostream& flux, Header const& h){
	flux<<h.get();
	return flux;
}
