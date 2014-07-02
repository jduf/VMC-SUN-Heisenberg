#include "Header.hpp"

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

void Header::init(std::string const& s){
	title(s,"=");
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

std::ostream& operator<<(std::ostream& flux, Header const& h){
	flux<<h.get();
	return flux;
}
