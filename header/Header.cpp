#include "Header.hpp"

std::string Header::when(){
	Time t;
	std::string s("Created on");
	s += " " + my::tostring(t.day());
	s += "." + my::tostring(t.month());
	s += "." + my::tostring(t.year());
	s += " at " + my::tostring(t.hour());
	s += ":" + my::tostring(t.min());
	s += ":" + my::tostring(t.sec());
	return s;
}

void Header::init(std::string const& s){
	title(s,'=');
	textit(when());
	nl();
}

void Header::add(std::string const& s){
	text(s);
}

void Header::add(std::string const& s, double const& d){
	def(s + " = " + my::tostring(d), "double");
}

void Header::add(std::string const& s, bool const& d){
	if(d){ def(s + " = true", "bool"); }
	else { def(s + " = false", "bool"); }
}

void Header::add(std::string const& s, std::string const& d){
	def(s + " = " + d, "string");
}

void Header::add(std::string const& s, unsigned int const& d){
	def(s + " = " + my::tostring(d), "unsigned int");
}

void Header::add(std::string const& s, int const& d){
	def(s + " = " + my::tostring(d), "int");
}

void Header::add(std::string const& s, std::complex<double> const& d){
	def(s + " = " + my::tostring(d), "complex");
}

std::ostream& operator<<(std::ostream& flux, Header const& h){
	flux<<h.get();
	return flux;
}
