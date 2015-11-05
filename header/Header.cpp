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

void Header::add(std::string const& name, double const& val){
	def(name + " = " + my::tostring(val), "double");
}

void Header::add(std::string const& name, bool const& val){
	if(val){ def(name + " = true", "bool"); }
	else { def(name + " = false", "bool"); }
}

void Header::add(std::string const& name, std::string const& val){
	def(name + " = " + val, "string");
}

void Header::add(std::string const& name, unsigned int const& val){
	def(name + " = " + my::tostring(val), "unsigned int");
}

void Header::add(std::string const& name, int const& val){
	def(name + " = " + my::tostring(val), "int");
}

void Header::add(std::string const& name, std::complex<double> const& val){
	def(name + " = " + my::tostring(val), "complex");
}

std::ostream& operator<<(std::ostream& flux, Header const& h){
	flux<<h.get();
	return flux;
}
