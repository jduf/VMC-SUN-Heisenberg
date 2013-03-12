#include "Header.hpp"

Header::Header():
	rst(new RST)
{ }

Header::~Header(){delete rst; }

void Header::add(std::string const& s, double const& d){
	rst->def(s + "=" + tostring(d), "double");
}

void Header::add(std::string const& s, bool const& d){
	if(d){rst->def(s + "=true", "bool");}
	else{rst->def(s + "=false", "bool");}
}

void Header::add(std::string const& s, std::string const& d){
	rst->def(s + "=" + d, "string");
}

void Header::add(std::string const& s, unsigned int const& d){
	rst->def(s + "=" + tostring(d), "unsigned int");
}

void Header::add(std::string const& s, std::complex<double> const& d){
	rst->def(s + "=" + tostring(d), "complex");
}

void Header::add(std::string const& s, Matrice<double> const& mat){
	rst->def(s + "(" + tostring(mat.size()) + "x" + tostring(mat.size())+ ")","Matrice<double>");
}

void Header::add(std::string const& s, Matrice<std::complex<double> > const& mat){
	rst->def(s + "(" + tostring(mat.size()) + "x" + tostring(mat.size())+ ")","Matrice<complex>");
}

void Header::add(std::string const& s, Array2D<unsigned int> const& arr){
	rst->def(s + "(" + tostring(arr.row()) + "x" + tostring(arr.col())+ ")","Array2D<unsigned int>");
}

std::string Header::when(){
	time_t t(time(0));
	struct tm* now(localtime(&t));
	std::string h;
	h = "Created on";
	h += " " + tostring(now->tm_mday);
	h += "." + tostring(now->tm_mon+1);
	h += "." + tostring(now->tm_year+1900);
	h += " at " + tostring(now->tm_hour);
	h += ":" + tostring(now->tm_min);
	h += ":" + tostring(now->tm_sec);
	return h;
}

void Header::init(std::string title){
	rst->title(title,"=");
	rst->textit(when());
	rst->np();
}

std::ostream& operator<<(std::ostream& flux, Header const& h){
	flux<<h.get();
	return flux;
}

void Header::hyperlink(std::string const& display, std::string const& link) {
	rst->hyperlink(display,link);
}
