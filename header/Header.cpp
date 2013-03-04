#include "Header.hpp"

Header::Header():
	rst(new RST)
{ }

Header::~Header(){delete rst; }

void Header::add(std::string const& s, double const& d){
	rst->def(s + "=" + tostring(d), "double");
}

void Header::add(std::string const& s, std::complex<double> const& c){
	rst->def(s + "=" + tostring(c), "complex");
}

void Header::add(std::string const& s, Matrice<double> const& mat){
	rst->def(s + "(" + tostring(mat.size()) + "x" + tostring(mat.size())+ ")","Matrice<double>");
}

void Header::add(std::string const& s, Matrice<std::complex<double> > const& mat){
	rst->def(s + "(" + tostring(mat.size()) + "x" + tostring(mat.size())+ ")","Matrice<complex>");
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
