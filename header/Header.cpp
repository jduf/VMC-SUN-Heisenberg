#include "Header.hpp"

Header::Header(std::string const& s):
	h("")
{ 
	when();
	intro(s);
}

Header::Header(Read& r):
	h("")
{
	r >> h;
}

void Header::add(std::string str, double const& d){
	h += rst.item + "double " + str + "=" + tostring(d) + "\n";
}

void Header::add(std::string str, std::complex<double> const& c){
	h += rst.item + "complex " + str + "=" +tostring(c) + "\n";
}

void Header::add(std::string str, Matrice<double>& m){
	h += rst.item + "Matrice<double> " + str + " of size " + tostring(m.size()) + "x" + tostring(m.size()) + "\n";
}

void Header::write(Write& w){
	std::string t(w.get_filename());
	h.insert(0,title(t));
	w<<h;
}	

void Header::when(){
	time_t t(time(0));
	struct tm* now(localtime(&t));
	h += rst.it + "Created on";
	h += " " + tostring(now->tm_mday);
	h += "." + tostring(now->tm_mon+1);
	h += "." + tostring(now->tm_year+1900);
	h += " at " + tostring(now->tm_hour);
	h += ":" + tostring(now->tm_min);
	h += ":" + tostring(now->tm_sec) + rst.it + rst.np;
}

void Header::intro(std::string const& s){
	h += rst.bf + s + rst.bf + rst.np;
}

std::string Header::title(std::string& t){
	unsigned int N(t.size());
	t += "\n";
	t.append(N,rst.title);
	return t+rst.np;
}
