#include "RST.hpp"
#include "Write.hpp"

RST::RST():
	RST_nl("\n"),
	RST_np("\n\n"),
	RST_title("="),
	RST_subtitle("-"),
	RST_item("+ "),
	rst(""),
	links(0),
	filename(""),
	w(NULL)
{}

RST::RST(std::string filename):
	RST_nl("\n"),
	RST_np("\n\n"),
	RST_title("="),
	RST_subtitle("-"),
	RST_item("+ "),
	rst(""),
	links(0),
	filename(filename),
	w(new Write(filename + ".rst"))
{ }

RST::~RST()
{ 
	if(w){ 
		rst += RST_np;
		for(unsigned int i(0); i<links.size();i++){
			rst += links[i] + RST_nl;
		}
		(*w)<<rst;
		std::string command("rst2html " + filename + ".rst " + filename + ".html");  
		system(command.c_str());
		delete w;
	}
}

void RST::title(std::string t){
	rst += RST_np + t + RST_nl;
	for(unsigned int i(0);i<t.size();i++){ rst +=  RST_title; }
	rst += RST_np;
}

void RST::subtitle(std::string t){
	rst += RST_np + t + RST_nl;
	for(unsigned int i(0);i<t.size();i++){ rst += RST_subtitle; }
	rst += RST_np;
}

void RST::text(std::string t){
	 rst += t + " ";
}

void RST::textit(std::string t){
	 rst += " *" + t + "* ";
}

void RST::textbf(std::string t){
	 rst += " **" + t + "** ";
}

void RST::item(std::string t){
	 rst += RST_item + t + RST_nl;
}

void RST::np(){
	 rst += RST_np;
}

void RST::hyperlink(std::string t, std::string l){
	 rst += "`" + t + "`_ ";
	links.push_back(".. _" + t + ": " + l);
}

std::ostream& operator<<(std::ostream& flux, RST const& rst){
	flux<<rst.get();
	return flux;
}

