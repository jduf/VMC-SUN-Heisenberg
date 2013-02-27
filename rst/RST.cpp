#include "RST.hpp"

RST::RST(std::string fname):
	RST_nl("\n"),
	RST_np("\n\n"),
	RST_title("="),
	RST_subtitle("-"),
	RST_item("+ "),
	w(fname + ".rst"),
	fname(fname),
	links(0)
{}

RST::~RST(){
	w<<RST_np;
	for(unsigned int i(0); i<links.size();i++){
		w<<links[i]<<RST_nl;
	}
	std::string command("rst2html " + fname + ".rst " + fname + ".html");  
	system(command.c_str());
}

void RST::title(std::string t){
	w<<RST_np<<t<<RST_nl;
	for(unsigned int i(0);i<t.size();i++){ w << RST_title; }
	w<<RST_np;
}

void RST::subtitle(std::string t){
	w<<RST_np<<t<<RST_nl;
	for(unsigned int i(0);i<t.size();i++){ w << RST_subtitle; }
	w<<RST_np;
}

void RST::text(std::string t){
	w<<t<<" ";
}

void RST::item(std::string t){
	w<<RST_item<<t<<RST_nl;
}

void RST::np(){
	w<<RST_np;
}

void RST::hyperlink(std::string t, std::string l){
	w<<"`"<<t<<"`_ ";
	links.push_back(".. _" + t + ": " + l);
}

