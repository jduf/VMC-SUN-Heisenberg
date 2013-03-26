#include "RST.hpp"
#include "Write.hpp"

RST::RST():
	RST_nl("\n"),
	RST_np("\n\n"),
	RST_item("+ "),
	rst(""),
	filename(""),
	w(NULL)
{}

RST::RST(std::string filename):
	RST_nl("\n"),
	RST_np("\n\n"),
	RST_item("+ "),
	rst(""),
	filename(filename),
	w(new Write(filename + ".rst"))
{ }

RST::~RST()
{ 
	if(w){ 
		rst += RST_np;
		//for(unsigned int i(0); i<links.size();i++){
			//rst += links[i] + RST_nl;
		//}
		(*w)<<rst;
		std::string command("rst2html " + filename + ".rst " + filename + ".html");  
		system(command.c_str());
		delete w;
	}
}

void RST::title(std::string t,std::string symb){
	rst += RST_nl + t + RST_nl;
	for(unsigned int i(0);i<t.size();i++){ rst +=  symb; }
	rst += RST_np;
}

void RST::text(std::string t){
	 rst += t + RST_nl;
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

void RST::def(std::string t, std::string def){
	if(def.size()>25){
		std::cerr<<"RST : too long"<<std::endl;
	}
	rst += ":" + t + ":" + " " + def + RST_nl;
}

void RST::np(){
	 rst += RST_np;
}

void RST::hyperlink(std::string display, std::string link){
	rst += "`" + display + " <" + link + ">`_ ";
	//links.push_back(".. _" + t + ": " + l);
}

std::ostream& operator<<(std::ostream& flux, RST const& rst){
	flux<<rst.get();
	return flux;
}

