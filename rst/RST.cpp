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

RST::~RST() { 
	if(w){ 
		rst += RST_np;
		//for(unsigned int i(0); i<links.size();i++){
			//rst += links[i] + RST_nl;
		//}
		(*w)<<rst;
		std::string command("rst2html " + filename + ".rst " + filename + ".html");  
		int status(0);
		status = system(command.c_str());
		if(status != 0){std::cerr<<"RST : ~RST() : the command function called returns an error"<<std::endl; }
		delete w;
	}
}

void RST::title(std::string t,std::string symb){
	rst += t + RST_nl;
	for(unsigned int i(0);i<t.size();i++){ rst +=  symb; }
	rst += RST_np;
}

void RST::text(std::string t){
	rst += t + RST_nl;
}

void RST::lineblock(std::string t){
	size_t pos0(0);
	size_t pos1(t.find("\n",pos0));
	while(t.find("\n",pos1+1)  != std::string::npos){
		rst += "| " +  t.substr(pos0,pos1-pos0) + RST_nl;
		pos0 = pos1+1;
		pos1 = t.find("\n",pos0);
	}
	rst += RST_nl;
}

void RST::textit(std::string t){
	 rst += "*" + t + "*" + RST_nl;
}

void RST::textbf(std::string t){
	 rst += "**" + t + "**" + RST_nl;
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
	rst += "`" + display + " <" + link + ">`_ " + RST_nl;
	//links.push_back(".. _" + t + ": " + l);
}

std::ostream& operator<<(std::ostream& flux, RST const& rst){
	flux<<rst.get();
	return flux;
}

