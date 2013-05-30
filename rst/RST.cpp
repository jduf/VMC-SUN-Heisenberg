#include "RST.hpp"
#include "Write.hpp"

RST::RST():
	RST_nl("\n"),
	RST_np("\n\n"),
	RST_item("+ "),
	rst(""),
	path(""),
	filename(""),
	w(NULL),
	create_pdf(false)
{}

RST::RST(std::string path, std::string filename):
	RST_nl("\n"),
	RST_np("\n\n"),
	RST_item("+ "),
	rst(""),
	path(path),
	filename(filename),
	w(new Write(path + filename + ".rst")),
	create_pdf(false)
{ }

RST::~RST() { 
	if(w){ 
		rst += RST_np;
		//for(unsigned int i(0); i<links.size();i++){
		//rst += links[i] + RST_nl;
		//}
		(*w)<<rst;
		std::string command("rst2html " + path  + filename + ".rst " + path + filename + ".html");  
		int status(0);
		status = system(command.c_str());
		if(status != 0){
			std::cerr<<"RST : ~RST() : the command function called returns an error for the html creation"<<std::endl; 
		} else {
			if(create_pdf){
				command="rst2latex " + path + filename + ".rst " + path + filename + ".tex";  
				status = system(command.c_str());
				command = "pdflatex -output-directory " + path + " "  + filename + ".tex"; 
				status = system(command.c_str());
				if(status != 0){std::cerr<<"RST : ~RST() : the command function called returns an error for the pdf creation"<<std::endl; }
			}
		}
		delete w;
	}
}

void RST::pdf(){
	create_pdf=true;
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
		std::cerr<<"RST : def(string t,string def) : t too long"<<std::endl;
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

void RST::figure(std::string image, std::string legend, unsigned int scale){
	rst += ".. figure:: " + image + RST_nl;
	rst += "   :scale: " + tostring(scale) + " %" + RST_nl;
	rst += "   :align: center" + RST_np;
	rst += "   " + legend + RST_np;
}

std::ostream& operator<<(std::ostream& flux, RST const& rst){
	flux<<rst.get();
	return flux;
}

