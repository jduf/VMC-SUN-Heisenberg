#include "RST.hpp"

RST::RST():
	RST_nl("\n"),
	RST_np("\n\n"),
	RST_item("+ "),
	rst("")
{}

RST::RST(std::string rst):
	RST_nl("\n"),
	RST_np("\n\n"),
	RST_item("+ "),
	rst(rst)
{}

void RST::title(std::string t,std::string symb){
	rst += RST_np + t + RST_nl;
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
	if(def.size()>25){std::cerr<<"RST : def(t,def) : def.size>25"<<std::endl;}
	rst += ":" + t + ":" + " " + def + RST_nl;
}

void RST::np(){
	rst += RST_np;
}

void RST::nl(){
	rst += RST_nl;
}

void RST::hyperlink(std::string display, std::string link){
	rst += "`" + display + " <" + link + ">`_ " + RST_nl;
	//links.push_back(".. _" + t + ": " + l);
}

void RST::figure(std::string image, std::string legend, unsigned int width){
	rst += RST_nl;
	rst += ".. figure:: " + image + RST_nl;
	rst += "   :width: " + tostring(width)  + RST_nl;
	rst += "   :align: center" + RST_np;
	rst += "   " + legend + RST_np;
}

void RST::link_figure(std::string image, std::string legend, std::string link, unsigned int width){
	rst += RST_nl;
	rst += ".. figure:: " + image + RST_nl;
	rst += "   :width: " + tostring(width)  + RST_nl;
	rst += "   :align: center" + RST_nl;
	rst += "   :target: " + link + RST_np;
	rst += "   " + legend + RST_np;
}

std::ostream& operator<<(std::ostream& flux, RST const& rst){
	flux<<rst.get();
	return flux;
}
