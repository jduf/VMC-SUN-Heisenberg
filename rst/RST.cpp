#include "RST.hpp"

RST::RST():
	RST_nl_("\n"),
	RST_np_("\n\n"),
	RST_item_("+ "),
	rst_("")
{}

RST::RST(std::string rst):
	RST_nl_("\n"),
	RST_np_("\n\n"),
	RST_item_("+ "),
	rst_(rst)
{}

void RST::title(std::string t,std::string symb){
	rst_ += RST_np_ + t + RST_nl_;
	for(unsigned int i(0);i<t.size();i++){ rst_ +=  symb; }
	rst_ += RST_np_;
}

void RST::text(std::string t){
	rst_ += t + RST_nl_;
}

void RST::lineblock(std::string t){
	size_t pos0(0);
	size_t pos1(t.find("\n",pos0));
	while(t.find("\n",pos1+1)  != std::string::npos){
		rst_ += "| " +  t.substr(pos0,pos1-pos0) + RST_nl_;
		pos0 = pos1+1;
		pos1 = t.find("\n",pos0);
	}
	rst_ += RST_nl_;
}

void RST::textit(std::string t){
	rst_ += "*" + t + "*" + RST_nl_;
}

void RST::textbf(std::string t){
	rst_ += "**" + t + "**" + RST_nl_;
}

void RST::item(std::string t){
	rst_ += RST_item_ + t + RST_nl_;
}

void RST::def(std::string t, std::string def){
	//if(def.size()>25){std::cerr<<"RST : def(t,def) : def.size>25"<<std::endl;}
	rst_ += ":" + t + ":" + " " + def + RST_nl_;
}

void RST::np(){
	rst_ += RST_np_;
}

void RST::nl(){
	rst_ += RST_nl_;
}

void RST::hyperlink(std::string display, std::string link){
	rst_ += "`" + display + " <" + link + ">`_ " + RST_nl_;
	//links.push_back(".. _" + t + ": " + l);
}

void RST::figure(std::string image, std::string legend, unsigned int width){
	rst_ += RST_nl_;
	rst_ += ".. figure:: " + image + RST_nl_;
	rst_ += "   :width: " + tostring(width)  + RST_nl_;
	rst_ += "   :align: center" + RST_np_;
	rst_ += "   " + legend + RST_np_;
}

void RST::link_figure(std::string image, std::string legend, std::string link, unsigned int width){
	rst_ += RST_nl_;
	rst_ += ".. figure:: " + image + RST_nl_;
	rst_ += "   :width: " + tostring(width)  + RST_nl_;
	rst_ += "   :align: center" + RST_nl_;
	rst_ += "   :target: " + link + RST_np_;
	rst_ += "   " + legend + RST_np_;
}

std::ostream& operator<<(std::ostream& flux, RST const& rst){
	flux<<rst.get();
	return flux;
}
