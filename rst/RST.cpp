#include "RST.hpp"

/*{static element*/
std::string const RST::nl_ = "\n";
std::string const RST::np_ = "\n\n";
std::string const RST::item_ = "+ ";

std::string RST::math(std::string const& t){ return ":math:`" + t + "`"; }
std::string RST::textit(std::string const& t){ return "*" + t + "*"; }
std::string RST::textbf(std::string const& t){ return "**" + t + "**"; }
/*}*/

RST::RST(std::string const& rst):
	rst_(rst)
{}

void RST::title(std::string const& t, char const& symb){
	rst_ += RST::np_ + t + RST::nl_;
	for(unsigned int i(0);i<t.size();i++){ rst_ +=  symb; }
	rst_ += RST::np_;
}

void RST::text(std::string const& t){
	rst_ += t + RST::nl_;
}

void RST::math_centered(std::string const& t){
	rst_ += RST::nl_;
	rst_ += ".. math::" + RST::nl_;
	rst_ += " " + t + RST::np_;
}

void RST::item(std::string const& t){
	rst_ += RST::item_ + t + RST::nl_;
}

void RST::lineblock(std::string const& t){
	size_t pos0(0);
	size_t pos1(t.find("\n",pos0));
	while(t.find("\n",pos1+1) != std::string::npos){
		rst_ += "| " +  t.substr(pos0,pos1-pos0) + RST::nl_;
		pos0 = pos1+1;
		pos1 = t.find("\n",pos0);
	}
	rst_ += RST::nl_;
}

void RST::def(std::string const& t, std::string const& def){
	//if(def.size()>25){std::cerr<<"RST : def(t,def) : def.size>25"<<std::endl;}
	rst_ += ":" + t + ":" + " " + def + RST::nl_;
}

void RST::hyperlink(std::string const& display, std::string const& link){
	rst_ += "`" + display + " <" + link + ">`_ " + RST::nl_;
}

void RST::figure(std::string const& image, std::string const& legend, unsigned int width){
	rst_ += RST::nl_;
	rst_ += ".. figure:: " + image + RST::nl_;
	rst_ += "   :width: " + my::tostring(width)  + RST::nl_;
	rst_ += "   :align: center" + RST::np_;
	rst_ += "   " + legend + RST::np_;
}

void RST::link_figure(std::string const& image, std::string const& legend, std::string const& link, unsigned int width){
	rst_ += RST::nl_;
	rst_ += ".. figure:: " + image + RST::nl_;
	rst_ += "   :width: " + my::tostring(width)  + RST::nl_;
	rst_ += "   :align: center" + RST::nl_;
	rst_ += "   :target: " + link + RST::np_;
	rst_ += "   " + legend + RST::np_;
}

void RST::comment(std::string const& t){
	rst_ += ".. " + t + RST::np_;
}

void RST::np(){
	rst_ += RST::np_; 
}

void RST::nl(){
	rst_ += RST::nl_; 
}

std::ostream& operator<<(std::ostream& flux, RST const& rst){
	flux<<rst.get();
	return flux;
}
