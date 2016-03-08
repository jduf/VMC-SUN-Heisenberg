#include "RST.hpp"

/*{static element*/
std::string const RST::nl_ = "\n";
std::string const RST::np_ = "\n\n";
std::string const RST::item_ = " + ";
std::string const RST::hashtag_line_ = std::string(35,'#');
std::string const RST::dash_line_ = std::string(35,'-');

std::string RST::textit(std::string const& t){ return " *" + t + "* "; }
std::string RST::textbf(std::string const& t){ return " **" + t + "** "; }
std::string RST::math(std::string const& t)  { return ":math:`" + t + "` "; }
std::string RST::width(std::string const& t) { return "   :width: " + t + RST::nl_; }
std::string RST::scale(std::string const& t) { return "   :scale: " + t + RST::nl_; }
std::string RST::target(std::string const& t){ return "   :target: "+ t + RST::nl_; }
/*}*/

void RST::title(std::string const& t, char const& symb){
	rst_ += RST::np_;
	rst_ += t + RST::nl_;
	rst_ += std::string(t.size(),symb) + RST::np_;
}

void RST::text(std::string const& t){
	rst_ += t + RST::nl_;
}

void RST::math_paragraph(std::string const& t){
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
	rst_ += ":" + t + ":" + " " + def + RST::nl_;
}

void RST::hyperlink(std::string const& display, std::string const& link){
	rst_ += "`" + display + " <" + link + ">`_ " + RST::nl_;
}

void RST::figure(std::string const& image, std::string const& legend, std::string const& option){
	rst_ += RST::nl_;
	rst_ += ".. figure:: " + image + RST::nl_;
	rst_ += option;
	rst_ += "   :align: center" + RST::np_;
	rst_ += "   " + legend + RST::np_;
}

void RST::comment(std::string const& t){
	rst_ += ".. " + t + RST::np_;
}

void RST::replace(std::string const& pattern, std::string const& replace){
	rst_ += RST::nl_;
	rst_ += ".. |" + pattern + "| replace:: " + replace + RST::np_;
}

void RST::change_text_onclick(std::string const& old_txt, std::string const& new_txt){
	rst_ += ".. raw:: html" + RST::np_;
	rst_ += "   <p style=\"color:green\" onclick=\"this.innerHTML = '" + new_txt + "'\">" + old_txt + "</p>" + RST::np_;
}

void RST::np(){ rst_ += RST::np_; }

void RST::nl(){ rst_ += RST::nl_; }

std::ostream& operator<<(std::ostream& flux, RST const& rst){
	flux<<rst.get();
	return flux;
}
