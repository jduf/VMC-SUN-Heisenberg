/*! @file load.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);

	unsigned int i;
	if(P.find("sim",i,false)){
		IOFiles read(P.get<std::string>(i),false);

		Vector<double> tmp(read);
		System sys(read);
		CreateSystem cs(&sys);
		cs.init(&tmp,NULL);

		RSTFile rst("/tmp/",cs.get_filename());
		IOSystem ios(cs.get_filename(),"","","","","/tmp/",&rst);
		cs.set_IOSystem(&ios);

		cs.display_results();

		rst.text(read.get_header());
		rst.save(false,true);
		Linux command;
		command(Linux::html_browser("/tmp/"+cs.get_filename()+".html"),true);
	}
}
