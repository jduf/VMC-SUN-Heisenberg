/*!  @file jdread.cpp */
#include "RSTFile.hpp"

int main(int argc, char* argv[]){
	if(argc == 2){
		std::string name(argv[1]);
		IOFiles r(name,false,false);
		name = r.get_filename();
		std::size_t pos(name.find_last_of('/'));
		if(pos != std::string::npos){ name = name.substr(pos,name.size()-6); }
		RSTFile html("/tmp/",name);
		html.text(r.get_header());
		html.save(false,false);
		Linux command;
		command(Linux::html_browser("/tmp/"+name+".html"),true);
		command("rm /tmp/"+name+"*",false);
	} else {
		std::cerr<<"jdhtml : takes exactly one input argument : ./jdhtml filename"<<std::endl;
	}
}	
