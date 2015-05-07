/*!  @file jdread.cpp */
#include "RSTFile.hpp"

int main(int argc, char* argv[]){
	if(argc == 2){
		std::string name(argv[1]);
		IOFiles r(name,false);
		name = r.get_filename();
		name = name.substr(name.find_last_of('/'),name.size()-6);
		RSTFile html("/tmp/",name);
		html.text(r.get_header());
		html.save(false);
		Linux command;
		command("firefox /tmp/"+name+".html");
		command("rm /tmp/"+name+"*");
	} else {
		std::cerr<<"jdhtml : takes exactly one input argument : ./jdhtml filename"<<std::endl;
	}
}	
