/*!  @file jdwrite.cpp */
#include "Read.hpp"
#include "Write.hpp"
#include "Rewrite.hpp"
#include "Linux.hpp"

void read(std::string filename){
	Read r(filename);
	Write w("/tmp/tmp.txt");
	w<<r.get_header();
}

void rewrite(std::string filename){
		Read tmp("/tmp/tmp.txt");
		std::string n_header;
		tmp>>n_header;
		Rewrite rw(filename);
		rw.rewrite_header(n_header);
}

int main(int argc, char* argv[]){
	if(argc==2){
		read(argv[1]);
		Linux command;
		command("vim /tmp/tmp.txt");
		rewrite(argv[1]);
	} else {
		std::cerr<<"jdwrite : take exactly one input argument"<<std::endl;
		std::cerr<<"usage -> ./jdwrite filename"<<std::endl;
		std::cerr<<"note  -> the file should be given without the .jdbin extension"<<std::endl;
	}
}
