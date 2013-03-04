#include "Read.hpp"
#include "Write.hpp"
#include "Rewrite.hpp"
#include <cstdlib>

int main(int argc, char* argv[]){
	if(argc==2){
		std::string filename(argv[1]);
		Read r(filename);
		Write w("/tmp/tmp.txt");
		w<<r.header();
		std::cout<<r.header()<<std::endl;

		system("vim /tmp/tmp.txt");

		Read tmp("/tmp/tmp.txt");
		std::string n_header;
		tmp>>n_header;

		Rewrite rw(filename);
		rw.rewrite_header(n_header);
	} else {
		std::cerr<<"jdwrite : take exactly one input argument"<<std::endl;
		std::cerr<<"usage -> ./jdwrite filename"<<std::endl;
		std::cerr<<"note  -> the file should be given without the .jdbin extension"<<std::endl;
	}
}
