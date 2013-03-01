#include "Read.hpp"
#include "Write.hpp"
#include "Rewrite.hpp"
#include <cstdlib>


int main(){
	std::string filename("data-5.jdbin");
	Read r(filename,true);
	Write w("/tmp/tmp.txt");
	w<<r.header();
	std::cout<<r.header()<<std::endl;

	system("vim /tmp/tmp.txt");

	Read tmp("/tmp/tmp.txt");
	std::string n_header;
	tmp>>n_header;

	Rewrite rw(filename);
	rw.rewrite_header(n_header);
}
