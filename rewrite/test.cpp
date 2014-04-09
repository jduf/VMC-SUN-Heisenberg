#include "Read.hpp"
#include "Write.hpp"
#include "Rewrite.hpp"
#include "Linux.hpp"

#include <cstdlib>

int main(){
	std::string filename("data-5.jdbin");
	Read r(filename);
	Write w("/tmp/tmp.txt");
	w<<r.get_header();
	std::cout<<r.get_header()<<std::endl;

	Linux command;
	command("vim /tmp/tmp.txt");

	Read tmp("/tmp/tmp.txt");
	std::string n_header;
	tmp>>n_header;

	Rewrite rw(filename);
	rw.rewrite_header(n_header);
}
