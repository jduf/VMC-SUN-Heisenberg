#include "Directory.hpp"
#include "Read.hpp"
#include "Write.hpp"
#include "RST.hpp"
#include "RSTfile.hpp"

#include <cstdlib>

int main(){
	Directory d;
	d.search_file_ext(".jdbin","/home/jdufour/travail/cpp-dev/",false);
	std::string save_in("/home/jdufour/travail/cpp-dev/rst/output/");
	std::string file;
	std::string command("");
	d.print();
	for(unsigned int i(0); i<d.size();i++){
		file = d.get_path(i) + "/" + d.get_name(i) + d.get_ext(i);
		Read r(file);
		file = d.get_name(i) + ".rst";
		Write w(save_in+file);
		w<<r.get_header();

		file = d.get_name(i);
		command = "rst2html " + save_in + file + ".rst " + save_in + file + ".html";  
		system(command.c_str());
	}

	RSTfile rst("index");
	rst.title("RST","=");
	rst.title("Possibilities","-");
	rst.text("Here I display all the bin files");
	rst.np();
	rst.item("this is an item");
	rst.item("and here another");
	rst.np();
	rst.text("this link should work");
	rst.hyperlink("output/data-1","output/data-1.html");
	rst.text("or at least, i hope");
	rst.np();
	rst.title("List of all the files","<");

	std::string file_link;
	for(unsigned int i(0);i<d.size();i++){
		file_link = save_in + d.get_name(i) + ".html";
		file = d.get_path(i) + "/" + d.get_name(i);
		rst.hyperlink(file,file_link);
		rst.np();
	}
}
