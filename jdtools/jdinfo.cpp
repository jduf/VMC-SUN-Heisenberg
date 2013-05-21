#include "Directory.hpp"
#include "Parseur.hpp"
#include "Read.hpp"
#include "Write.hpp"
#include "RST.hpp"

#include <string>
#include <vector>
#include <cstdlib> // system(std::string commad)


void update_readme(Directory const& d, std::string const& directory_name);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	if(argc == 2){
		char buff[PATH_MAX];
		getcwd(buff,PATH_MAX);
		std::string directory_name("");
		std::string save_in("");
		std::string command("");

		P.set(directory_name);
		if(directory_name == "."){
			directory_name = buff;
			if(directory_name[directory_name.size()-1] != '/'){
				directory_name = directory_name + "/";
			}
			save_in = directory_name + "info/";
		} else {
			directory_name ="/"+directory_name;
			directory_name = buff+directory_name;
			if(directory_name[directory_name.size()-1] != '/'){
				directory_name = directory_name + "/";
			}
			save_in = directory_name + "info/";
		}
		Directory d;
		d.search_file_ext(".png",directory_name);
		d.sort();
		
		update_readme(d,directory_name);

		command = "firefox " + save_in + "README.html &";
		system(command.c_str());
	} else {
		std::cerr<<"need to give a directory"<<std::endl;
	 }
}

void update_readme(Directory const& d, std::string const& directory_name){
	Read r(directory_name + "README");
	std::string h("");
	r>>h;

	RST rst(directory_name + "info/","README");
	rst.text(h);
	rst.hyperlink("List of all simulations", directory_name + "info/index.html");
	rst.np();
	rst.title("PLOTS","=");
	for(unsigned int i(0);i<d.size();i++){
		rst.figure(d[i],d[i],80);
	}
	rst.pdf();
}


