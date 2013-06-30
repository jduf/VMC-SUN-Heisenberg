/*!  @file jdtags.cpp */

#include "Linux.hpp"
#include "Read.hpp"
#include "RSTfile.hpp"
#include "Array2D.hpp"
#include "Parseur.hpp"

#include <dirent.h>// not useful but need to find the correct include to let getcwd work
#include <string>
#include <vector>
#include <iostream>

void search_tag(std::string save_in, int argc, char* argv[]);

int main(int argc, char* argv[]){
	if(argc > 1){
		Linux command;
		char buff[PATH_MAX];
		getcwd(buff,PATH_MAX);
		std::string save_in(std::string(buff) + "/info/");
		std::string list_tags("");

		search_tag(save_in,argc,argv);

		command("firefox " +save_in + "TAGS.html");
	} else {
		std::cerr<<"jdtags : take at least one input argument"<<std::endl;
		std::cerr<<"usage -> ./jdtag tag1 tag2 ..."<<std::endl;
	}
}	

void search_tag(std::string save_in, int argc, char* argv[]){
	std::string list_tags("");
	for(int i(1);i<argc-1;i++){
		list_tags += std::string(argv[i])+", ";
	}
	list_tags += std::string(argv[argc-1]);
	RSTfile rst("TAGS",save_in);
	rst.title("File with tag *"+list_tags+"*","+");

	Array2D<std::string> DF;
	Read r(save_in + "TAGS.bin");
	r>>DF;
	for(int i(1);i<argc;i++){
		for(unsigned int j(0);j<DF.row();j++){
			if(DF(j,0).find(argv[i]) != std::string::npos){
				Read r_tag(DF(j,1) + ".rst");
				std::string file_tagged;
				r_tag>>file_tagged;
				rst.text(file_tagged);
			}
		}
	}
}
