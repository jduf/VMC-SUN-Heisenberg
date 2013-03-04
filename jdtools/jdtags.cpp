#include "Read.hpp"
#include "RST.hpp"
#include "Array2D.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <cstdlib> // system(std::string commad)

void search_tag(std::string tag, Array2D<std::string> const& DF, std::vector<bool>& is_tag){
	for(unsigned int i(0);i<DF.row();i++){
		if(DF(i,0).find(tag) == std::string::npos){
			is_tag[i] = false;
		}
	}
}

int main(int argc, char* argv[]){
	if(argc > 1){
		std::string list_tags("");
		std::vector<std::string> list_file;
		Array2D<std::string> DF;
		std::string save_in("/home/jdufour/travail/cpp-dev/jdtools/rst-output/");
		Read r(save_in + "TAGS.bin");
		r>>DF;

		std::vector<bool> is_tag(DF.row(),true);
		for(int i(1);i<argc;i++){
			search_tag(argv[i],DF,is_tag);
		}
		for(int i(1);i<argc-1;i++){
			list_tags += argv[i];
			list_tags += ", ";
		}
		list_tags += argv[argc-1];

		RST rst(save_in + "TAGS");
		rst.title("File with tag *"+list_tags+"*","-");
		for(unsigned int i(0);i<DF.row();i++){
			if(is_tag[i]){
				rst.np();
				Read r(DF(i,1) + ".rst");
				std::string file_tagged;
				r>>file_tagged;
				rst.text(file_tagged);
			}
		}

		std::string command("");
		command = "firefox " +save_in + "TAGS.html";
		system(command.c_str());
	} else {
		std::cerr<<"jdtags : take at least one input argument"<<std::endl;
		std::cerr<<"usage -> ./jdtag tag1 tag2 ..."<<std::endl;
	}
}	
