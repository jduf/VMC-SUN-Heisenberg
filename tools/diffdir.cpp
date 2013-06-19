#include "Directory.hpp"
#include "Parseur.hpp"

#include <iostream>
#include <string>
#include <cstdlib> // system(std::string commad)

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string directory_name_1(P.get<std::string>("0"));
	std::string directory_name_2(P.get<std::string>("1"));
	std::string ext(P.get<std::string>("e"));
	if(!P.status()){
		std::string command("");
		std::string answer("");
		bool match(false);

		if(directory_name_1[directory_name_1.size()-1] != '/'){
			directory_name_1 += "/";
		}
		if(directory_name_2[directory_name_2.size()-1] != '/'){
			directory_name_2 += "/";
		}
		Directory d1;
		d1.search_file_ext(ext,directory_name_1,true);
		d1.sort();
		d1.print();
		Directory d2;
		d2.search_file_ext(ext,directory_name_2,true);
		d2.sort();
		d2.print();
		
		for(unsigned int i(0);i<d1.size();i++){
			for(unsigned int j(0);j<d2.size();j++){
				if(d1.get_name(i)==d2.get_name(j)){
					std::cout<<"Do you want to compare " << d1.get_name(i) << d1.get_ext(i) << " ? [y/N] ";
					std::getline(std::cin,answer);
					if(answer=="y"){	
						command = "vimdiff " + d1[i] + " " + d2[j];
						system(command.c_str());
					}
					match = true;
				}
			}
			if(match){ match=false; }
			else { std::cerr<<"the file "<<d1[i]<<" doesn't exist in "<< directory_name_2<<std::endl; }
		}
	} else {
		std::cerr<<"diffdir : You need to precise the extension and the two directories that you want to compare"<<std::endl;
	}
}
