/*! @file diffdir.cpp*/

#include "Directory.hpp"
#include "Parseur.hpp"
#include "Linux.hpp"

#include <iostream>
#include <string>

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string directory_name_1(P.get<std::string>("0"));
	std::string directory_name_2(P.get<std::string>("1"));
	std::string ext("");
	std::string keyword("");
	std::string scp("");
	unsigned int i(0);
	if(P.search("e",i)){ ext = P.get<std::string>(i); }
	if(P.search("k",i)){ keyword = P.get<std::string>(i); }
	if(P.search("scp",i)){ scp = P.get<std::string>(i); }

	if(!P.status()){
		Linux command;
		std::string answer("");

		if(directory_name_1[directory_name_1.size()-1] != '/'){ directory_name_1 += "/"; }
		if(directory_name_2[directory_name_2.size()-1] != '/'){ directory_name_2 += "/"; }
		Directory d1;
		Directory d2;
		if(ext!=""){
			d1.search_file_ext(ext,directory_name_1,true,false);
			d2.search_file_ext(ext,directory_name_2,true,false);
		}
		if(keyword!=""){
			d1.search_file(keyword,directory_name_1,true,false);
			d2.search_file(keyword,directory_name_2,true,false);
		}
		d1.sort();
		d2.sort();

		std::cout<<"A : "<< directory_name_1<<" contains "<<d1.size()<<" files matching '"<<ext+keyword<<"'"<<std::endl;
		std::cout<<"B : "<< directory_name_2<<" contains "<<d2.size()<<" files matching '"<<ext+keyword<<"'"<<std::endl;
		unsigned int i(0);
		unsigned int j(0);
		unsigned int old_j(0);
		while (i<d1.size()) {
			if(d1.get_name(i) == d2.get_name(j)){
				for(unsigned int k(old_j+1); k<j; k++){
					std::cout<<"     B : "<< d2.get_name(k) << d2.get_ext(k)<<std::endl;
				}
				if(scp == ""){ 
					command("diff -q " + d1[i] + " " + d2[j] + " >> /tmp/diff.txt");
					if(command.status()){
						std::cout<<"A != B : "<< d1.get_name(i) << d1.get_ext(i) << " ? [y/N] ";
						std::getline(std::cin,answer);
						if(answer=="y"){ command("vimdiff " + d1[i] + " " + d2[j]); }
					} else {
						std::cout<<"A == B : "<<d1.get_name(i)<<d1.get_ext(i)<<std::endl;
					}
				} else { 
					command("ssh "+ scp +" 'cat " + d2[j] + "' | diff -q - " + d1[i] + " >> /tmp/diff.txt");
					if(command.status()){
						std::cout<<"A != B : "<< d1.get_name(i) << d1.get_ext(i) << " ? [y/N] ";
						std::getline(std::cin,answer);
						if(answer=="y"){ command("vimdiff " + d1[i] + " scp://"+scp+"/"+ d2[j]);  }
					} else {
						std::cout<<"A == B : "<<d1.get_name(i)<<d1.get_ext(i)<<std::endl;
					}
				}
				old_j = j;
				i++;
			}
			j++;
			if(j >= d2.size() && i < d1.size()) {
				std::cout<<"A      : "<< d1.get_name(i) << d1.get_ext(i)<<std::endl;
				j=old_j;
				i++;
			}
		}
		for(unsigned int k(old_j+1); k<d2.size(); k++){
			std::cout<<"     B : "<< d2.get_name(k) << d2.get_ext(k)<<std::endl;
		}
	} else {
		std::cerr<<"diffdir : You need to precise the extension and the two directories that you want to compare"<<std::endl;
	}
}
