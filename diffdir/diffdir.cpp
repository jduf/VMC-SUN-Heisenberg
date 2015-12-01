/*! @file diffdir.cpp*/

#include "Directory.hpp"
#include "Parseur.hpp"
#include "Linux.hpp"

#include <iostream>
#include <string>

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string directory_name_1(P.get<std::string>(0));
	std::string directory_name_2(P.get<std::string>(1));
	unsigned int i(0);
	std::string ext("");
	std::string keyword("");
	std::string host("");
	std::string dft_answer(P.find("y",i,false)?"y":"n");
	std::string dft_string(dft_answer=="y"?" ? [Y/n]:":" ? [y/N]:");
	bool ssh(P.find("ssh",i,false));
	if(ssh){ 
		std::cerr<<"diffdir : will only take into account the local directory"<<std::endl;
		host = P.get<std::string>(i); 
	}

	if(!P.locked()){
		Linux command;
		std::string answer("");

		if(directory_name_1[directory_name_1.size()-1] != '/'){ directory_name_1 += "/"; }
		if(directory_name_2[directory_name_2.size()-1] != '/'){ directory_name_2 += "/"; }
		Directory d1;
		Directory d2;
		if(P.find("e",i,false)){
			ext = P.get<std::string>(i);
			d1.search_file_ext(ext,directory_name_1,true,false);
			d2.search_file_ext(ext,directory_name_2,true,false);
		}
		if(P.find("k",i,false)){
			keyword = P.get<std::string>(i);
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
		while(i<d1.size()){
			if(d1.get_name(i) == d2.get_name(j)){
				for(unsigned int k(old_j+1); k<j; k++){
					std::cout<<"     B : "<< d2.get_name(k) << d2.get_ext(k)<<std::endl;
				}
				if(ssh){ command("ssh " + host + " cat '"+ d2[j] + "' | diff -q " + d1[i] + " - > /dev/null",false); } 
				else {   command("diff -q " + d1[i] + " " + d2[j]+ "> /dev/null",false); }
				if(command.status()){
					std::cout<<"A != B : "<< d1.get_name(i) << d1.get_ext(i) << dft_string;
					std::getline(std::cin,answer);
					if(answer==""){ answer = dft_answer; }
					if(answer=="y"){ 
						if(ssh){ command("vimdiff " + d1[i] + " scp://" + host + "/" + d2[j],false); } 
						else {   command("vimdiff " + d1[i] + " " + d2[j],false); }
					}
				} else { std::cout<<"A == B : "<<d1.get_name(i)<<d1.get_ext(i)<<std::endl; }
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
