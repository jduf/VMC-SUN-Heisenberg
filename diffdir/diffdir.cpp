#include "Directory.hpp"
#include "Parseur.hpp"
#include "Linux.hpp"

#include <iostream>
#include <string>

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string directory_name_1(P.get<std::string>("0"));
	std::string directory_name_2(P.get<std::string>("1"));
	std::string ext(P.get<std::string>("e"));
	if(!P.status()){
		Linux command;
		std::string answer("");

		if(directory_name_1[directory_name_1.size()-1] != '/'){
			directory_name_1 += "/";
		}
		if(directory_name_2[directory_name_2.size()-1] != '/'){
			directory_name_2 += "/";
		}
		Directory d1;
		d1.search_file_ext(ext,directory_name_1,true);
		d1.sort();
		Directory d2;
		d2.search_file_ext(ext,directory_name_2,true);
		d2.sort();
		
		//bool match(false);
		//for(unsigned int i(0);i<d1.size();i++){
			//for(unsigned int j(0);j<d2.size();j++){
				//if(d1.get_name(i)==d2.get_name(j)){
					//command("diff -q " + d1[i] + " " + d2[j]);
					//if(command.status()){
						//std::cout<<"Do you want to compare " << d1.get_name(i) << d1.get_ext(i) << " ? [y/N] ";
						//std::getline(std::cin,answer);
						//if(answer=="y"){	
							//command("vimdiff " + d1[i] + " " + d2[j]);
						//}
					//}
					//else { std::cerr<<"the file "<<d1[i]<<" and "<<d2[i]<<" are identical"<<std::endl; }
					//j=d2.size();
					//match = true;
				//}
			//}
			//if(match){ match=false; }
			//else { std::cerr<<"----> the file "<<d1[i]<<" doesn't exist in "<< directory_name_2<<std::endl; }
		//}
		std::cout<<"A : "<< directory_name_1<<std::endl;
		std::cout<<"B : "<< directory_name_2<<std::endl;
		unsigned int i(0);
		unsigned int j(0);
		unsigned int old_j(0);
		do {
			if(j >= d2.size()) {
				j = old_j+1;
				std::cout<<"A      : "<< d1.get_name(i) << d1.get_ext(i)<<std::endl;
				i++;
			}
			if(d1.get_name(i) == d2.get_name(j)){
				for(unsigned int k(old_j+1); k<j; k++){
					std::cout<<"     B : "<< d2.get_name(k) << d2.get_ext(k)<<std::endl;
				}
				old_j = j;
				command("diff -q " + d1[i] + " " + d2[j] + " >> /tmp/diff.txt");
				if(command.status()){
					std::cout<<"A != B : "<< d1.get_name(i) << d1.get_ext(i) << " ? [y/N] ";
					std::getline(std::cin,answer);
					if(answer=="y"){	
						command("vimdiff " + d1[i] + " " + d2[j]);
					}
				} else {
					std::cout<<"A == B : "<<d1.get_name(i)<<d1.get_ext(i)<<std::endl;
				}
				i++;
			}
			j++;
		} while (i<d1.size());
	} else {
		std::cerr<<"diffdir : You need to precise the extension and the two directories that you want to compare"<<std::endl;
	}
}
