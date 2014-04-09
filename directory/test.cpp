#include "Directory.hpp"

int main(){
	Directory d;
	d.search_file_ext(".dat","/home/jdufour/travail/cpp-dev/",true,true);
	d.print();
}

