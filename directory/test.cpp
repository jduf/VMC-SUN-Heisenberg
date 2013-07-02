#include "Directory.hpp"

int main(){
	Directory d;
	d.search_file_ext(".bin","/home/jdufour/travail/cpp-dev/",true);
	d.print();
}

