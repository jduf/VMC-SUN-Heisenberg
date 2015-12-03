#include "Directory.hpp"

int main(){
	Directory d;
	d.search_file_ext(".cpp","./",true,true);
	d.print();
}
