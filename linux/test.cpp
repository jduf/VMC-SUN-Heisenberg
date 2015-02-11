#include "Linux.hpp" 
#include <iostream>

int main(){
	Linux command;
	std::string dir(command.pwd());
	std::cout<<dir<<std::endl;
	command("ls");
	command("mpg123 /home/media/miuz/s/*.mp3");
}
