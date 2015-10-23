#include "Linux.hpp" 
#include "Time.hpp" 
#include <iostream>
#include <vector>

int main(){
	//std::vector<double> a(1e9,1);
	std::string path1("aa/bb/cc/kk/ff/ee/");
	Linux command;
	unsigned int N(1e2);
	Time chrono;
	chrono.set();
	for(unsigned int i(0);i<N;i++){
		command.mkdir("aa/bb/cc/kk/ff/ee/");
	}
	std::cout<<chrono.elapsed()<<std::endl;
	chrono.set();
	for(unsigned int i(0);i<N;i++){
		system("mkdir -p AA/BB/CC/KK/FF/EE/"); 
	}
	std::cout<<chrono.elapsed()<<std::endl;
	command.mkdir("AA");

	command.open("gp1.log");

	command("gnuplot -e 'plot exp(x)'",false);
	std::cout<<command.status()<<std::endl;
	Linux command_bis;

	command_bis.open("gp1.log");
	command.close(false);
	command_bis.open("gp2.log");

	command_bis("gnuplot -p -e 'plot cos(x)'",false);
	command_bis.close(true);
}
