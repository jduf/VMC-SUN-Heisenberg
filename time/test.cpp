#include "Time.hpp"

#include<iostream>
#include<unistd.h>

int main(){
	Time t;
	std::cout<<t.day()<<" "<<t.month()<<" "<<t.year()<<" "<<t.hour()<<" "<<t.min()<<" "<<t.sec()<<std::endl;
	std::string d;
	std::cout<<t.date()<<std::endl;
	unsigned int i(0);
	do{
		std::cout<<i<<std::endl;
		sleep(3);
		i++;
	} while (!t.limit_reached(5));
}

