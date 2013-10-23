#include "Time.hpp"

#include<iostream>

int main(){
	Time t;
	std::cout<<t.day()<<" "<<t.month()<<" "<<t.year()<<" "<<t.hour()<<" "<<t.min()<<" "<<t.sec()<<std::endl;
	unsigned int i(0);
	do{
		std::cout<<i<<std::endl;
		sleep(3);
		i++;
	} while (!t.time_limit_reached(5));
}

