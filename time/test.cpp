#include "Time.hpp"

#include<iostream>

int main(){
	Time t;
	std::cout<<t.day()<<" "<<t.month()<<" "<<t.year()<<" "<<t.hour()<<" "<<t.min()<<" "<<t.sec()<<std::endl;
}

