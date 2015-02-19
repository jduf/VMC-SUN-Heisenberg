#include "Matrix.hpp"
#include "Container.hpp"

int main(){
	Container c;
	c.set("int",1);
	c.set("string",std::string("string"));
	c.set("double",1.2);
	std::cout<<c.get<unsigned int>("int")<<std::endl;
	std::cout<<c.get<std::string>("string")<<std::endl;
	double d;
	c.get("double",d);
	std::cout<<d<<std::endl;
}
