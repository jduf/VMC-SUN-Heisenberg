#include "List.hpp"

int main(){
	List<double> a;
	for(unsigned int i(0);i<10;i++){ a.append(i); }
	std::cout<<a<<std::endl;
	std::cout<<"size of the list "<<a.size()<<std::endl;
	std::cout<<"will remove 4th entry"<<std::endl;
	a.remove(4);
	std::cout<<a<<std::endl;
	std::cout<<"will remove last entry"<<std::endl;
	a.remove(a.size()-1);
	std::cout<<a<<std::endl;
	std::cout<<"will remove 1st entry"<<std::endl;
	a.remove(0);
	std::cout<<a<<std::endl;
	std::cout<<"new list with 1 entry"<<std::endl;
	List<double> b(0);
	std::cout<<b<<std::endl;
	std::cout<<"will remove 1ss entry"<<std::endl;
	b.remove(0);
	std::cout<<b<<std::endl;
	std::cout<<"new list copied from the first one between enty 3 and 5"<<std::endl;
	List<double> c(a.sublist(3,5));/*understand why it doesn't call a copy constructor*/
	std::cout<<c<<std::endl;
	std::cout<<"will remove the 2nd entry"<<std::endl;
	c.remove(1);
	std::cout<<c<<std::endl;
}
