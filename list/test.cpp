#include "List.hpp"

class A{
	public:
		A(double a):a_(a){
			N_++;
			std::cout<<"normal"<<std::endl;
		}

		A(A const& a):a_(a.a_){
			N_++;
			std::cout<<"copy"<<std::endl;
		}

		A(A&& a):a_(std::move(a.a_)){
			std::cout<<"move"<<std::endl;
		}

		void print(std::ostream& flux) const { flux<<a_; }

		static unsigned int N_;
		double a_;
};

unsigned int A::N_ = 0;

std::ostream& operator<<(std::ostream& flux, A const& a){
	a.print(flux);
	return flux;
}

int main(){
	List<A> a;
	for(unsigned int i(0);i<10;i++){ a.add_end(A(i)); }
	std::cout<<a<<std::endl;
	std::cout<<A::N_<<std::endl;
	std::cout<<"size of the list "<<a.size()<<std::endl;
	A b(std::move(A(0)));
	std::cout<<A::N_<<std::endl;
	//std::cout<<"will remove 4th entry with pop(4)"<<std::endl;
	//a.pop(4);
	//std::cout<<a<<std::endl;
	//std::cout<<"will remove last entry with pop(idx)"<<std::endl;
	//a.pop(a.size()-1);
	//std::cout<<"will remove last entry with pop_end"<<std::endl;
	//a.pop_end();
	//std::cout<<a<<std::endl;
	//std::cout<<"will remove 1st entry with pop(0)"<<std::endl;
	//a.pop(0);
	//std::cout<<a<<std::endl;
	//std::cout<<"new list with 1 entry set to 3"<<std::endl;
	//List<double> b(3);
	//std::cout<<b<<std::endl;
	//std::cout<<"will remove 1st entr with pop_start()"<<std::endl;
	//b.pop_start();
	//std::cout<<b<<std::endl;
	//std::cout<<"new list copied from the first one between enty 3 and 5"<<std::endl;
	//List<double> c(a.sublist(3,5));/*understand why it doesn't call a copy constructor*/
	//std::cout<<c<<std::endl;
	//std::cout<<"will remove the 2nd entry"<<std::endl;
	//c.pop(1);
	//std::cout<<c<<std::endl;
	//std::cout<<"will swap the 2nd and 5th entry of the first list"<<std::endl;
	//a.swap(2,5);
	//std::cout<<a<<std::endl;
	//std::cout<<"will add 3 in at the 3rd position"<<std::endl;
	//a.add(3,3);
	//std::cout<<a<<std::endl;
	//std::cout<<"will add 3 in at the first ordered position"<<std::endl;
	//a.add_sort(3,std::less<double>());
	//std::cout<<a<<std::endl;
	//std::cout<<"will add 0 in at the 0th position"<<std::endl;
	//a.add(0,0);
	//std::cout<<a<<std::endl;
//
	//std::cout<<"CHECK THAT I REALL MINIMIZE THE COPY/CREATION OF TYPE"<<std::endl;
}
