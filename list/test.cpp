#include "List.hpp"
#include "Vector.hpp"
#include "Rand.hpp"

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

		~A(){ 
			std::cout<<"destroyed"<<std::endl;
			N_--; 
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
	IOFiles a("bla",true);
	//{
		//List<A> a;
		//for(unsigned int i(0);i<10;i++){ a.add_end(new A(i)); }
		//std::cout<<a<<std::endl;
		//std::cout<<A::N_<<std::endl;
		//std::cout<<"size of the list "<<a.size()<<std::endl;
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
		//List<A> b(new A(3));
		//std::cout<<b<<std::endl;
		//std::cout<<"will remove 1st entr with pop_start()"<<std::endl;
		//b.pop_start();
		//std::cout<<b<<std::endl;
		//std::cout<<"will swap the 2nd and 5th entry"<<std::endl;
		//a.swap(2,5);
		//std::cout<<a<<std::endl;
		//std::cout<<"will add 3 in at the 3rd position"<<std::endl;
		//a.add(new A(3),3);
		//std::cout<<a<<std::endl;
		//std::cout<<"will add 3 in at the first ordered position"<<std::endl;
		//auto func = [](A const* a, A const* b){ return a->a_<b->a_;};
		//a.add_sort(new A(3),func);
		//std::cout<<a<<std::endl;
		//std::cout<<"will add 10 in at the first ordered position"<<std::endl;
		//a.add_sort(new A(10),func);
		//std::cout<<a<<std::endl;
		//std::cout<<"will add 9 in at the first ordered position"<<std::endl;
		//a.add_sort(new A(9),func);
		//std::cout<<a<<std::endl;
		//std::cout<<"will add 0 in at the first ordered position"<<std::endl;
		//a.add_sort(new A(0),func);
		//std::cout<<a<<std::endl;
		//
		//std::cout<<"will add 0 in at the 0th position"<<std::endl;
		//a.add(new A(-1),0);
		//std::cout<<a<<std::endl;
//
		//std::cout<<"new list copied from the first one between enty 3 and 5"<<std::endl;
		//List<A> c(a.sublist(3,5));/*understand why it doesn't call a copy constructor*/
		//std::cout<<c<<std::endl;
		//std::cout<<"will remove the 2nd entry with pop(1)"<<std::endl;
		//c.pop(1);
		//std::cout<<c<<std::endl;
		//std::cout<<"will remove the 1st entry with pop(0)"<<std::endl;
		//c.pop(0);
		//std::cout<<c<<std::endl;
		//std::cout<<"will add at the 1st entry with add_start()"<<std::endl;
		//c.add_start(new A(8));
		//std::cout<<c<<std::endl;
		//std::cout<<"will add at the 1st entry with add_start()"<<std::endl;
		//c.add_start(new A(9));
		//std::cout<<c<<std::endl;
//
		//a.check();
	//}
	{
		List<int> a;
		Rand<int> rnd(0,100);
		for(unsigned int i(0);i<30;i++){ a.add_sort(new int(rnd.get()),[](int* a, int* b){ return *a<*b;} ); }
		std::cout<<a<<std::endl;
		std::cout<<a.last()<<std::endl;
	}
	{
		auto func = [](Vector<int>* a, Vector<int>* b) { 
			unsigned int i(0);
			while(i<a->size()){
				if((*a)(i) > (*b)(i)){ return false; }
				if((*a)(i) < (*b)(i)){ return true; }
				if((*a)(i) ==(*b)(i)){ i++; }
			}
			return false;
		};
		Rand<int> rnd(0,10);
		List<Vector<int> > a;
		Vector<int> tmp(2);
		for(unsigned int i(0);i<10;i++){
			tmp(0) = rnd.get();
			tmp(1) = rnd.get();
			a.add_sort(new Vector<int>(tmp), func);
			std::cout<<a<<std::endl;
		}
	}
	//{
		//auto cmp_for_fuse = [](Vector<int>* a, Vector<int>* b) { 
			//unsigned int i(0);
			//while(i<a->size()){
				//if((*a)(i) > (*b)(i)){ return 0; }
				//if((*a)(i) < (*b)(i)){ return 1; }
				//if((*a)(i)== (*b)(i)){ i++; }
			//}
			//std::cout<<"ca"<<std::endl;
			//return 2;
		//};
		//auto fuse = [](Vector<int>* a, Vector<int>* b) { 
			//std::cout<<"should fuse"<<std::endl;
			//std::cout<<*a<<std::endl;
			//std::cout<<*b<<std::endl;
		//};
		//Rand<int> rnd(0,10);
		//List<Vector<int> > a;
		//Vector<int> tmp(2);
		//for(unsigned int i(0);i<10;i++){
			//tmp(0) = rnd.get();
			//tmp(1) = rnd.get();
			//a.add_or_fuse_sort(new Vector<int>(tmp), cmp_for_fuse, fuse);
			//std::cout<<a<<std::endl;
		//}
	//}
	std::cout<<A::N_<<std::endl;
}
