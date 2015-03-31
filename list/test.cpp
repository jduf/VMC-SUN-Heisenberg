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

		A const& operator*=(double x){ 
			a_ *= x;
			return (*this);
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
	{
		std::cout<<"#### general tests ####"<<std::endl;
		List<A> a;
		for(unsigned int i(0);i<10;i++){ a.add_end(std::make_shared<A>(i)); }
		std::cout<<a<<std::endl;
		std::cout<<A::N_<<std::endl;
		std::cout<<"size of the list "<<a.size()<<std::endl;
		std::cout<<"will remove 4th entry with pop(4)"<<std::endl;
		a.pop(4);
		std::cout<<a<<std::endl;
		std::cout<<"will remove last entry with pop(idx)"<<std::endl;
		a.pop(a.size()-1);
		std::cout<<a<<std::endl;
		std::cout<<"will remove last entry with pop_end"<<std::endl;
		a.pop_end();
		std::cout<<a<<std::endl;
		std::cout<<"will remove 1st entry with pop(0)"<<std::endl;
		a.pop(0);
		std::cout<<a<<std::endl;
		std::cout<<"will swap entry 2 end 5"<<std::endl;
		a.swap(2,5);
		std::cout<<a<<std::endl;
		std::cout<<"will add 3 in at the 3rd position"<<std::endl;
		a.add(std::make_shared<A>(3),3);
		std::cout<<a<<std::endl;
		std::cout<<"will add 3 in at the first ordered position"<<std::endl;
		auto func = [](A const& a, A const& b){ return a.a_<b.a_;};
		a.add_sort(std::make_shared<A>(3),func);
		std::cout<<a<<std::endl;
		std::cout<<"will add 10 in at the first ordered position"<<std::endl;
		a.add_sort(std::make_shared<A>(10),func);
		std::cout<<a<<std::endl;
		std::cout<<"will add 0 in at the first ordered position"<<std::endl;
		a.add_sort(std::make_shared<A>(0),func);
		std::cout<<a<<std::endl;
		std::cout<<"will add -1 in at the first ordered position"<<std::endl;
		a.add_sort(std::make_shared<A>(-1),func);
		std::cout<<a<<std::endl;
		std::cout<<"will add 9 in at the first ordered position"<<std::endl;
		a.add_sort(std::make_shared<A>(9),func);
		std::cout<<a<<std::endl;
		std::cout<<"will add 0 in at the first ordered position"<<std::endl;
		a.add_sort(std::make_shared<A>(0),func);
		std::cout<<a<<std::endl;
		std::cout<<"will print using get_next must be identical to the previous line"<<std::endl;
		do{
			std::cout<<a.get()<<" ";
		} while ( a.move_forward());
		std::cout<<std::endl;
		std::cout<<"will copy the whole first list without actually copying the value."<<std::endl;
		List<A> c;
		do{
			c.add_end(a.get_ptr());
		} while ( a.move_forward());
		do{
			std::cout<<c.get()<<" ";
		} while ( c.move_forward());
		std::cout<<std::endl;

		std::cout<<"multiply the first entry bigger than 2 by pi"<<std::endl;
		while( c.get().a_ < 2.0 && c.move_forward());
		c.get() *= M_PI;
		std::cout<<c<<std::endl;
		std::cout<<"the other list should also be affected"<<std::endl;
		std::cout<<a<<std::endl;


		std::cout<<"list copied from the first one between enty [3,7)"<<std::endl;
		List<A> b(a.sublist(3,7));/*understand why it doesn't call a copy constructor*/
		std::cout<<b<<std::endl;
		std::cout<<"will remove the 2nd entry with pop(1)"<<std::endl;
		b.pop(1);
		std::cout<<b<<std::endl;
		std::cout<<"will remove the 1st entry with pop(0)"<<std::endl;
		b.pop(0);
		std::cout<<b<<std::endl;
		std::cout<<"will add 8 at the 1st entry with add_start()"<<std::endl;
		b.add_start(std::make_shared<A>(8));
		std::cout<<b<<std::endl;
		std::cout<<"will remove last entry with pop(idx)"<<std::endl;
		b.pop(b.size()-1);
		std::cout<<b<<std::endl;
		std::cout<<"will remove last entry with pop(idx)"<<std::endl;
		b.pop(b.size()-1);
		std::cout<<b<<std::endl;
		std::cout<<"will remove last entry with pop(idx)"<<std::endl;
		b.pop(b.size()-1);
		std::cout<<b<<std::endl;
		std::cout<<"will remove last entry with pop(idx)"<<std::endl;
		b.pop(b.size()-1);
		std::cout<<b<<std::endl;
		std::cout<<"will add at the last entry with add_end()"<<std::endl;
		b.add_end(std::make_shared<A>(9));
		std::cout<<b<<std::endl;


	}
	{
		std::cout<<"#### test add_sort with int ####"<<std::endl;
		List<int> a;
		Rand<int> rnd(0,100);
		auto cmp = [](int const& a, int const& b){ return a<b;};
		for(unsigned int i(0);i<30;i++){ a.add_sort(std::make_shared<int>(rnd.get()),cmp ); }
		std::cout<<a<<std::endl;
	}
	{
		std::cout<<"#### test add_sort with vector ####"<<std::endl;
		auto func = [](Vector<int> const& a, Vector<int> const& b) { 
			unsigned int i(0);
			while(i<a.size()){
				if(a(i) > b(i)){ return false; }
				if(a(i) < b(i)){ return true; }
				if(a(i) ==b(i)){ i++; }
			}
			return false;
		};
		Rand<int> rnd(0,10);
		List<Vector<int> > a;
		Vector<int> tmp(2);
		for(unsigned int i(0);i<10;i++){
			tmp(0) = rnd.get();
			tmp(1) = rnd.get();
			a.add_sort(std::make_shared<Vector<int> >(tmp), func);
			std::cout<<a<<std::endl;
		}
	}
	{
		std::cout<<"#### test add_or_fuse_sort ####"<<std::endl;
		auto cmp_for_fuse = [](Vector<int> const& a, Vector<int> const& b) { 
			unsigned int i(0);
			while(i<a.size()){
				if(a(i) > b(i)){ return 0; }
				if(a(i) < b(i)){ return 1; }
				if(a(i)== b(i)){ i++; }
			}
			return 2;
		};
		auto fuse = [](Vector<int>& a, Vector<int> const& b) { 
			std::cout<<"should fuse"<<std::endl;
			(void)(a);
			(void)(b);
		};
		Rand<int> rnd(0,10);
		List<Vector<int> > a;
		Vector<int> tmp(2);
		std::shared_ptr<Vector<int> > tmp_shared;
		for(unsigned int i(0);i<10;i++){
			tmp(0) = rnd.get();
			tmp(1) = rnd.get();
			tmp_shared = std::make_shared<Vector<int> >(tmp);
			if( a.find(tmp_shared, cmp_for_fuse) ){ a.fuse_with_move(tmp_shared,fuse); }
			else { a.add_after_move(tmp_shared); }

			std::cout<<a<<std::endl;
		}
	}
	std::cout<<"#(constructor calls)-#(destructor calls)="<<A::N_<<std::endl;
}
