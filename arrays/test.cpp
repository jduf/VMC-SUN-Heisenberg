#include "Array1D.hpp"
#include "Array2D.hpp"
#include "Write.hpp"
#include "Read.hpp"

#include <complex>

template <typename T> void out(Array2D<T> const& t);

template <typename T> void in(Array2D<T>& t);

int main(){
	Array2D<int> t(5,2);
	Array2D<int> t_read(3,2);
	t(0,0) = 1;
	t(0,0) = 1;
	t(0,1) = 2;
	t(1,0) =-3;
	t(1,1) = 4;
	t(2,0) = 5;
	t(2,1) =-6;
	t(3,0) = 7;
	t(3,1) = 8;
	t(4,0) =-9;
	t(4,1) = 0;

	out(t);

	in(t_read);
	Array2D<std::complex<double> > c(2,6,std::complex<double>(3.1,4.3));	
	Array2D<std::complex<double> > c_read(2,6);	
	
	out(c);
	in(c_read);	
}


template <typename T>
void out(Array2D<T> const& t){
	Write w("blalba",false);
	w<<t;
}

template <typename T>
void in(Array2D<T>& t){
	Read r("blalba",false);
	r>>t;
	std::cout<<t<<std::endl;
}
