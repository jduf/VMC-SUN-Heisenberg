#include "Miscellaneous.hpp"

int main(){
	std::cout<<my::tostring(my::norm_squared(std::complex<double>(1,2)))<<std::endl;
	Vector<std::complex<double> > a(30,1.0);
	Vector<std::complex<double> > b(30,1.0);
	std::cout<<BLAS::dot(a.ptr(),1,0,b.ptr(),1,0)<<std::endl;
}
