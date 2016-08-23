#include "Vector.hpp"
#include "Rand.hpp"

int main(){
	//Rand<double> r(0,100);
	//Rand<unsigned int> ur(10,100);
	//Vector<std::complex<double> > a(ur.get());
	//Vector<std::complex<double> > b(a.size());
	//for(unsigned int i(0);i<a.size();i++){
		//a(i) = std::complex<double>(r.get(),r.get());
		//b(i) = std::complex<double>(r.get(),r.get());
	//}

	Matrix<std::complex<double> > a(3,31);
	Vector<std::complex<double> > b(31);
	for(unsigned int i(0);i<b.size();i++){
		a(1,i) = std::complex<double>(0,1.*i);
		b(i) = std::complex<double>(1.*i,0);
	}
	//a(1,0) = std::complex<double>(10,11);
	//a(1,1) = std::complex<double>(12,13);

	std::cout<<BLAS::dot(b.size(),a.ptr(),true,a.row(),1,b.ptr(),true,1,0)<<std::endl;
}
