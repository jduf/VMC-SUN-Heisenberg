#include "Matrice.hpp"
template<>
Matrice<std::complex<double> > Matrice<std::complex<double> >::trans_conj() const{
	Matrice<std::complex<double> > tmp(N);
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			tmp.m[i+j*N] = conj(m[j+i*N]);
		}
	}
	return tmp;
}

template<>
Matrice<double> Matrice<double>::trans_conj() const{
	std::cerr<<"conjugate : matrix conjugate of a real matrix doesn't exist"<<std::endl;
	return (*this);
}
