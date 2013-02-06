#include "Lapack.hpp"
/*private methods that depend on the type, used to call lapack*/
/*{*/
/*general matrices*/
/*{*/
/*compute lu factorization : dgetrf zgetrf*/
/*{*/
template<>
void Lapack<double>::getrf(int ipiv[]) const {
	int info;
	dgetrf_(&N, &N, m, &N, ipiv, &info);
}
template<>
void Lapack<std::complex<double> >::getrf(int ipiv[]) const {
	int info;
	zgetrf_(&N, &N, m, &N, ipiv, &info);
}
/*}*/
/*compute inverse of a matrix after using a lu decomposition : dgetri zgetri*/
/*{*/
template<>
void Lapack<double>::getri(int ipiv[]) const {
	unsigned int const lwork(3*N-1);
	double work[lwork];
	int info;
	dgetri_(&N, m, &N, ipiv, work,  &lwork, &info);
}
template<>
void Lapack<std::complex<double> >::getri(int ipiv[]) const {
	unsigned int const lwork(3*N-1);
	std::complex<double> work[lwork];
	int info;
	zgetri_(&N, m, &N, ipiv, work,  &lwork, &info);
}
/*}*/
/*}*/
/*symetric matrices*/
/*{*/
template<>
void Lapack<double>::syev(Vecteur<double> & EVal) const {
	char jobz('V'),uplo('U');
	unsigned int const lwork(3*N-1);
	double work[lwork];
	int info;
	dsyev_(&jobz, &uplo, &N, m, &N, EVal.ptr(), work ,&lwork, &info);
}
/*}*/
/*}*/
