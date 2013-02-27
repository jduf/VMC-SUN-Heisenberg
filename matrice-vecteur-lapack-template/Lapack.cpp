#include "Lapack.hpp"
/*private methods that depend on the type, used to call lapack*/
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

/*public methods that depend on the type, used to call lapack*/
template<> 
void Lapack<double>::eigensystem(Vecteur<double>& EVal, bool EVec) const {
	char jobz('V');
	if(!EVec){jobz='N';}
	switch(matrix_type){
		case 'S':
			{
				char uplo('U');
				unsigned int const lwork(3*N-1);
				double work[lwork];
				int info;
				dsyev_(&jobz, &uplo, &N, m, &N, EVal.ptr(), work ,&lwork, &info);
			}
			break;
		default:
			{
				std::cerr<<"eigensystem : Matrix type "<<matrix_type<<" not implemented for real matrix"<<std::endl;
				std::cerr<<"eigensystem : the only matrix type implemented is S"<<std::endl;
			}
			break;
	}
}

template<> 
void Lapack<std::complex<double> >::eigensystem(Vecteur<double>& EVal, bool EVec) const {
	char jobz('V');
	if(!EVec){jobz='N';}
	switch(matrix_type){
		case 'H':
			{
				char uplo('U');
				unsigned int const lwork(2*N-1);
				std::complex<double> work[lwork];
				double rwork[3*N-2];
				int info;
				zheev_(&jobz, &uplo, &N, m, &N, EVal.ptr(), work, &lwork, rwork, &info);
			}
			break;
		default:
			{
				std::cerr<<"eigensystem : Matrix type "<<matrix_type<<" not implemented for complex matrix"<<std::endl;
				std::cerr<<"eigensystem : the only matrix type implemented is H"<<std::endl;
			}
			break;
	}
}
