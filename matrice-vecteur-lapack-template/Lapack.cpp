#include "Lapack.hpp"
/*private methods that depend on the type, used to call lapack*/
/*general matrices*/
/*{*/
/*compute lu factorization : dgetrf zgetrf*/
/*{*/
template<>
void Lapack<double>::getrf(int *ipiv) const {
	int info(1);
	dgetrf_(N, N, m, N, ipiv, info);
	if(info !=0){std::cerr<<"Lapack : getrf<double> : info="<<info<<std::endl; }
}

template<>
void Lapack<std::complex<double> >::getrf(int *ipiv) const {
	int info(1);
	zgetrf_(N, N, m, N, ipiv, info);
	if(info !=0){std::cerr<<"Lapack : getrf<complex> : info="<<info<<std::endl; }
}
/*}*/

/*compute inverse of a matrix after using a lu decomposition : dgetri zgetri*/
/*{*/
template<>
void Lapack<double>::getri(int *ipiv) const {
	unsigned int const lwork(3*N-1);
	double* work(new double[lwork]);
	int info(1);
	dgetri_(N, m, N, ipiv, work, lwork, info);
	if(info !=0){std::cerr<<"Lapack : getri<double> : info="<<info<<std::endl;}
	delete[] work;
}

template<>
void Lapack<std::complex<double> >::getri(int *ipiv) const {
	unsigned int const lwork(3*N-1);
	std::complex<double>* work(new std::complex<double>[lwork]);
	int info(1);
	zgetri_(N, m, N, ipiv, work, lwork, info);
	if(info !=0){std::cerr<<"Lapack : getri<complex> : info="<<info<<std::endl;}
	delete[] work;
}
/*}*/

/*compute the condition number : dgecon zgecon*/
/*{*/
template<>
double Lapack<double>::gecon(double anorm) const {
	int info(1);
	double rcond(0);
	double* work(new double[4*N]);
	int* iwork(new int[N]);
	dgecon_('0', N, m, N, anorm, rcond, work, iwork, info);
	if(info !=0){
		std::cerr<<"Lapack : getrf<double> : info="<<info<<std::endl;
		return 0;
	} else {
		return rcond;
	}
}

template<>
double Lapack<std::complex<double> >::gecon(double anorm) const {
	std::cerr<<"Lapack : gecon<std::complex<double> > : not implemented"<<anorm<<std::endl;
	return 0;
}
/*}*/

/*compute the norm : dlange zlange*/
/*{*/
template<>
double Lapack<double>::lange() const {
	double *work(new double[4*N]);
	return dlange_('0', N, N, m, N, work);
}

template<>
double Lapack<std::complex<double> >::lange() const {
	std::cerr<<"Lapack : lange<std::complex<double> > : not implemented"<<std::endl;
	return 0;
}
/*}*/
/*}*/

/*public methods that depend on the type, used to call lapack*/
template<> 
void Lapack<double>::eigensystem(Vecteur<double>& EVal, bool EVec) const {
	char jobz('N');
	if(EVec){jobz='V';}
	switch(matrix_type){
		case 'S':
			{
				int const lwork(3*N-1);
				double* work(new double[lwork]);
				int info(1);
				dsyev_(jobz, 'U', N, m, N, EVal.ptr(), work, lwork, info);
				if(info !=0) { std::cerr<<"Lapack : eigensystem<double> : info="<<info<<std::endl; }
				delete[] work;
				break;
			}
		default:
			{
				std::cerr<<"eigensystem : Matrix type "<<matrix_type<<" not implemented for real matrix"<<std::endl;
				std::cerr<<"eigensystem : the only matrix type implemented is S"<<std::endl;
				break;
			}
	}
}

template<> 
void Lapack<std::complex<double> >::eigensystem(Vecteur<double>& EVal, bool EVec) const {
	char jobz('N');
	if(EVec){jobz='V';}
	switch(matrix_type){
		case 'H':
			{
				int const lwork(2*N-1);
				std::complex<double>* work(new std::complex<double>[lwork]);
				double* rwork(new double[3*N-2]);
				int info(1);
				zheev_(jobz, 'U', N, m, N, EVal.ptr(), work, lwork, rwork, info);
				delete[] work;
				delete[] rwork;
				if(info !=0) { std::cerr<<"Lapack : eigensystem<complex> : info="<<info<<std::endl; }
				break;
			}
		default:
			{
				std::cerr<<"eigensystem : Matrix type "<<matrix_type<<" not implemented for complex matrix"<<std::endl;
				std::cerr<<"eigensystem : the only matrix type implemented is H"<<std::endl;
				break;
			}
	}
}
