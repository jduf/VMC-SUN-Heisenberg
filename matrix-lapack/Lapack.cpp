#include "Lapack.hpp"

/*private methods that depend on the type, used to call lapack*/
/*general Matrixs*/
/*{*/
/*compute lu factorization : dgetrf zgetrf*/
/*{*/
template<>
void Lapack<double>::getrf(int *ipiv){
	int info(1);
	dgetrf_(mat->row(), mat->col(), mat->ptr(), mat->row(), ipiv, info);
	if(info !=0){std::cerr<<"Lapack : getrf<double> : info="<<info<<std::endl; }
}

template<>
void Lapack<std::complex<double> >::getrf(int *ipiv){
	int info(1);
	zgetrf_(mat->row(), mat->col(), mat->ptr(), mat->row(), ipiv, info);
	if(info !=0){std::cerr<<"Lapack : getrf<complex> : info="<<info<<std::endl; }
}
/*}*/

/*compute qr factorization : dgeqp3 zgeqp3*/
/*{*/
template<>
void Lapack<double>::geqp3(Matrix<double>& Q, Matrix<double>& R){
	int info(1);
	int const lwork(4*mat->col());
	double* work(new double[lwork]);
	unsigned int k(std::min(mat->row(),mat->col()));
	double *tau(new double[k]);
	int *jptv(new int[mat->col()]);
	//for(unsigned int j(0);j<mat->col();j++){
		//jptv[j] = k;
	//}
	
	dgeqp3_(mat->row(), mat->col(), mat->ptr(), mat->row(), jptv, tau, work, lwork, info);

	for(unsigned int i(0); i<mat->row(); i++){ Q(i,i) = 1.0;}
	Matrix<double> v(mat->row(),1,0);
	Matrix<double> H(mat->row(),mat->row());
	for(unsigned int i(0);i<k;i++){
		v(i)=1;
		for(unsigned int j(i+1);j<mat->row();j++){
			v(j) = (*mat)(j,i);
		}
		H = v * v.transpose();
		H *= -tau[i];
		for(unsigned int j(0);j<mat->row();j++){
			H(j,j) = 1+H(j,j);
		}
		Q *= H;
		v(i) = 0;
	}
	for(unsigned int i(0);i<R.col();i++){
		for(unsigned int j(i);j<R.col();j++){
			R(i,j) = (*mat)(i,j);
		}
	}
	//for(unsigned int i(0);i<k;i++){
		//std::cout<<jptv[i]<<" ";
	//}
	//std::cout<<std::endl;	
	if(info !=0){std::cerr<<"Lapack : getrf<double> : info="<<info<<std::endl;}
}

template<>
void Lapack<std::complex<double> >::geqp3(Matrix<std::complex<double> >& Q, Matrix<std::complex<double> >& R){
	//int info(1);
	//zgetrf_(mat->row(), mat->col(), mat->ptr(), mat->row(), jptv, info);
	//if(info !=0){std::cerr<<"Lapack : getrf<complex> : info="<<info<<std::endl; }
	std::cerr<<"Lapack : geqp3 : not implemented for Matrix<complex>"<<std::endl;
}
/*}*/

/*compute inverse of a matrix after using a lu decomposition : dgetri zgetri*/
/*{*/
template<>
void Lapack<double>::getri(int *ipiv){
	unsigned int N(mat->row());
	unsigned int const lwork(3*N);
	double* work(new double[lwork]);
	int info(1);
	dgetri_(N, mat->ptr(), N, ipiv, work, lwork, info);
	if(info !=0){std::cerr<<"Lapack : getri<double> : info="<<info<<std::endl;}
	delete[] work;
}

template<>
void Lapack<std::complex<double> >::getri(int *ipiv) {
		unsigned int N(mat->row());
	unsigned int const lwork(3*N);
	std::complex<double>* work(new std::complex<double>[lwork]);
	int info(1);
	zgetri_(N, mat->ptr(), N, ipiv, work, lwork, info);
	if(info !=0){std::cerr<<"Lapack : getri<complex> : info="<<info<<std::endl;}
	delete[] work;
}
/*}*/

/*compute the condition number : dgecon zgecon*/
/*{*/
template<>
double Lapack<double>::gecon(double anorm) {
		unsigned int N(mat->row());
	int info(1);
	double rcond(0);
	double* work(new double[4*N]);
	int* iwork(new int[N]);
	dgecon_('0', N, mat->ptr(), N, anorm, rcond, work, iwork, info);
	if(info !=0){
		std::cerr<<"Lapack : gecon<double> : info="<<info<<std::endl;
		return 0;
	} else {
		return rcond;
	}
}

template<>
double Lapack<std::complex<double> >::gecon(double anorm) {
	std::cerr<<"Lapack : gecon<std::complex<double> > : not implemented"<<anorm<<std::endl;
	return 0;
}
/*}*/

/*compute the norm : dlange zlange*/
/*{*/
template<>
double Lapack<double>::lange() {
	double *work(new double[4*mat->row()]);
	return dlange_('0', mat->row() , mat->col(), mat->ptr(), mat->row(), work);
}

template<>
double Lapack<std::complex<double> >::lange() {
	std::cerr<<"Lapack : lange<std::complex<double> > : not implemented"<<std::endl;
	return 0;
}
/*}*/
/*}*/

/*public methods that depend on the type, used to call lapack*/
template<> 
void Lapack<double>::eigensystem(Matrix<double>& EVal, bool EVec) {
	char jobz('N');
	if(EVec){jobz='V';}
	switch(matrix_type){
		case 'S':
			{
		unsigned int N(mat->row());
				int const lwork(3*N);
				double* work(new double[lwork]);
				int info(1);
				dsyev_(jobz, 'U', N, mat->ptr(), N, EVal.ptr(), work, lwork, info);
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
void Lapack<std::complex<double> >::eigensystem(Matrix<double>& EVal, bool EVec) {
	char jobz('N');
	if(EVec){jobz='V';}
	switch(matrix_type){
		case 'H':
			{
		unsigned int N(mat->row());
				int const lwork(2*N);
				std::complex<double>* work(new std::complex<double>[lwork]);
				double* rwork(new double[3*N-2]);
				int info(1);
				zheev_(jobz, 'U', N, mat->ptr(), N, EVal.ptr(), work, lwork, rwork, info);
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
