#include "Lapack.hpp"

/*private methods that depend on the type, used to call lapack*/
/*general Matrixs*/
/*{*/
/*compute lu factorization : dgetrf zgetrf*/
/*{*/
template<>
void Lapack<double>::getrf(Matrix<int>& ipiv){
	int info(1);
	dgetrf_(mat->row(), mat->col(), mat->ptr(), mat->row(), ipiv.ptr(), info);
	if(info !=0){std::cerr<<"Lapack : getrf<double> : info="<<info<<std::endl; }
}

template<>
void Lapack<std::complex<double> >::getrf(Matrix<int>& ipiv){
	int info(1);
	zgetrf_(mat->row(), mat->col(), mat->ptr(), mat->row(), ipiv.ptr(), info);
	if(info !=0){std::cerr<<"Lapack : getrf<complex> : info="<<info<<std::endl; }
}
/*}*/

/*compute qr factorization : dgeqp3 dorgqr ?? ??*/
/*{*/
template<>
void Lapack<double>::geqp3(double* tau, int* jpvt){
	int info(1);
	int const lwork(4*mat->col());
	double* work(new double[lwork]);
	
	dgeqp3_(mat->row(), mat->col(), mat->ptr(), mat->row(), jpvt,tau, work, lwork, info);
	if(info !=0){std::cerr<<"Lapack : getrf<double> : info="<<info<<std::endl;}
	delete[] work;
}

template<>
void Lapack<double>::gqr(unsigned int k, double* tau){
	int info(1);
	//if(mat->row()>=mat->col()){ 
		int const lwork(mat->col());
		double* work(new double[lwork]);
		unsigned int row(mat->row());
		unsigned int col(mat->col());
		if(col>row){ col=row; }
		dorgqr_(row, col, k, mat->ptr(), row, tau, work, lwork, info); 
		delete[] work;
		if(info !=0){std::cerr<<"Lapack : gqr<complex> : info="<<info<<std::endl;}
}

template<>
void Lapack<std::complex<double> >::geqp3(double* tau, int* jpvt){
	//int info(1);
	//zgetrf_(mat->row(), mat->col(), mat->ptr(), mat->row(), jptv, info);
	//if(info !=0){std::cerr<<"Lapack : getrf<complex> : info="<<info<<std::endl; }
	std::cerr<<"Lapack : geqp3 : not implemented for Matrix<complex>"<<std::endl;
}

template<>
void Lapack<std::complex<double> >::gqr(unsigned int k, double* tau){
	//int info(1);
	//zgetrf_(mat->row(), mat->col(), mat->ptr(), mat->row(), jptv, info);
	//if(info !=0){std::cerr<<"Lapack : getrf<complex> : info="<<info<<std::endl; }
	std::cerr<<"Lapack : gqr : not implemented for Matrix<complex>"<<std::endl;
}
/*}*/
/*}*/

/*compute inverse of a matrix after using a lu decomposition : dgetri zgetri*/
/*{*/
template<>
void Lapack<double>::getri(Matrix<int>& ipiv){
	unsigned int N(mat->row());
	unsigned int const lwork(3*N);
	double* work(new double[lwork]);
	int info(1);
	dgetri_(N, mat->ptr(), N, ipiv.ptr(), work, lwork, info);
	if(info !=0){std::cerr<<"Lapack : getri<double> : info="<<info<<std::endl;}
	delete[] work;
}

template<>
void Lapack<std::complex<double> >::getri(Matrix<int>& ipiv) {
	unsigned int N(mat->row());
	unsigned int const lwork(3*N);
	std::complex<double>* work(new std::complex<double>[lwork]);
	int info(1);
	zgetri_(N, mat->ptr(), N, ipiv.ptr(), work, lwork, info);

	if(info !=0){std::cerr<<"Lapack : getri<complex> : info="<<info<<std::endl;}
	delete[] work;
}
/*}*/

/*compute the condition number : dgecon zgecon*/
/*{*/
template<>
double Lapack<double>::gecon(double anorm) {
	unsigned int N(mat->col());
	int info(1);
	double rcond(0);
	double* work(new double[4*N]);
	int* iwork(new int[N]);
	dgecon_('1', N, mat->ptr(), N, anorm, rcond, work, iwork, info);
	delete[] work;
	delete[] iwork;
	if(info !=0){
		std::cerr<<"Lapack : gecon<double> : info="<<info<<std::endl;
		return 0;
	} else {
		return rcond;
	}
}

template<>
double Lapack<std::complex<double> >::gecon(double anorm) {
	unsigned int N(mat->col());
	int info(1);
	double rcond(0);
	std::complex<double> *work(new std::complex<double>[2*N]);
	double *rwork(new double[2*N]);
	zgecon_('1', N, mat->ptr(), N, anorm, rcond, work, rwork, info);
	delete[] work;
	delete[] rwork;
	if(info !=0){
		std::cerr<<"Lapack : gecon<double> : info="<<info<<std::endl;
		return 0;
	} else {
		return rcond;
	}
}
/*}*/

/*compute the norm : dlange zlange*/
/*{*/
template<>
double Lapack<double>::lange(){
	double *work(new double[mat->row()]);
	double anorm(dlange_('1', mat->row(), mat->col(), mat->ptr(), mat->row(), work));
	delete[] work;
	return anorm;
}

template<>
double Lapack<std::complex<double> >::lange() {
	double *work(new double[mat->row()]);
	double anorm(zlange_('1', mat->row(), mat->col(), mat->ptr(), mat->row(), work));
	delete[] work;
	return anorm;
}
/*}*/

/*compute the eigensystem : syev heev geev*/
/*{*/
template<>
void Lapack<double>::syev(Matrix<double>& EVal, char job){
	unsigned int N(mat->row());
	EVal.set(N,1);
	int lwork(-1);
	double wopt;
	int info(1);
	dsyev_(job, 'U', N, mat->ptr(), N, EVal.ptr(), &wopt, lwork, info);
	lwork = int(wopt);
	double* work(new double[lwork]);
	dsyev_(job, 'U', N, mat->ptr(), N, EVal.ptr(), work, lwork, info);
	if(info!=0) { std::cerr<<"Lapack : eigensystem<double> : info="<<info<<std::endl; }
	delete[] work;
}

template<>
void Lapack<std::complex<double> >::syev(Matrix<double>& EVal, char job){
	EVal.set();
	std::cerr<<"Lapack<double> : syev : complex symmetric not implemented"<<std::endl;
}

template<>
void Lapack<std::complex<double> >::heev(Matrix<double>& EVal, char job){
	unsigned int N(mat->row());
	EVal.set(N,1);
	int lwork(-1);
	std::complex<double> wopt;
	double* rwork(new double[3*N-2]);
	int info(1);
	zheev_(job, 'U', N, mat->ptr(), N, EVal.ptr(), &wopt, lwork, rwork, info);
	lwork = int(wopt.real());
	std::complex<double>* work(new std::complex<double>[lwork]);
	zheev_(job, 'U', N, mat->ptr(), N, EVal.ptr(), work, lwork, rwork, info);
	delete[] work;
	delete[] rwork;
	if(info!=0) { std::cerr<<"Lapack : eigensystem<complex> : info="<<info<<std::endl; }
}

template<>
void Lapack<double>::heev(Matrix<double>& EVal, char job){
	std::cerr<<"Lapack<double> : heev : real Hermitian is not implemented"<<std::endl;
}

template<>
void Lapack<double>::geev(Matrix<std::complex<double> >& EVal, char job, Matrix<std::complex<double> >& EVec){
	unsigned int N(mat->row());
	EVal.set(N,1);
	double* wr(new double[N]);
	double* wi(new double[N]);
	double* vl(new double[N]); /*not used if jobvl='N'*/
	Matrix<double> vr(N,N); /*not used if jobvr='N'*/
	int lwork(-1);
	double wopt;
	int info(1);

	dgeev_('N', job, N, mat->ptr(), N, wr, wi, vl, 1, vr.ptr(), N, &wopt, lwork, info); 
	lwork = int(wopt);
	std::cout<<lwork<<std::endl;
	double* work(new double[lwork]);
	dgeev_('N', job, N, mat->ptr(), N, wr, wi, vl, 1, vr.ptr(), N, work, lwork, info); 

	for(unsigned int i(0);i<N;i++){
		EVal(i) = std::complex<double>(wr[i],wi[i]);
	}

	EVec.set(N,N);
	for(unsigned int i(0);i<N;i++){
		if( std::abs(EVal(i).imag())<1e-14){
			for(unsigned int j(0);j<N;j++){
				EVec(j,i) = vr(j,i);
			}
		} else {
			for(unsigned int j(0);j<N;j++){
				EVec(j,i) = std::complex<double>(vr(j,i),vr(j,(1+i)));
				EVec(j,i+1) = std::complex<double>(vr(j,i),-vr(j,(1+i)));
			}
			i++;
		}
	}

	delete[] wr;
	delete[] wi;
	delete[] vl;
	delete[] work;
}

template<>
void Lapack<std::complex<double> >::geev(Matrix<std::complex<double> >& EVal, char job, Matrix<std::complex<double> >& EVec){
	std::cerr<<"Lapack<double> : geev : complex general is not implemented"<<std::endl;
}
/*}*/
/*}*/
