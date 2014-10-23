#include "Lapack.hpp"

/*private methods that depend on the type, used to call lapack*/
/*compute lu factorization : dgetrf zgetrf*/
/*{*/
template<>
void Lapack<double>::getrf(Vector<int>& ipiv){
	int info(1);
	dgetrf_(mat_->row(), mat_->col(), mat_->ptr(), mat_->row(), ipiv.ptr(), info);
	if(info){std::cerr<<"Lapack : getrf<double> : info="<<info<<std::endl;}
}

template<>
void Lapack<std::complex<double> >::getrf(Vector<int>& ipiv){
	int info(1);
	zgetrf_(mat_->row(), mat_->col(), mat_->ptr(), mat_->row(), ipiv.ptr(), info);
	if(info){std::cerr<<"Lapack : getrf<complex> : info="<<info<<std::endl;}
}
/*}*/

/*compute qr factorization : dgeqp3 dorgqr ?? ??*/
/*{*/
template<>
void Lapack<double>::geqp3(double* tau, int* jpvt){
	int info(1);
	int const lwork(4*mat_->col());
	double* work(new double[lwork]);
	
	dgeqp3_(mat_->row(), mat_->col(), mat_->ptr(), mat_->row(), jpvt, tau, work, lwork, info);
	if(info){std::cerr<<"Lapack : getrf<double> : info="<<info<<std::endl;}
	delete[] work;
}

template<>
void Lapack<double>::gqr(unsigned int k, double* tau){
	int info(1);
	//if(mat_->row()>=mat_->col()){ 
		int const lwork(mat_->col());
		double* work(new double[lwork]);
		unsigned int row(mat_->row());
		unsigned int col(mat_->col());
		if(col>row){ col=row; }
		dorgqr_(row, col, k, mat_->ptr(), row, tau, work, lwork, info); 
		delete[] work;
		if(info){std::cerr<<"Lapack : gqr<complex> : info="<<info<<std::endl;}
}

template<>
void Lapack<std::complex<double> >::geqp3(double* tau, int* jpvt){
	//int info(1);
	//zgetrf_(mat_->row(), mat_->col(), mat_->ptr(), mat_->row(), jptv, info);
	//if(info){std::cerr<<"Lapack : getrf<complex> : info="<<info<<std::endl; }
	std::cerr<<"Lapack : geqp3 : not implemented for Matrix<complex>"<<tau<<" "<<jpvt<<std::endl;
}

template<>
void Lapack<std::complex<double> >::gqr(unsigned int k, double* tau){
	//int info(1);
	//zgetrf_(mat_->row(), mat_->col(), mat_->ptr(), mat_->row(), jptv, info);
	//if(info){std::cerr<<"Lapack : getrf<complex> : info="<<info<<std::endl; }
	std::cerr<<"Lapack : gqr : not implemented for Matrix<complex>"<<k<<" "<<tau<<std::endl;
}
/*}*/
/*}*/

/*compute inverse of a matrix after using a lu decomposition : dgetri zgetri*/
/*{*/
template<>
void Lapack<double>::getri(Vector<int>& ipiv){
	unsigned int N(ipiv.size());
	unsigned int const lwork(3*N);
	double* work(new double[lwork]);
	int info(1);
	dgetri_(N, mat_->ptr(), N, ipiv.ptr(), work, lwork, info);
	if(info){std::cerr<<"Lapack : getri<double> : info="<<info<<std::endl;}
	delete[] work;
}

template<>
void Lapack<std::complex<double> >::getri(Vector<int>& ipiv){
	unsigned int N(ipiv.size());
	unsigned int const lwork(3*N);
	std::complex<double>* work(new std::complex<double>[lwork]);
	int info(1);
	zgetri_(N, mat_->ptr(), N, ipiv.ptr(), work, lwork, info);

	if(info){std::cerr<<"Lapack : getri<complex> : info="<<info<<std::endl;}
	delete[] work;
}
/*}*/

/*compute the condition number : dgecon zgecon*/
/*{*/
template<>
double Lapack<double>::gecon(double anorm){
	unsigned int N(mat_->col());
	int info(1);
	double rcond(0);
	double* work(new double[4*N]);
	int* iwork(new int[N]);
	dgecon_('1', N, mat_->ptr(), N, anorm, rcond, work, iwork, info);
	delete[] work;
	delete[] iwork;
	if(info){
		std::cerr<<"Lapack : gecon<double> : info="<<info<<std::endl;
		return 0;
	} else {
		return rcond;
	}
}

template<>
double Lapack<std::complex<double> >::gecon(double anorm){
	unsigned int N(mat_->col());
	int info(1);
	double rcond(0);
	std::complex<double> *work(new std::complex<double>[2*N]);
	double *rwork(new double[2*N]);
	zgecon_('1', N, mat_->ptr(), N, anorm, rcond, work, rwork, info);
	delete[] work;
	delete[] rwork;
	if(info){
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
	double *work(new double[mat_->row()]);
	double anorm(dlange_('1', mat_->row(), mat_->col(), mat_->ptr(), mat_->row(), work));
	delete[] work;
	return anorm;
}

template<>
double Lapack<std::complex<double> >::lange(){
	double *work(new double[mat_->row()]);
	double anorm(zlange_('1', mat_->row(), mat_->col(), mat_->ptr(), mat_->row(), work));
	delete[] work;
	return anorm;
}
/*}*/

/*compute the eigensystem : syev heev geev*/
/*{*/
template<>
void Lapack<double>::syev(Vector<double>& EVal, char job){
	unsigned int N(mat_->row());
	EVal.set(N);
	int lwork(-1);
	double wopt;
	int info(1);

	dsyev_(job, 'U', N, mat_->ptr(), N, EVal.ptr(), &wopt, lwork, info);
	lwork = int(wopt);
	double* work(new double[lwork]);
	dsyev_(job, 'U', N, mat_->ptr(), N, EVal.ptr(), work, lwork, info);

	if(info){ std::cerr<<"Lapack : syev : info="<<info<<std::endl; }
	delete[] work;
}

template<>
void Lapack<std::complex<double> >::syev(Vector<double>& EVal, char job){
	std::cerr<<"Lapack<double> : syev : complex symmetric not implemented"<<EVal<<" "<<job<<std::endl;
}

template<>
void Lapack<std::complex<double> >::heev(Vector<double>& EVal, char job){
	unsigned int N(mat_->row());
	EVal.set(N);
	int lwork(-1);
	std::complex<double> wopt;
	double* rwork(new double[3*N-2]);
	int info(1);

	zheev_(job, 'U', N, mat_->ptr(), N, EVal.ptr(), &wopt, lwork, rwork, info);
	lwork = int(wopt.real());
	std::complex<double>* work(new std::complex<double>[lwork]);
	zheev_(job, 'U', N, mat_->ptr(), N, EVal.ptr(), work, lwork, rwork, info);

	delete[] work;
	delete[] rwork;
	if(info){ std::cerr<<"Lapack : heev : info="<<info<<std::endl; }
}

template<>
void Lapack<double>::heev(Vector<double>& EVal, char job){
	std::cerr<<"Lapack<double>::heev() : heev : real Hermitian is not implemented"<<EVal<<" "<<job<<std::endl;
}

template<>
void Lapack<double>::geev(Vector<std::complex<double> >& EVal, Matrix<std::complex<double> >* REVec, Matrix<std::complex<double> >* LEVec){
	unsigned int N(mat_->row());
	unsigned int ldvl(1);
	unsigned int ldvr(1);
	double wopt;
	double* wr(new double[N]);
	double* wi(new double[N]);
	double* vr(NULL);
	double* vl(NULL);
	int lwork(-1);
	int info(1);
	char jobvr('N');
	char jobvl('N');
	if(REVec){
		REVec->set(N,N); 
		jobvr = 'V'; 
		ldvr = N;
		vr = new double[N*N];
	}
	if(LEVec){
		LEVec->set(N,N);
		jobvl = 'V';
		ldvl = N;
		vl = new double[N*N];
	}
	EVal.set(N);

	dgeev_(jobvl, jobvr, N, mat_->ptr(), N, wr, wi, vl, ldvl, vr, ldvr, &wopt, lwork, info); 
	lwork = int(wopt);
	double* work(new double[lwork]);
	dgeev_(jobvl, jobvr, N, mat_->ptr(), N, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info); 

	for(unsigned int i(0);i<N;i++){
		EVal(i) = std::complex<double>(wr[i],wi[i]);
	}

	if(REVec){
		REVec->set(N,N);
		for(unsigned int j(0);j<N;j++){
			if(j==N-1 || !(wi[j]*wi[j+1]<0 && are_equal(wr[j],wr[j+1]))){
				//std::cout<<j<<" real    "<<wr[j]<<" "<<wi[j]<<std::endl;
				for(unsigned int i(0);i<N;i++){ (*REVec)(i,j) = vr[i+j*N]; }
			} else {
				//std::cout<<j<<" complex "<<wi[j]<<" "<<wi[j+1]<<" "<<wr[j]<<" "<<wr[j+1]<<std::endl;
				for(unsigned int i(0);i<N;i++){
					(*REVec)(i,j)  =std::complex<double>(vr[i+j*N], vr[i+(1+j)*N]);
					(*REVec)(i,j+1)=std::complex<double>(vr[i+j*N],-vr[i+(1+j)*N]);
				}
				j++;
			}
		}
		delete[] vr;
	}
	if(LEVec){
		std::cerr<<"void Lapack<double>::geev(Vector<std::complex<double> >& EVal, Matrix<std::complex<double> >* REVec, Matrix<std::complex<double> >* LEVec) : might be a proble because LU*RU != 1 but LU*A*RU = diag"<<std::endl; 
		LEVec->set(N,N);
		for(unsigned int j(0);j<N;j++){
			if(j==N-1 || !(wi[j]*wi[j+1]<0 && are_equal(wr[j],wr[j+1]))){
				for(unsigned int i(0);i<N;i++){ (*LEVec)(j,i) = vl[i+j*N]; }
			} else {
				for(unsigned int i(0);i<N;i++){
					(*LEVec)(j,i)  =std::complex<double>(vl[i+j*N],-vl[i+(1+j)*N]);
					(*LEVec)(j+1,i)=std::complex<double>(vl[i+j*N], vl[i+(1+j)*N]);
				}
				j++;
			}
		}
		delete[] vl;
	}

	delete[] wr;
	delete[] wi;
	delete[] work;
	if(info){ std::cerr<<"Lapack : geev<double> : info="<<info<<std::endl; }
}

template<>
void Lapack<std::complex<double> >::geev(Vector<std::complex<double> >& EVal, Matrix<std::complex<double> >* REVec, Matrix<std::complex<double> >* LEVec){
	unsigned int N(mat_->row());
	unsigned int ldvl(1);
	unsigned int ldvr(1);
	std::complex<double> wopt;
	std::complex<double>* vr(NULL);
	std::complex<double>* vl(NULL);
	int lwork(-1);
	double* rwork(new double[2*N]);
	int info(1);
	char jobvr('N');
	char jobvl('N');
	if(REVec){ 
		REVec->set(N,N); 
		jobvr = 'V';
		ldvr = N;
		vr = REVec->ptr();
	}
	if(LEVec){ 
		LEVec->set(N,N); 
		jobvl = 'V';
		ldvl = N;
		vl = LEVec->ptr();
	}
	EVal.set(N);

	zgeev_(jobvl, jobvr, N, mat_->ptr(), N, EVal.ptr(), vl, ldvl, vr, ldvr, &wopt, lwork, rwork, info); 
	lwork = int(wopt.real());
	std::complex<double>* work(new std::complex<double>[N*N]);
	zgeev_(jobvl, jobvr, N, mat_->ptr(), N, EVal.ptr(), vl, ldvl, vr, ldvr, work, lwork, rwork, info); 

	delete[] work;

	if(info){ std::cerr<<"Lapack : geev<complex> : info="<<info<<std::endl; }
}
/*}*/
