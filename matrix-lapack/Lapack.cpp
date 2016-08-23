#include "Lapack.hpp"

/*private methods that depend on the type, used to call lapack*/
/*compute lu factorization : dgetrf zgetrf*/
/*{*/
template<>
void Lapack<double>::getrf(Vector<int>& ipiv){
	int info(1);
	dgetrf_(mat_->row(), mat_->col(), mat_->ptr(), mat_->row(), ipiv.ptr(), info);
	if(info<0){
		std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl;
		ipiv.set();
	}
}

template<>
void Lapack<std::complex<double> >::getrf(Vector<int>& ipiv){
	int info(1);
	zgetrf_(mat_->row(), mat_->col(), mat_->ptr(), mat_->row(), ipiv.ptr(), info);
	if(info<0){
		std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl;
		ipiv.set();
	}
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
	if(info){ std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl; }
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
	if(info){ std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl; }
}

template<>
void Lapack<std::complex<double> >::geqp3(double* tau, int* jpvt){
	(void)(tau);
	(void)(jpvt);
	std::cerr<<__PRETTY_FUNCTION__<<" : not implemented for complex"<<std::endl;
}

template<>
void Lapack<std::complex<double> >::gqr(unsigned int k, double* tau){
	(void)(k);
	(void)(tau);
	std::cerr<<__PRETTY_FUNCTION__<<" : not implemented for complex"<<std::endl;
}
/*}*/
/*}*/

/*compute Bunch-Kaufman factorization : dsytrf*/
/*{*/
template<>
void Lapack<double>::sytrf(Vector<int>& ipiv){
	int info(1);
	int lwork(-1);
	double wopt;

	dsytrf_('U', mat_->col(), mat_->ptr(), mat_->row(), ipiv.ptr(), &wopt, lwork, info);
	if(info){
		std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl;
		ipiv.set();
	} else {
		lwork = int(wopt);
		double* work(new double[lwork]);
		dsytrf_('U', mat_->col(), mat_->ptr(), mat_->row(), ipiv.ptr(), work, lwork, info);

		if(info<0){
			std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl;
			ipiv.set();
		}
		delete[] work;
	}
}

template<>
void Lapack<std::complex<double> >::sytrf(Vector<int>& ipiv){
	(void)(ipiv);
	std::cerr<<__PRETTY_FUNCTION__<<" : not implemented for complex"<<std::endl;
}
/*}*/

/*compute inverse of a matrix after using a lu factorization : dgetri zgetri*/
/*{*/
template<>
void Lapack<double>::getri(Vector<int>& ipiv){
	unsigned int N(ipiv.size());
	unsigned int const lwork(3*N);
	double* work(new double[lwork]);
	int info(1);
	dgetri_(N, mat_->ptr(), N, ipiv.ptr(), work, lwork, info);
	if(info){ std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl; }
	delete[] work;
}

template<>
void Lapack<std::complex<double> >::getri(Vector<int>& ipiv){
	unsigned int N(ipiv.size());
	unsigned int const lwork(3*N);
	std::complex<double>* work(new std::complex<double>[lwork]);
	int info(1);
	zgetri_(N, mat_->ptr(), N, ipiv.ptr(), work, lwork, info);

	if(info){ std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl; }
	delete[] work;
}
/*}*/

/*compute inverse of a matrix after using a Bunch-Kaufman factorization : dsytri*/
/*{*/
template<>
void Lapack<double>::sytri(Vector<int>& ipiv){
	int info(1);
	unsigned int N(mat_->col());
	double* work(new double[N]);
	dsytri_('U', N, mat_->ptr(), mat_->row(), ipiv.ptr(), work, info);

	if(info){
		std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl;
		ipiv.set();
	}
	delete[] work;
}

template<>
void Lapack<std::complex<double> >::sytri(Vector<int>& ipiv){
	(void)(ipiv);
	std::cerr<<__PRETTY_FUNCTION__<<" : not implemented for complex"<<std::endl;
}
/*}*/

/*compute the condition number : dgecon zgecon sycon*/
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
		std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl;
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
		std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl;
		return 0;
	} else {
		return rcond;
	}
}

template<>
double Lapack<double>::sycon(Vector<int> const& ipiv, double anorm){
	unsigned int N(mat_->col());
	int info(1);
	double rcond(0);
	double* work(new double[2*N]);
	int* iwork(new int[N]);
	dsycon_('U', N, mat_->ptr(), N, ipiv.ptr(), anorm, rcond, work, iwork, info);
	delete[] work;
	delete[] iwork;
	if(info){
		std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl;
		return 0;
	} else {
		return rcond;
	}
}

template<>
double Lapack<std::complex<double> >::sycon(Vector<int> const& ipiv, double anorm){
	(void)(anorm);
	(void)(ipiv);
	std::cerr<<__PRETTY_FUNCTION__<<" : not implemented for complex"<<std::endl;
	return 0;
}
/*}*/

/*compute the norm : dlange zlange dlansy*/
/*{*/
template<>
double Lapack<double>::lange(){
	double *work(NULL);
	return dlange_('1', mat_->row(), mat_->col(), mat_->ptr(), mat_->row(), work);
}

template<>
double Lapack<std::complex<double> >::lange(){
	double *work(NULL);
	return zlange_('1', mat_->row(), mat_->col(), mat_->ptr(), mat_->row(), work);
}

template<>
double Lapack<double>::lansy(){
	unsigned int N(mat_->row());
	double *work(new double[N]);
	double anorm(dlansy_('1', 'U', N, mat_->ptr(), N, work));
	delete[] work;
	return anorm;
}

template<>
double Lapack<std::complex<double> >::lansy(){
	std::cerr<<__PRETTY_FUNCTION__<<" : not implemented for complex"<<std::endl;
	return 0;
}
/*}*/

/*compute the eigensystem : syev heev geev*/
/*{*/
template<>
void Lapack<double>::syev(Vector<double>& EVal, char job){
	unsigned int N(mat_->row());
	EVal.set(N);
	int info(1);
	int lwork(-1);
	double wopt;

	dsyev_(job, 'U', N, mat_->ptr(), N, EVal.ptr(), &wopt, lwork, info);
	lwork = int(wopt);
	double* work(new double[lwork]);
	dsyev_(job, 'U', N, mat_->ptr(), N, EVal.ptr(), work, lwork, info);

	if(info){ std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl; }
	delete[] work;
}

template<>
void Lapack<std::complex<double> >::syev(Vector<double>& EVal, char job){
	(void)(EVal);
	(void)(job);
	std::cerr<<__PRETTY_FUNCTION__<<" : not implemented for complex"<<std::endl;
}

template<>
void Lapack<std::complex<double> >::heev(Vector<double>& EVal, char job){
	unsigned int N(mat_->row());
	EVal.set(N);
	int info(1);
	int lwork(-1);
	std::complex<double> wopt;
	double* rwork(new double[3*N-2]);

	zheev_(job, 'U', N, mat_->ptr(), N, EVal.ptr(), &wopt, lwork, rwork, info);
	lwork = int(wopt.real());
	std::complex<double>* work(new std::complex<double>[lwork]);
	zheev_(job, 'U', N, mat_->ptr(), N, EVal.ptr(), work, lwork, rwork, info);

	delete[] work;
	delete[] rwork;
	if(info){ std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl; }
}

template<>
void Lapack<double>::heev(Vector<double>& EVal, char job){
	(void)(EVal);
	(void)(job);
	std::cerr<<__PRETTY_FUNCTION__<<" : not implemented for real"<<std::endl;
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
			if(j==N-1 || !(wi[j]*wi[j+1]<0 && my::are_equal(wr[j],wr[j+1]))){
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
		std::cerr<<__PRETTY_FUNCTION__<<" : might be a proble because LU*RU != 1 but LU*A*RU = diag"<<std::endl;
		LEVec->set(N,N);
		for(unsigned int j(0);j<N;j++){
			if(j==N-1 || !(wi[j]*wi[j+1]<0 && my::are_equal(wr[j],wr[j+1]))){
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
	if(info){ std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl; }
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
		jobvr = 'V';
		ldvr = N;
		REVec->set(ldvr,N);
		vr = REVec->ptr();
	}
	if(LEVec){
		jobvl = 'V';
		ldvl = N;
		LEVec->set(ldvl,N);
		vl = LEVec->ptr();
	}
	EVal.set(N);

	zgeev_(jobvl, jobvr, N, mat_->ptr(), N, EVal.ptr(), vl, ldvl, vr, ldvr, &wopt, lwork, rwork, info);
	lwork = int(wopt.real());
	std::complex<double>* work(new std::complex<double>[N*N]);
	zgeev_(jobvl, jobvr, N, mat_->ptr(), N, EVal.ptr(), vl, ldvl, vr, ldvr, work, lwork, rwork, info);

	delete[] work;

	if(info){ std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl; }
}
/*}*/

/*compute the generalized eigensystem : dsygv*/
/*{*/
template<>
void Lapack<double>::sygv(Matrix<double>& B, Vector<double>& EVal){
	unsigned int N(mat_->row());
	EVal.set(N);
	double wopt;
	int lwork(-1);
	int info(1);
	char jobl('N');

	dsygv_(1,jobl, 'U', N, mat_->ptr(), N, B.ptr(), N, EVal.ptr(), &wopt, lwork, info);
	lwork = int(wopt);
	double* work(new double[lwork]);
	dsygv_(1,jobl, 'U', N, mat_->ptr(), N, B.ptr(), N, EVal.ptr(), work, lwork, info);

	delete[] work;
	if(info){
		std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl;
		EVal.set();
	}
}

template<>
void Lapack<std::complex<double> >::sygv(Matrix<std::complex<double> >& B, Vector<double>& EVal){
	(void)(B);
	(void)(EVal);
	std::cerr<<__PRETTY_FUNCTION__<<" : not implemented for complex"<<std::endl;
}
/*}*/

/*solve Ax=b : dposv*/
/*{*/
template<>
void Lapack<double>::posv(Vector<double>& b){
	unsigned int N(mat_->row());
	int info(1);
	dposv_('U',N,1,mat_->ptr(),N,b.ptr(),N,info);
	if(info){ std::cerr<<__PRETTY_FUNCTION__<<" : info="<<info<<std::endl; }
}

template<>
void Lapack<std::complex<double> >::posv(Vector<std::complex<double> >& b){
	(void)(b);
	std::cerr<<__PRETTY_FUNCTION__<<" : not implemented for complex"<<std::endl;
}
/*}*/
