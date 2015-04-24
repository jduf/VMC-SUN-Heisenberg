#ifndef DEF_MLAPACK
#define DEF_MLAPACK

#include "Vector.hpp"

/*{*/
extern "C" void dsyev_( /*eigensystem of a real symmetric matrix*/
		char const& jobz,
		char const& uplo,
		unsigned int const& n,
		double *m,
		unsigned int const& ldm,
		double *w,
		double *work,
		int const& lwork,
		int& info
		);
extern "C" void zheev_( /*eigensystem of a complex hermitian matrix*/
		char const& jobz,
		char const& uplo,
		unsigned int const& n,
		std::complex<double> *m,
		unsigned int const& ldm,
		double *w,
		std::complex<double> *work,
		int const& lwork,
		double *rwork,
		int& info
		);
extern "C" void dgeev_( /*eigensystem of a real general matrix*/
		char const& jobvl,
		char const& jobvr,
		unsigned int const& n,
		double *m,
		unsigned int const& ldm,
		double *wr,
		double *wi,
		double *vl,
		unsigned int const& ldvl,
		double *vr,
		unsigned int const& ldvr,
		double *work,
		int const& lwork,
		int& info
		);
extern "C" void zgeev_( /*eigensystem of a complex general matrix*/
		char const& jobvl,
		char const& jobvr,
		unsigned int const& n,
		std::complex<double> *m,
		unsigned int const& ldm,
		std::complex<double> *w,
		std::complex<double> *vl,
		unsigned int const& ldvl,
		std::complex<double> *vr,
		unsigned int const& ldvr,
		std::complex<double> *work,
		int const& lwork,
		double *rwork,
		int& info
		);

extern "C" void dgetrf_(
		unsigned int const& row,
		unsigned int const& col,
		double *m,
		unsigned int const& ldm,
		int *ipiv,
		int& info
		);
extern "C" void zgetrf_(
		unsigned int const& row,
		unsigned int const& col,
		std::complex<double> *m,
		unsigned int const& ldm,
		int *ipiv,
		int& info
		);

extern "C" void dgetri_(
		unsigned int const& n,
		double *m,
		unsigned int const& ldm,
		int *ipiv,
		double *work,
		int const& lwork,
		int& info
		);
extern "C" void zgetri_(
		unsigned int const& n,
		std::complex<double> *m,
		unsigned int const& ldm,
		int *ipiv,
		std::complex<double> *work,
		int const& lwork,
		int& info
		);

extern "C" void dgecon_(
		char const& norm,
		unsigned int const& n,
		double const *m, 
		unsigned int const& ldm,
		double const& anorm, 
		double& rcond, 
		double *work,
		int const *iwork,
		int& info
		);
extern "C" void zgecon_(
		char const& norm,
		unsigned int const& n,
		std::complex<double> const *m, 
		unsigned int const& ldm,
		double const& anorm, 
		double& rcond, 
		std::complex<double> *work,
		double *rwork,
		int& info
		);

extern "C" double dlange_(
		char const& norm,
		unsigned int const& row,
		unsigned int const& col,
		double const *m,
		unsigned int const& ldm,
		double *work
		);
extern "C" double zlange_(
		char const& norm,
		unsigned int const& row,
		unsigned int const& col,
		std::complex<double> const *m,
		unsigned int const& ldm,
		double *work
		);

/*!start the A=QRP factorisation, gives directly R and P*/
extern "C" void dgeqp3_(
		unsigned int const& row,
		unsigned int const& col,
		double *m,
		unsigned int const& ldm,
		int *jpvt,
		double *tau,
		double *work,
		int const& lwork,
		int& info
		);
/*!ends the A=QRP factorisation, extract the Q matrix*/
extern "C" void dorgqr_(
		unsigned int const& row,
		unsigned int const& col,
		unsigned int const& k,
		double *m,
		unsigned int const& ldm,
		double *tau,
		double *work,
		int const& lwork,
		int& info
		);
/*}*/

/*{Description*/
/*!Class that allows an easy use of the LAPACK routines
 * 
 * - compute determinant
 * - compute inverse
 * - compute eigensystem
 * - compute LU decomposition
 * - compute QR decomposition
 */
/*}*/
template<typename Type>
class Lapack{
	public:
		/*{Description*/
		/*!Constructor that copy the input matrix, if use_new_matrix=false, the
		 * input matrix will be modified*/
		/*}*/
		Lapack(Matrix<Type>& mat, bool use_new_matrix, char matrix_type);
		/*!Destructor (delete mat_ if use_new_matrix_==true)*/
		~Lapack();
		/*!Forbids constructors*/
		Lapack() = delete;
		Lapack(Lapack const&) = delete;
		Lapack(Lapack&&) = delete;

		/*{Description*/
		/*!Specialized routine that computes the eigensystem. If the
		 * eigenvectors are computed, they are stored in column in mat_,
		 * otherwise mat_ is overwritten*/
		/*}*/
		void eigensystem(Vector<double>& EVal, bool compute_EVec);
		/*{Description*/
		/*!Specialized routine that computes the eigensystem. If the
		 * eigenvectors are computed (EVec != NULL), they are stored in column
		 * in EVec. In any case, mat_ is overwritten*/
		/*}*/
		void eigensystem(Vector<std::complex<double> >& EVal, Matrix<std::complex<double> >* REVec=NULL, Matrix<std::complex<double> >* LEVec=NULL);
		/*!Computes the determinant*/
		Type det();
		/*!Computes the LU decomposition*/
		void lu(Matrix<Type>& L, Matrix<Type>& U);
		/*!Computes the QR decomposition*/
		void qr(Matrix<Type>& Q, Matrix<Type>& R, bool permutation=false);
		/*!Computes the inverse*/
		void inv();
		/*{Description*/
		/*!Computes the inverse when Vector<int> is_singular(double& rcn=0.0)
		 * has allready been called and if the condition number is not too
		 * small*/
		/*}*/
		void inv(Vector<int>& ipiv);
		/*!Computes the norm of the matrix*/
		double norm();
		/*!Checks if a matrix is singular, the smaller rcn is, the more the
		 * matrix is ill-defined*/
		Vector<int> is_singular(double& rcn=0.0);

		void set_mat(Matrix<Type>& mat) const { mat = *mat_; }

	protected:
		/*!Forbids assertion operator*/
		Lapack& operator=(Lapack const& l);

		char const matrix_type_; //!< type of matrix (symmetric, hermitian, general,...)
		bool use_new_matrix_; //!< false if the original matrix will be overwritten by the LAPACK routines
		Matrix<Type>* mat_; //!< pointer on a Matrix

		/*!Specialized subroutine that calls a LAPACK routine to compute the
		 * LU decomposition*/
		void getrf(Vector<int>& ipiv);
		/*!Specialized subroutine that calls a LAPACK routine to compute the
		 * inverse after the call of getrf*/
		void getri(Vector<int>& ipiv);
		/*!Specialized subroutine that calls a LAPACK routine to compute the
		 * QR decomposition*/
		void geqp3(double* tau, int* jpvt);
		/*!Specialized subroutine that calls a LAPACK routine to compute the Q
		 * matrix of the QR decomposition after the call of geqp3*/
		void gqr(unsigned int k, double* tau);
		/*!Specialized subroutine that calls a LAPACK routine to compute the
		 * eigensystem of a symmetric real matrix, if the eigenvectors are
		 * required, they are stored in column in mat_, otherwise mat_ is
		 * overwritten*/
		void syev(Vector<double>& EVal, char job);
		/*!Specialized subroutine that calls a LAPACK routine to compute the
		 * eigensystem of an hermitian complex matrix, if the eigenvectors are
		 * required, they are stored in column in mat_, otherwise mat_ is
		 * overwritten*/
		void heev(Vector<double>& EVal, char job); 
		/*!Specialized subroutine that calls a LAPACK routine to compute the
		 * eigensystem of a general matrix, if the eigenvectors are
		 * required, they are stored in column in EVec, In any case, mat_ is
		 * overwritten*/
		void geev(Vector<std::complex<double> >& EVal, Matrix<std::complex<double> >* REVec, Matrix<std::complex<double> >* LEVec);
		/*!Specialized subroutine that calls a LAPACK routine to compute the
		 * norm of an general real matrix*/
		double lange(); 
		/*!Specialized subroutine that calls a LAPACK routine to compute the
		 * condition number of an general real matrix*/
		double gecon(double anorm); 
};

/*Constructors and destructor*/
/*{*/
template<typename Type>
Lapack<Type>::Lapack(Matrix<Type>& mat, bool use_new_matrix, char matrix_type):
	matrix_type_(matrix_type),
	use_new_matrix_(use_new_matrix),
	mat_(use_new_matrix_?new Matrix<Type>(mat):&mat)
{
	if(matrix_type_ == 'S' && mat_->row() != mat_->col()){ std::cerr<<"Lapack::constructor() : the matix is not symmetric"<<std::endl; }
}

template<typename Type>
Lapack<Type>::~Lapack(){
	if(use_new_matrix_){ delete mat_; } 
}
/*}*/

/*methods used to call lapack*/
/*{*/
template<typename Type> 
Type Lapack<Type>::det(){
	if (matrix_type_ != 'G'){
		std::cerr<<"Lapack::det() : Matrix type "<<matrix_type_<<" not implemented"<<std::endl;
		std::cerr<<"                the only matrix type implemented is G"<<std::endl;
		return 0;
	} else {
		/*! \warning the determinant needs to be rewritten for rectangles matrices */
		Type d(1.0);
		unsigned int N(mat_->row());
		Vector<int> ipiv(N);
		getrf(ipiv);
		for(unsigned int i(0); i<N; i++){
			if(ipiv(i) != int(i+1)){ //! \note ipiv return value between [1,N]
				d *= - (*mat_)(i,i);
			} else {
				d *= (*mat_)(i,i);
			}
		}
		return d;
	}
}

template<typename Type>
void Lapack<Type>::inv(Vector<int>& ipiv){
	if (matrix_type_ != 'G'){
		std::cerr<<"Lapack::inv(Vector<int>& ipiv) : Matrix type "<<matrix_type_<<" not implemented"<<std::endl;
		std::cerr<<"                                 the only matrix type implemented is G"<<std::endl;
	} else {
		if(ipiv.ptr()){ getri(ipiv); }
		else { std::cerr<<"Lapack::inv(Vector<int>& ipiv) : the matrix is singular"<<std::endl;}
	}
}

template<typename Type>
void Lapack<Type>::inv(){
	if (mat_->row() != mat_->col()) {
		std::cerr<<"void Lapack<Type>::inv() : no inverse for a rectangular matrix"<<std::endl;
	} else {
		if (matrix_type_ != 'G'){
			std::cerr<<"void Lapack<Type>::inv() : Matrix type "<<matrix_type_<<" not implemented"<<std::endl;
			std::cerr<<"                           the only matrix type implemented is G"<<std::endl;
		} else {
			Vector<int> ipiv(std::min(mat_->row(),mat_->col()));
			getrf(ipiv);
			if(ipiv.ptr()){ getri(ipiv); }
			else { std::cerr<<"Lapack::inv() : illegal value"<<std::endl; }
		}
	}
}

template<typename Type>
void Lapack<Type>::lu(Matrix<Type>& L, Matrix<Type>& U){
	if (matrix_type_ != 'G'){
		std::cerr<<"Lapack::lu(Matrix<Type>& L, Matrix<Type>& U) : Matrix type "<<matrix_type_<<" not implemented"<<std::endl;
		std::cerr<<"                                               the only matrix type implemented is G"<<std::endl;
	} else {
		if (mat_->row()!=mat_->col()){
			std::cerr<<"Lapack::lu(Matrix<Type>& L, Matrix<Type>& U) : need to be implemented for rectangle matrix"<<std::endl;
		} else {
			std::cerr<<"Lapack::lu(Matrix<Type>& L, Matrix<Type>& U) : could return P instead of having a void method"<<std::endl;
			unsigned int N(mat_->row());
			L=Matrix<Type>(N,N);
			U=Matrix<Type>(N,N);
			Vector<int> ipiv(N);
			getrf(ipiv);
			for(unsigned int i(0); i<N; i++){
				L(i,i)=1.0;
				U(i,i)=(*mat_)(i,i);
				for(unsigned int j(i+1); j<N; j++){
					U(i,j) = (*mat_)(i,j);
					L(j,i) = (*mat_)(j,i);
				}
			}
		}
	}
}

template<typename Type>
void Lapack<Type>::qr(Matrix<Type>& Q, Matrix<Type>& R, bool permutation){
	if (matrix_type_ != 'G'){
		std::cerr<<"Lapack::qr(Matrix<Type>& Q, Matrix<Type>& R, bool permutation) : Matrix type "<<matrix_type_<<" not implemented"<<std::endl;
		std::cerr<<"                                                               : the only matrix type implemented is G"<<std::endl;
	} else {
		unsigned int k(std::min(mat_->row(),mat_->col()));
		double *tau(new double[k]);
		int *jpvt(new int[mat_->col()]); //! \warning contains element in [1,N]
		if(!permutation){
			for(unsigned int i(0); i<mat_->col();i++){
				jpvt[i] = 1;
			}
		}

		geqp3(tau,jpvt);
		R.set(k,mat_->col());
		for(unsigned int i(0);i<R.row();i++){
			for(unsigned int j(0);j<i;j++){
				R(i,j) = 0.0;
			}
			for(unsigned int j(i);j<R.col();j++){
				R(i,j) = (*mat_)(i,j);
			}
		}
		gqr(k,tau);

		if(mat_->col() > mat_->row()){
			Q.set(mat_->row(),mat_->row());
			for(unsigned int i(0);i<Q.row();i++){
				for(unsigned int j(0);j<Q.col();j++){
					Q(i,j) = (*mat_)(i,j);
				}
			}
		} else { 
			Q = *mat_; 
		}

		if(permutation){
			mat_->set(mat_->col(),mat_->col(),0.0);
			for(unsigned int i(0);i<mat_->col();i++){
				(*mat_)(i,jpvt[i]-1) = 1;
			}
		}

		delete[] tau;
		delete[] jpvt;
		std::cerr<<"Lapack::qr(Matrix<Type>& Q, Matrix<Type>& R, bool permutation) : delete jpvt may cause trouble"<<std::endl;
	}
}

template<typename Type>
Vector<int> Lapack<Type>::is_singular(double& rcn){
	Vector<int> ipiv;
	if (matrix_type_ != 'G'){
		std::cerr<<"Lapack::is_singular(double& rcn) : Matrix type "<<matrix_type_<<" not implemented"<<std::endl;
		std::cerr<<"                                   the only matrix type implemented is G"<<std::endl;
	} else if (mat_->row() != mat_->col()){
		std::cerr<<"Lapack::is_singular(double& rcn) : need to be checked for rectangle matrix"<<std::endl;
	} else {
		ipiv.set(std::min(mat_->row(),mat_->col()));
		double m_norm(lange());
		getrf(ipiv);
		rcn = gecon(m_norm);
		if(rcn<1e-5){ ipiv.set(); }//!set the output matrix to NULL pointer 
	}
	return ipiv;
}

template<typename Type> 
double Lapack<Type>::norm(){
	if (matrix_type_ != 'G'){
		std::cerr<<"Lapack::norm() : Matrix type "<<matrix_type_<<" not implemented"<<std::endl;
		std::cerr<<"                 the only matrix type implemented is G"<<std::endl;
		return 0;
	} else {
		return lange();
	}
}

template<typename Type>
void Lapack<Type>::eigensystem(Vector<double>& EVal, bool compute_EVec){
	if(mat_->row() != mat_->col()){std::cerr<<"Lapack::eigensystem(Vector<double>& EVal, bool compute_EVec) : matrix is not square"<<std::endl;}
	char jobvr('N');
	if(compute_EVec){jobvr = 'V';}
	switch(matrix_type_){
		case 'S':
			{ syev(EVal,jobvr); } break;
		case 'H':
			{ heev(EVal,jobvr); } break;
		default:
			{
				std::cerr<<"Lapack::eigensystem<double>(Vector<double>& EVal, bool compute_EVec) : Matrix type "<<matrix_type_<<" not implemented for real matrix"<<std::endl;
				std::cerr<<"                                                                       the only matrix type implemented are G,S"<<std::endl;
			} break;
	}
}

template<typename Type>
void Lapack<Type>::eigensystem(Vector<std::complex<double> >& EVal, Matrix<std::complex<double> >* REVec, Matrix<std::complex<double> >* LEVec){
	if(mat_->row() != mat_->col()){std::cerr<<"Lapack::eigensystem(Vector<std::complex<double> >& EVal, Matrix<std::complex<double> >* REVec, Matrix<std::complex<double> >* LEVec) : matrix is not square"<<std::endl;}
	switch(matrix_type_){
		case 'G':
			{
				geev(EVal,REVec,LEVec); 
				mat_->set();
			} break;
		default:
			{
				std::cerr<<"Lapack<Type>::eigensystem(Vector<std::complex<double> >& EVal, Matrix<std::complex<double> >* EVec)() : Matrix type "<<matrix_type_<<" not implemented for real matrix"<<std::endl;
				std::cerr<<"                                                                                                   () : the only matrix type implemented are G"<<std::endl;
			} break;
	}
}
/*}*/
#endif
