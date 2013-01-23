#ifndef DEF_MLAPACK
#define DEF_MLAPACK

#include "Matrice.hpp"
#include "Vecteur.hpp"

class Lapack{
	public:
/*Constructors and destructor*/
		Lapack(Matrice const& mat, char matrice_type);
		Lapack(double* m, unsigned int N, char matrice_type);
		~Lapack();

		Vecteur eigensystem() const;
		void lu(Matrice& L, Matrice& U) const;
		double det() const;
		void inv() const;

		void reset(); 
		Matrice LapToMat() const;
		
	private:
		Lapack();
		Lapack(Lapack const& l);

		double *m;
		double const *m_original;
		unsigned int const N;
		char const matrice_type; //implémenter les routines lapack selon G,S,H...
		bool const preserve_matrix;
		
/*private methods used to call lapack*/
		void dgetrf(int ipiv[]) const;
		void dgetri(int ipiv[]) const;
		void dsyev(Vecteur const& EVal) const;
};

//work for a symmetric matrix
extern "C" void dsyev_(char *jobz,
		char *uplo,
		unsigned int const *n,
		double *m,
		unsigned int const *lda, 
		double *w, 
		double work[], 
		unsigned int const *lwork, 
		int *info);

extern "C" void dgetrf_(unsigned int const *row,
		unsigned int const *rol,
		double *m,
		unsigned int const *lda,
		int ipiv[],
		int *info);

extern "C" void dgetri_(unsigned int const *n,
		double *m,
		unsigned int const *lda, 
		int ipiv[],
		double work[], 
		unsigned int const *lwork, 
		int *info);

/* //{ for a symmetric matrix, doesn't give me the lu decomposition
extern "C" void dsytrf_(char *uplo,
		unsigned int const *n,
		double *m,
		unsigned int const *lda,
		int ipiv[],
		double work[], 
		unsigned int const *lwork, 
		int *info);
	//}	*/
#endif
