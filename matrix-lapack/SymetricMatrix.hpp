#ifndef DEF_SYMETRICMATRIX
#define DEF_SYMETRICMATRIX

#include "Matrix.hpp"

class SymetricMatrix: public Matrix<double>{
	public:
		/*!Default constructor that initializes *m to NULL and N to 0*/
		SymetricMatrix() = default;
		/*!Initializes a static array of double of size N_row*N_col*/
		SymetricMatrix(unsigned int N_row);
		/*!Initializes a static array of double of size N_row*N_col to a value val*/
		SymetricMatrix(unsigned int N_row, double val);
		/*!Constructor that reads from file*/
		SymetricMatrix(IOFiles& r);
		/*!Deep copy constructor*/
		SymetricMatrix(SymetricMatrix const& mat) = default;
		/*!Move constructor (can't be default because of mat_)*/
		SymetricMatrix(SymetricMatrix&& mat) = default;
		/*!Delete the static array*/
		~SymetricMatrix() = default;

		/*!Accesses the (i,j)th entry of the matrix*/
		double const& operator()(unsigned int const& i, unsigned int const& j) const { 
			assert(i<row_ && j<col_); 
			if(i<j){ return mat_[i+j*(j-1)/2]; }
			else   { return mat_[j+i*(i-1)/2]; }
		}
		/*!Sets the (i,j)th entry of the matrix*/
		double& operator()(unsigned int const& i, unsigned int const& j){
			assert(i<row_ && j<col_); 
			if(i<j){ return mat_[i+j*(j-1)/2]; }
			else   { return mat_[j+i*(i-1)/2]; }
		}
};
#endif
