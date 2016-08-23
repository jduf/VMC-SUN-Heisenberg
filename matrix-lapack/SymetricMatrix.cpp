#include "SymetricMatrix.hpp"

/*constructors and destructor*/
/*{*/
SymetricMatrix::SymetricMatrix(unsigned int N){
	row_=N;
	col_=N;
	size_=N*(N+1)/2;
	mat_=size_?new double[size_]:NULL;
} 

SymetricMatrix::SymetricMatrix(unsigned int N, double val){
	row_ = N;
	col_ = N;
	size_= N*(N+1)/2;
	mat_ = size_?new double[size_]:NULL;
	for(unsigned int i(0);i<size_;i++){ mat_[i] = val; }
}

SymetricMatrix::SymetricMatrix(IOFiles& r){
	col_ = row_ = r.read<unsigned int>();
	size_= row_*col_;
	mat_ = size_?new double[size_]:NULL;
	r.read(mat_,size_,sizeof(double));
}
/*}*/
