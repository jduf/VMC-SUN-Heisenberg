#ifndef DEF_BOSONIC
#define DEF_BOSONIC

#include "Matrix.hpp"
#include "Vector.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class Bosonic{
	public:
		/*!Creates a Bosonic without any parameters set*/
		Bosonic();
		/*!delete all the variables dynamically allocated*/
		virtual ~Bosonic();

		Bosonic<Type>* get_bosonic(){ return this; }

	protected:
		Matrix<unsigned int> nn_; //!< nn_(i,j):jth neighbour of the ith site
		Matrix<unsigned int> cc_;
		Matrix<double> nu_;
		Vector<unsigned int> sl_;
		Matrix<Type> omega_;
};

/*constructors and destructor*/
/*{*/
template<typename Type>
Bosonic<Type>::Bosonic()
{}

template<typename Type>
Bosonic<Type>::~Bosonic(){}
/*}*/
#endif

