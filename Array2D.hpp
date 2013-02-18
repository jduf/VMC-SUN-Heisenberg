#ifndef DEF_ARRAY2D
#define DEF_ARRAY2D

#include <iostream>
#include <cassert>

template<typename A>
class Array2D{
	public:
		/*!Default constructor that initializes *m to NULL and N to 0*/
		Array2D();
		/*!Initializes a static array of M of size N*N to a value val*/
		Array2D(unsigned int N_row, unsigned int N_col);
		/*!Initializes a static array of M of size N*N to a value val*/
		Array2D(unsigned int N_row, unsigned int N_col, A val);
		/*!Deep copy*/
		Array2D(Array2D const& ar);
		/*!Delete the static array*/
		~Array2D();

		/*!Does a deep copie*/
		Array2D<A>& operator=(Array2D<A> const& arr);
		/*!Accesses the (i,j)th entry of the vector*/
		inline A const& operator()(unsigned int const& i, unsigned int const& j) const { assert(i<N_row); assert(j<N_col); return a[i+j*N_row]; };
		/*!Sets the (i,j)th entry of the vector*/
		inline A& operator()(unsigned int const& i, unsigned int const& j) { assert(i<N_row); assert(j<N_col); return a[i+j*N_row]; };
		/*!Additions this Array2D with another*/

		/*!Returns the pointer to the array*/
		inline A* ptr() const { return a; };
		/*!Returns the size of the arrix*/
		inline unsigned int col() const {return N_col;};
		inline unsigned int row() const {return N_row;};

	private:
		A *a; //!< pointer to a static array of the form m = [[ligne0],[ligne1],...]
		unsigned int N_row; //!< number of line
		unsigned int N_col; //!< number of col

		/*methods that modify the class*/
		void fill_Array2D(A val);
};

template<typename A>
std::ostream& operator<<(std::ostream& flux, Array2D<A> const& arr);
template<typename A>
std::istream& operator>>(std::istream& flux, Array2D<A>& arr);

/*Constructors and destructor*/
/*{*/
template<typename A>
Array2D<A>::Array2D():
	a(NULL),
	N_row(0),
	N_col(0)
{ }

template<typename A>
Array2D<A>::Array2D(unsigned int N_row, unsigned int N_col):
	a(new A[N_row*N_col]),
	N_row(N_row),
	N_col(N_col)
{ }

template<typename A>
Array2D<A>::Array2D(unsigned int N_row, unsigned int N_col, A val):
	a(new A[N_row*N_col]),
	N_row(N_row),
	N_col(N_col)
{
	fill_Array2D(val);
}

template<typename A>
Array2D<A>::Array2D(Array2D<A> const& arr):
	a(new A[arr.row()*arr.col()]),
	N_row(arr.row()),
	N_col(arr.col())
{
	for(unsigned int i(0);i<N_row*N_col;i++){
		a[i] = arr.a[i];
	}
}

template<typename A>
Array2D<A>::~Array2D(){
	delete[] a;
}
/*}*/

/*operators*/
/*{*/
template<typename A>
Array2D<A>& Array2D<A>::operator=(Array2D<A> const& arr){
	if(this->N_row!=arr.N_row || this->N_col!=arr.N_col ){
		delete[] this->a;
		this->a = new A[arr.N_row*arr.N_col];
		this->N_row = arr.N_row;
		this->N_col = arr.N_col;
	}
	for(unsigned int i(0); i<this->N_row*this->N_col; i++){
		this->a[i] = arr.a[i];
	}
	return (*this);
}

template<typename A>
std::ostream& operator<<(std::ostream& flux, Array2D<A> const& arr){
	for(unsigned int i(0);i<arr.row();i++){
		for(unsigned int j(0);j<arr.col();j++){
			flux<<arr(i,j)<<" ";
		}
		flux<<std::endl;
	}
	return flux;
}

template<typename A>
std::istream& operator>>(std::istream& flux, Array2D<A>& arr){
	for(unsigned int i(0);i<arr.row();i++){
		for(unsigned int j(0);j<arr.col();j++){
			flux>>arr(i,j);
		}
	}
	return flux;
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename A>
void Array2D<A>::fill_Array2D(A val){
	for(unsigned int i(0);i<N_row*N_col;i++){
		a[i]=val;
	}
}
/*}*/

#endif
