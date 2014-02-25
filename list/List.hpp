#ifndef DEF_List
#define DEF_List

#include<iostream>

template<typename Type>
class List{
	public:
		List();
		List(Type const& t);
		~List();

		Type operator[](unsigned int i) const;

		unsigned int size() const { return N_;}
		void append(Type const& t);

		void remove(unsigned int i);

	protected:
		List(Type const& t, unsigned int N);
		
		Type t_;
		List* next_;
		unsigned int N_;
};

template<typename Type>
std::ostream& operator<<(std::ostream& flux, List<Type> const& l);

/*constructors and destructor*/
/*{*/
template<typename Type>
List<Type>::List():
	next_(0),
	N_(0)
{ }

template<typename Type>
List<Type>::List(Type const& t):
	t_(t),
	next_(0),
	N_(1)
{ }

template<typename Type>
List<Type>::List(Type const& t, unsigned int N):
	t_(t),
	next_(0),
	N_(N)
{}

template<typename Type>
List<Type>::~List(){ 
	//std::cout<<"kill"<<t_<<" "<<N_<<std::endl;
	if(this->next_){
		delete this->next_;
		this->next_ = 0;
	}
}
/*}*/

/*operators*/
/*{*/
template<typename Type>
std::ostream& operator<<(std::ostream& flux, List<Type> const& l){
	for(unsigned int i(0); i<l.size();i++){
		flux<<l[i]<<" ";
	}
	return flux;
}


template<typename Type>
Type List<Type>::operator[](unsigned int i) const{
	if(i != 0){
		return (*next_)[--i];
	} else {
		return t_;
	}
}
/*}*/

/*methods*/
/*{*/
template<typename Type>
void List<Type>::append(Type const& t){
	if(N_== 0){
		t_ = t;
	} else {
		if(next_){
			next_->append(t);
		} else {
			next_ = new List(t,N_);
		}
	}
	N_++;
}

template<typename Type>
void List<Type>::remove(unsigned int i){
	N_--;
	if(i<2){
		if(this->next_){
			if(next_->next_){
				List* tmp(next_->next_);
				next_->next_ = 0;
				if(i==0){ this->t_ = next_->t_; }
				delete this->next_;
				this->next_ = tmp;
			} else {
				delete this->next_;
				this->next_ = NULL;
			}
		}
		/*!\warning the senario where i remove the first and only entry may not
		 * be correctly handled*/
	} else { 
		next_->remove(--i); 
	}

}
/*}*/

#endif
