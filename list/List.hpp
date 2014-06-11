#ifndef DEF_LIST
#define DEF_LIST

#include<iostream>

template<typename Type>
class List{
	public:
		List();
		List(Type const& t);
		List(List<Type> const& l);
		~List(){ clear(); }

		Type operator[](unsigned int idx) const;

		unsigned int size() const { return N_;}
		Type* get() const { return t_; }
		Type* last() const { return last_->t_; }

		void append(Type const& t);
		void pop();
		void remove(unsigned int idx);
		void clear();
		void swap(unsigned int a, unsigned int b);

		List<Type> sublist(unsigned int a, unsigned int b) const;
		void print(std::ostream& flux) const;

	protected:
		Type* t_;
		List<Type>* previous_;
		List<Type>* next_;
		List<Type>* last_;
		unsigned int N_;
};

template<typename Type>
std::ostream& operator<<(std::ostream& flux, List<Type> const& l);

/*constructors and destructor*/
/*{*/
template<typename Type>
List<Type>::List():
	t_(NULL),
	previous_(NULL),
	next_(NULL),
	last_(this),
	N_(0)
{}

template<typename Type>
List<Type>::List(Type const& t):
	t_(new Type(t)),
	previous_(NULL),
	next_(NULL),
	last_(this),
	N_(1)
{}

template<typename Type>
List<Type>::List(List<Type> const& l):
	t_(l.t_),
	previous_(NULL),
	next_(NULL),
	last_(this),
	N_(1)
{ 
	if(l.N_>1){
		List<Type> const* tmp(l.next_);
		while(tmp){
			append(*tmp->t_);
			tmp = tmp->next_;
		}
	}
}

/*}*/

/*operators*/
/*{*/
template<typename Type>
std::ostream& operator<<(std::ostream& flux, List<Type> const& l){
	l.print(flux);
	return flux;
}

template<typename Type>
Type List<Type>::operator[](unsigned int idx) const{
	List<Type> const* tmp(this);
	for(unsigned int i(0);i<idx;i++){ tmp = tmp->next_;}
	return (*tmp->t_);
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
void List<Type>::append(Type const& t){
	if(N_==0){ t_ = new Type(t); }
	else {
		last_->next_ = new List<Type>(t);
		last_->next_->previous_ = last_;		
		last_ = last_->next_;
	}
	N_++;
}

template<typename Type>
void List<Type>::pop(){
	switch(N_){
		case 0:{}break;
		case 1:{
				   last_ = this;
				   next_ = NULL;
				   delete t_;
				   t_ = NULL;
				   N_--;
			   }break;
		default:{
					List<Type>* tmp(last_);
					last_ = last_->previous_;
					last_->next_ = NULL;
					delete tmp;
					N_--;
				}
	}
}

template<typename Type>
void List<Type>::remove(unsigned int idx){
	if(idx==0){
		if(N_!=1){ 
			t_ = next_->t_; 
			List<Type>* tmp(next_);
			next_ = tmp->next_;
			tmp->next_ = NULL;
			tmp->t_ =  NULL;
			delete tmp;
			N_--;
		} else {
			clear();
		}
	} else {
		if(idx!=N_-1){
			List<Type>* tmp(this);
			for(unsigned int i(0);i<idx;i++){ tmp = tmp->next_; }
			List<Type>* next(tmp->next_);
			List<Type>* previous(tmp->previous_);
			tmp->next_ = NULL;
			delete tmp;
			previous->next_ = next;
			next->previous_ = previous;
			N_--;
		} else {
			pop();
		}
	}
}

template<typename Type>
void List<Type>::clear(){
	if(next_){ 
		delete next_;
		next_ = NULL;
	}
	if(t_){
		delete t_;
		t_ = NULL;
	}
	last_ = this;
	N_ = 0;
}

template<typename Type>
void List<Type>::swap(unsigned int a, unsigned int b){
	if(a>b){ swap(b,a); }
	if(a<b){
		List<Type>* list_a(this);
		for(unsigned int i(0);i<a;i++){ list_a = list_a->next_; }
		List<Type>* list_b(list_a);
		for(unsigned int i(0);i<b-a;i++){ list_b = list_b->next_; }
		Type* tmp(list_a->t_);
		list_a->t_ = list_b->t_; 
		list_b->t_ = tmp; 
	}
}
/*}*/

/*methods that return something*/
/*{*/
template<typename Type>
List<Type> List<Type>::sublist(unsigned int a, unsigned int b) const {
	List<Type> const* tmp(this);
	for(unsigned int i(0);i<a;i++){ tmp = tmp->next_; }
	List<Type> l(*(tmp->t_));
	for(unsigned int i(0);i<b-a;i++){
		tmp = tmp->next_; 
		l.append(*(tmp->t_));
	}
	return l;
}

template<typename Type>
void List<Type>::print(std::ostream& flux) const {
	List<Type> const* tmp(this);
	while(tmp && t_){
		flux<<(*tmp->t_)<<" "; 
		tmp = tmp->next_;
	}
}
/*}*/
#endif
