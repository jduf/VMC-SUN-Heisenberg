#ifndef DEF_LIST
#define DEF_LIST

#include <iostream>
#include <cassert>
#include <algorithm>
#include <utility>

template<typename Type>
class List{
	public:
		List();
		List(Type const& t);
		List(List<Type> const& l);
		~List(){ clear(); }

		Type const& operator[](unsigned int idx) const;
		Type& operator[](unsigned int idx);
		Type& first(){ return *t_; }
		Type const& first() const { return *t_; }
		Type& last() { return *(previous_->t_); }
		Type const& last() const { return *(previous_->t_); }

		unsigned int size() const { return N_;}

		void add_start(Type const& t);
		void add_end(Type const& t);
		void add(Type const& t, unsigned int const& idx);
		template<typename Function>
			void add_sort(Type const& t, Function cmp);

		void pop_start();
		void pop_end();
		void pop(unsigned int const& idx);
		void pop_range(unsigned int const& a, unsigned int const& b);

		void clear();
		void swap(unsigned int const& a, unsigned int const& b);

		List<Type> sublist(unsigned int const& a, unsigned int const& b) const;
		void print(std::ostream& flux) const;

	protected:
		Type* t_;
		List<Type>* previous_;
		List<Type>* next_;
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
	N_(0)
{}

template<typename Type>
List<Type>::List(Type const& t):
	t_(new Type(t)),
	previous_(NULL),
	next_(NULL),
	N_(1)
{

}


template<typename Type>
List<Type>::List(List<Type> const& l):
	t_(l.t_),
	previous_(NULL),
	next_(NULL),
	N_(1)
{ 
	if(l.N_>1){
		List<Type> const* tmp(l.next_);
		while(tmp){
			add_end(*tmp->t_);
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
Type const& List<Type>::operator[](unsigned int idx) const{
	assert(idx<N_);
	List<Type> const* tmp(this);
	for(unsigned int i(0);i<idx;i++){ tmp = tmp->next_;}
	return (*tmp->t_);
}

template<typename Type>
Type& List<Type>::operator[](unsigned int idx) {
	assert(idx<N_);
	List<Type> const* tmp(this);
	for(unsigned int i(0);i<idx;i++){ tmp = tmp->next_;}
	return (*tmp->t_);
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
void List<Type>::add_start(Type const& t){
	switch(N_){
		case 0:
			{ t_ = new Type(t); }break;
		case 1:
			{ 
				next_ = new List<Type>(t);
				next_->previous_ = this;
				previous_ = next_;
			}break;
		default:
			{
				List<Type>* tmp(next_);
				next_ = new List<Type>(*t_);
				tmp->previous_ = next_;
				next_->next_ = tmp;
				/*not cool delete+new*/
				delete t_;
				t_ = new Type(t);
			}
	}
	N_++;
}

template<typename Type>
void List<Type>::add_end(Type const& t){
	switch(N_){
		case 0:
			{ t_ = new Type(t); }break;
		case 1:
			{ 
				next_ = new List<Type>(t);
				next_->previous_ = this;
				previous_ = next_;
			}break;
		default:
			{
				List<Type>* old_last(previous_);
				previous_ = new List<Type>(t);
				previous_->previous_ = old_last;
				old_last->next_ = previous_;
			}
	}
	N_++;
}

template<typename Type>
void List<Type>::add(Type const& t, unsigned int const& idx){
	if(idx == 0){ add_start(t); }
	else{
		if( idx == N_ ){ add_end(t); }
		else {
			List<Type>* tmp(this);
			for(unsigned int i(0);i<idx-1;i++){ tmp = tmp->next_; }
			List<Type>* old_idx(tmp->next_);
			tmp->next_ = new List<Type>(t);
			tmp->next_->next_ = old_idx;
			old_idx->previous_ = tmp->next_;
		}
	}
}

template<typename Type>
template<typename Function>
void List<Type>::add_sort(Type const& t, Function cmp){
	List<Type> const* tmp(this);
	unsigned int i(0);
	while( tmp && cmp(*tmp->t_,t)){ 
		tmp = tmp->next_; 
		i++;
	}
	add(t,i);
}

template<typename Type>
void List<Type>::pop_end(){
	switch(N_){
		case 0:
			{}break;
		case 1:
			{
				next_ = NULL;
				delete t_;
				t_ = NULL;
				N_--;
			}break;
		default:
			{
				List<Type>* tmp(previous_);
				previous_ = previous_->previous_;
				previous_->next_ = NULL;
				delete tmp;
				N_--;
			}
	}
}

template<typename Type>
void List<Type>::pop_start(){
	if(N_!=1){ 
		List<Type>* tmp(next_);
		t_ = tmp->t_; 
		next_ = tmp->next_;
		if(N_!=2){ tmp->next_->previous_ = this; }
		tmp->next_ = NULL;
		tmp->t_ =  NULL;
		delete tmp;
		N_--;
	} else {
		clear();
	}
}

template<typename Type>
void List<Type>::pop(unsigned int const& idx){

	if(idx == 0){ pop_start(); }
	else {
		if(idx == N_-1) { pop_end(); }
		else {
			List<Type>* tmp(this);
			for(unsigned int i(0);i<idx;i++){ tmp = tmp->next_; }
			List<Type>* next(tmp->next_);
			List<Type>* previous(tmp->previous_);
			tmp->next_ = NULL;
			delete tmp;
			previous->next_ = next;
			next->previous_ = previous;
			N_--;
		} 
	}
}

template<typename Type>
void List<Type>::pop_range(unsigned int const& a, unsigned int const& b){
	if(a==b){ pop(a); }
	if(a>b){ pop_range(b,a); }
	if(a<b){
		if(a!=0){
			List<Type>* list_a(this);
			for(unsigned int i(0);i<a-1;i++){ list_a = list_a->next_; }
			List<Type>* tmp(list_a->next_);
			List<Type>* list_b(list_a);
			for(unsigned int i(0);i<b-a-1;i++){ list_b = list_b->next_; }
			list_a->next_ = list_b->next_;
			list_b->next_ = NULL;
			delete tmp;
		} else {
			pop_range(a+1,b);
			pop_start();
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
	N_ = 0;
}

template<typename Type>
void List<Type>::swap(unsigned int const& a, unsigned int const& b){
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
List<Type> List<Type>::sublist(unsigned int const& a, unsigned int const& b) const {
	List<Type> const* tmp(this);
	for(unsigned int i(0);i<a;i++){ tmp = tmp->next_; }
	List<Type> l(*(tmp->t_));
	for(unsigned int i(0);i<b-a;i++){
		tmp = tmp->next_; 
		l.add_end(*(tmp->t_));
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
