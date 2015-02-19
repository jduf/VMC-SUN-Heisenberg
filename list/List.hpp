#ifndef DEF_LIST
#define DEF_LIST

#include <iostream>
#include <cassert>
#include <functional>

template<typename Type>
class List{
	public:
		List();
		List(Type* t);
		List(List<Type> const& l);
		~List(){ clear(); }
		void clear();

		Type const& operator[](unsigned int idx) const;
		Type& operator[](unsigned int idx);
		Type& first(){ return *t_; }
		Type const& first() const { return *t_; }
		Type& last() { return (*prev_->t_); }
		Type const& last() const { return (*prev_->t_); }

		unsigned int size() const { return N_;}

		void add_start(Type* t);
		void add_end(Type* t);
		void add(Type* t, unsigned int const& idx);
		void add_sort(Type* t,  std::function<bool (Type*,Type*)> cmp);
		void add_or_fuse_sort(Type* t, std::function<unsigned int (Type*,Type*)> cmp, std::function<void (Type*,Type*)> fuse);

		void pop_start();
		void pop_end();
		void pop(unsigned int const& idx);
		void pop_range(unsigned int const& a, unsigned int const& b);

		List<Type> sublist(unsigned int const& a, unsigned int const& b) const;
		void swap(unsigned int const& a, unsigned int const& b);
		void print(std::ostream& flux) const;

	protected:
		Type* t_;
		List<Type>* prev_;
		List<Type>* next_;
		unsigned int N_;
};

template<typename Type>
std::ostream& operator<<(std::ostream& flux, List<Type> const& l){
	l.print(flux);
	return flux;
}

/*constructors and destructor*/
/*{*/
template<typename Type>
List<Type>::List():
	t_(NULL),
	prev_(this),
	next_(NULL),
	N_(0)
{}

template<typename Type>
List<Type>::List(Type* t):
	t_(t),
	prev_(this),
	next_(NULL),
	N_(1)
{}


template<typename Type>
List<Type>::List(List<Type> const& l):
	t_(l.t_),
	prev_(l.prev_),
	next_(l.next_),
	N_(1)
{ 
	if(l.N_>1){
		List<Type> const* tmp(l.next_);
		while(tmp){
			add_end(tmp->t_);
			tmp = tmp->next_;
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
/*}*/

/*operators*/
/*{*/
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

/*add an element to the list*/
/*{*/
template<typename Type>
void List<Type>::add_start(Type* t){
	switch(N_){
		case 0:
			{ t_ = t; }break;
		case 1:
			{ 
				next_ = new List<Type>(t_);
				next_->prev_ = this;
				prev_ = next_;
				t_ = t;
			}break;
		default:
			{
				next_->prev_ = new List<Type>(t_);
				next_->prev_->next_ = next_;
				next_= next_->prev_;
				next_->prev_ = this;
				t_ = t;
			}
	}
	N_++;
}

template<typename Type>
void List<Type>::add_end(Type* t){
	switch(N_){
		case 0:
			{ t_ = t; }break;
		case 1:
			{
				next_ = new List<Type>(t);
				next_->prev_ = this;
				prev_ = next_;
			}break;
		default:
			{
				prev_->next_ = new List<Type>(t);
				prev_->next_->prev_ = prev_;
				prev_ = prev_->next_;
			}
	}
	N_++;
}

template<typename Type>
void List<Type>::add(Type* t, unsigned int const& idx){
	assert(idx<=N_);
	if(idx == 0){ add_start(t); }
	else{
		if( idx == N_ ){ add_end(t); }
		else {
			List<Type>* tmp(this);
			for(unsigned int i(0);i<idx-1;i++){ tmp = tmp->next_; }
			tmp->next_->prev_ = new List<Type>(t);
			tmp->next_->prev_->next_ = tmp->next_;
			tmp->next_ = tmp->next_->prev_;
			tmp->next_->prev_ = tmp;
			N_++;
		}
	}
}

template<typename Type>
void List<Type>::add_sort(Type* t, std::function<bool (Type*,Type*)> cmp){
	switch(N_){
		case 0:
			{
				t_ = t;
				N_++;
			}break;
		case 1:
			{
				if( cmp(t_,t) ){ add_end(t); } 
				else { add_start(t); }
			}break;
		default:
			{
				if(cmp(prev_->t_,t)){ add_end(t); }
				else {
					List<Type> const* tmp(this);
					unsigned int i(0);
					while( cmp(tmp->t_,t) ){ /*removed the condition tmp->next_*/
						tmp = tmp->next_; 
						i++;
					}
					add(t,i);
				}
			}
	}
}

template<typename Type>
void List<Type>::add_or_fuse_sort(Type* t, std::function<unsigned int (Type*,Type*)> cmp, std::function<void (Type*,Type*)> fuse){
	switch(N_){
		case 0:
			{
				t_ = t;
				N_++;
			}break;
		case 1:
			{
				switch(cmp(t_,t)){
					case 0: { add_start(t); } break;
					case 1: { add_end(t); } break;
					case 2: { fuse(t_,t); } break;
				}
			}break;
		default:
			{
				switch(cmp(prev_->t_,t)){
					case 0:
						{
							List<Type> const* tmp(this);
							unsigned int i(0);
							unsigned int c(cmp(tmp->t_,t));
							while(c==1){
								tmp = tmp->next_; 
								c = cmp(tmp->t_,t);
								i++;
							}
							if(c==0){ add(t,i); }
							else{ fuse(tmp->t_,t); }
						}break;
					case 1:
						{ add_end(t); }break;
					case 2:
						{ fuse(prev_->t_,t); }break;
				}
			}
	}
}
/*}*/

/*remove an element from the list*/
/*{*/
template<typename Type>
void List<Type>::pop_start(){
	switch(N_){
		case 0:
			{}break;
		case 1:
			{clear();}break;
		case 2:
			{
				delete t_;
				t_ = next_->t_; 
				next_->t_ =  NULL;
				delete next_;
				next_ = NULL;
				N_--;
			}break;
		default:
			{
				delete t_;
				t_ = next_->t_; 
				next_->t_ =  NULL;
				next_ = next_->next_;
				next_->prev_->next_ = NULL;
				delete next_->prev_;
				next_->prev_ = this;
				N_--;
			}break;
	}
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
				List<Type>* tmp(prev_);
				prev_ = prev_->prev_;
				prev_->next_ = NULL;
				delete tmp;
				N_--;
			}break;
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
			List<Type>* previous(tmp->prev_);
			tmp->next_ = NULL;
			delete tmp;
			previous->next_ = next;
			next->prev_ = previous;
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
/*}*/

/*methods that modifies the list*/
/*{*/
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
	List<Type> l(new Type(*tmp->t_));
	for(unsigned int i(0);i<b-a;i++){
		tmp = tmp->next_; 
		l.add_end(new Type(*tmp->t_));
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
	flux<<"("<<N_<<")";
}
/*}*/
#endif
