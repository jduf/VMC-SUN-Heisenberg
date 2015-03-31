#ifndef DEF_LIST
#define DEF_LIST

#include <iostream>
#include <functional>
#include <memory>

/*{Description*/
/*! The list is constructed according this scheme
\verbatim
    _______________________
   |   *t_  *t_  *t_  *t_  |
   |    ↑    ↑    ↑    ↑   |
   | o─→o←──→o←──→o←──→o─→x|
   | │  │              ↑   |
   | │  └──────────────┘   |
   | └──"free pointer"     |
   |_______________________|
\endverbatim
 *
 * The first element is particular because t_=NULL and the next_ pointer points
 * to the first element of the list. The move_ pointer is free. This
 * construction allows the "free pointer" to point to any element of the list
 * and can therefore access quickly to the next or previous element.
 *
 * The size() method should not be called too often because it needs to count
 * the number of element by moving trough the whole chain.
 */
/*}*/
template<typename Type>
class List{
	public:
		List();
		List(List<Type> const& l);
		~List(){ set(); }
		void set();

		Type& get() { return (*move_->t_); }
		Type const& get() const { return (*move_->t_); }
		std::shared_ptr<Type> const& get_ptr() const { return move_->t_; }
		Type& first(){ return *next_->t_; }
		Type const& first() const { return *next_->t_; }
		Type& last() { return (*next_->move_->t_); }
		Type const& last() const { return (*next_->move_->t_); }

		void add_start(std::shared_ptr<Type> t);
		void add_end(std::shared_ptr<Type> t);
		void add(std::shared_ptr<Type> t, unsigned int const& idx);
		void add_sort(std::shared_ptr<Type> t,  std::function<bool (Type const&, Type const&)> cmp);
		void add_after_move(std::shared_ptr<Type> t);

		void pop_start();
		void pop_end();
		void pop(unsigned int const& idx);
		void pop_range(unsigned int const& a, unsigned int const& b);


		bool find(std::shared_ptr<Type> t, std::function<unsigned int (Type const&, Type const&)> cmp);
		void fuse_with_move(std::shared_ptr<Type> t, std::function<void (Type&, Type&)> fuse);

		unsigned int size() const;
		bool move_forward() const;
		void swap(unsigned int const& a, unsigned int const& b);
		List<Type> sublist(unsigned int const& a, unsigned int const& b) const;
		void print(std::ostream& flux) const;

	protected:
		List(std::shared_ptr<Type> t);

		std::shared_ptr<Type> t_;
		mutable List<Type>* move_;
		List<Type>* next_;
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
	move_(NULL),
	next_(NULL)
{}

template<typename Type>
List<Type>::List(std::shared_ptr<Type> t):
	t_(t),
	move_(NULL),
	next_(NULL)
{}


template<typename Type>
List<Type>::List(List<Type> const& l):
	t_(l.t_),
	move_(l.move_),
	next_(l.next_)
{ 
	List<Type> const* tmp(l.next_);
	while(tmp){
		add_end(tmp->t_);
		tmp = tmp->next_;
	}
}

template<typename Type>
void List<Type>::set(){
	if(next_){ 
		delete next_;
		next_ = NULL;
	}
}
/*}*/

/*add an element to the list*/
/*{*/
template<typename Type>
void List<Type>::add_start(std::shared_ptr<Type> t){
	if(next_){
		List<Type>* tmp(new List<Type>(t));
		if(move_ == next_){ move_ = tmp; }
		tmp->move_ = next_->move_;
		tmp->next_ = next_;
		next_->move_ = tmp;
		next_ = tmp;
	} else {
		next_ = new List<Type>(t);
		next_->move_ = next_;
		move_ = next_;
	}
}

template<typename Type>
void List<Type>::add_end(std::shared_ptr<Type> t){
	if(next_){
		next_->move_->next_ = new List<Type>(t);
		next_->move_->next_->move_ = next_->move_;
		next_->move_ = next_->move_->next_;
	} else {
		next_ = new List<Type>(t);
		next_->move_ = next_;
		move_ = next_;
	}
}

template<typename Type>
void List<Type>::add(std::shared_ptr<Type> t, unsigned int const& idx){
	if(idx == 0){ add_start(t); }
	else{
		List<Type>* tmp(next_);
		unsigned int i(0);
		while ( ++i<idx && tmp )
		{ tmp = tmp->next_; }
		if( tmp ){
			if(tmp->next_){
				tmp->next_->move_ = new List<Type>(t);
				tmp->next_->move_->next_ = tmp->next_;
				tmp->next_ = tmp->next_->move_;
				tmp->next_->move_ = tmp;
			} else {
				add_end(t);
			}
		}
	}
}

template<typename Type>
void List<Type>::add_sort(std::shared_ptr<Type> t, std::function<bool (Type const&, Type const&)> cmp){
	if(next_){
		if(next_->next_){
			if(cmp(*next_->move_->t_,*t)){ add_end(t); }
			else {
				List<Type> const* tmp(next_);
				unsigned int i(0);
				while( cmp(*tmp->t_,*t) ){ /*removed the condition tmp->next_*/
					tmp = tmp->next_; 
					i++;
				}
				add(t,i);
			}
		} else {
			if( cmp(*next_->t_,*t) ){ add_end(t); } 
			else { add_start(t); }
		}
	} else {
		add_end(t);
	}
}

template<typename Type>
void List<Type>::add_after_move(std::shared_ptr<Type> t){
	if(move_==this){ add_start(t); }
	else {
		if(move_->next_){
			move_->next_->move_ = new List<Type>(t);
			move_->next_->move_->next_ = move_->next_;
			move_->next_= move_->next_->move_;
			move_->next_->move_ = move_;
		} else {
			add_end(t);
		}
	}
	move_ = next_;
}
/*}*/

/*remove an element from the list*/
/*{*/
template<typename Type>
void List<Type>::pop_start(){
	if(next_){
		if(next_->next_){
			if(move_ == next_){ move_ = next_->next_; }
			List<Type>* tmp(next_->next_);
			tmp->move_ = next_->move_;
			next_->next_ = NULL;
			delete next_;
			next_ = tmp;
		} else {
			delete next_;
			next_ = NULL;
			move_ = NULL;
		}
	}
}

template<typename Type>
void List<Type>::pop_end(){
	if(next_){
		if(next_->next_){
			next_->move_ = next_->move_->move_;
			delete next_->move_->next_;
			next_->move_->next_ = NULL;
		} else {
			delete next_;
			next_ = NULL;
			move_ = NULL;
		}
	}
}

template<typename Type>
void List<Type>::pop(unsigned int const& idx){
	if(idx == 0){ pop_start(); }
	else {
		List<Type>* tmp(next_);
		unsigned int i(0);
		while ( i++<idx && tmp )
		{ tmp = tmp->next_; }
		if(tmp){
			if(tmp->next_){
				tmp->move_->next_ = tmp->next_;
				tmp->next_->move_ = tmp->move_;
				tmp->next_ = NULL;
				delete tmp;
			} else {
				pop_end();
			}
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
		List<Type>* list_a(next_);
		for(unsigned int i(0);i<a;i++){ list_a = list_a->next_; }
		List<Type>* list_b(list_a);
		for(unsigned int i(0);i<b-a;i++){ list_b = list_b->next_; }
		std::shared_ptr<Type> tmp(list_a->t_);
		list_a->t_ = list_b->t_; 
		list_b->t_ = tmp; 
	}
}

template<typename Type>
void List<Type>::fuse_with_move(std::shared_ptr<Type> t, std::function<void (Type&, Type&)> fuse){
	fuse(*move_->t_,*t);
}
/*}*/

/*methods that return something*/
/*{*/
template<typename Type>
bool List<Type>::find(std::shared_ptr<Type> t, std::function<unsigned int (Type const&, Type const&)> cmp){
	if(next_){
		if(next_->next_){
			if(cmp(*next_->move_->t_,*t)==2){
				std::cout<<"A"<<*t<<std::endl;
				move_ = next_->move_;
				return true;
			} else {
				move_=next_;
				unsigned int c(cmp(*move_->t_,*t));
				while(c==1){
					move_ = move_->next_;
					c = cmp(*move_->t_,*t);
				}
				if(c==0){
					std::cout<<"B"<<*t<<std::endl;
					return false; 
				}
				else{
					std::cout<<"C"<<*t<<std::endl;
					return true; 
				}
			}
		} else {
			switch(cmp(*next_->t_,*t)){
				case 0:
					{ 
						std::cout<<"D"<<*t<<std::endl;
						move_ = this;
						return false;
					} break;
				case 1: 
					{ 
						std::cout<<"E"<<*t<<std::endl;
						move_ = next_;
						return false;
					} break;
				default: 
					{ 
						std::cout<<"F"<<*t<<std::endl;
						return true;
					} break;
			}
		}
	} else {
		std::cout<<"G"<<*t<<std::endl;
		move_ = this;
		return false;
	}
}

template<typename Type>
bool List<Type>::move_forward() const { 
	move_ = move_->next_; 
	if(move_){ return true; }
	else { 
		move_ = next_;
		return false; 
	}
}

template<typename Type>
List<Type> List<Type>::sublist(unsigned int const& a, unsigned int const& b) const {
	List<Type> const* tmp(this);
	for(unsigned int i(0);i<a;i++){ tmp = tmp->next_; }
	List<Type> l(std::make_shared<Type>(*tmp->t_));
	for(unsigned int i(0);i<b-a;i++){
		tmp = tmp->next_; 
		l.add_end(std::make_shared<Type>(*tmp->t_));
	}
	return l;
}

template<typename Type>
unsigned int List<Type>::size() const{
	unsigned int N(0);
	List<Type> const* tmp(next_);
	while(tmp){ 
		N++;
		tmp = tmp->next_;
	}
	return N;
}

template<typename Type>
void List<Type>::print(std::ostream& flux) const {
	List<Type> const* tmp(next_);
	while(tmp){
		flux<<(*tmp->t_)<<" "; 
		tmp = tmp->next_;
	}
	flux<<"("<<size()<<")"; 
}
/*}*/
#endif
