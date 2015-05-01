#ifndef DEF_LIST
#define DEF_LIST

#include <memory>
#include "IOFiles.hpp"

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
 * to the first element of the list. The free_ pointer is free to point
 * anywhere. This construction allows the "free pointer" to point to any
 * element of the list and can therefore access quickly to the next or previous
 * element.
 *
 * The size() method should not be called too often because it needs to count
 * the number of element by moving trough the whole chain.
 */
/*}*/
template<typename Type>
class List{
	public:
		List() = default;
		List(List&&) = default;
		~List(){ set(); }
		/*{Forbidden*/
		List(List<Type> const& l) = delete;
		List& operator=(List) = delete;
		/*}*/

		void set();
		void set_free() const { free_ = const_cast<List<Type>* const>(this); }
		void set_free(List<Type>* free) const { free_ = free; }
		List<Type>* get_free() { return free_; }

		Type& get() { return (*free_->t_); }
		Type const& get() const { return (*free_->t_); }
		std::shared_ptr<Type> const& get_ptr() const { return free_->t_; }
		Type& first(){ return *next_->t_; }
		Type const& first() const { return *next_->t_; }
		Type& last() { return (*next_->free_->t_); }
		Type const& last() const { return (*next_->free_->t_); }

		void add_start(std::shared_ptr<Type> t);
		void add_end(std::shared_ptr<Type> t);
		void add(std::shared_ptr<Type> t, unsigned int const& idx);
		void add_sort(std::shared_ptr<Type> t,  std::function<bool (Type const&, Type const&)> cmp);
		void add_after_free(std::shared_ptr<Type> t);

		void pop_start();
		void pop_end();
		void pop(unsigned int const& idx);
		void pop_range(unsigned int const& a, unsigned int const& b);

		bool find_sorted(std::shared_ptr<Type> t, std::function<unsigned int (Type const&, Type const&)> cmp);
		void fuse_with_free(std::shared_ptr<Type> t, std::function<void (Type&, Type&)> fuse);

		unsigned int size() const;
		bool go_to_next() const;
		void swap(unsigned int const& a, unsigned int const& b);
		List<Type> sublist(unsigned int const& a, unsigned int const& b) const;

		void header_rst(std::string const& s, RST& rst) const;

	protected:
		List(std::shared_ptr<Type> t);

		std::shared_ptr<Type> t_  = NULL;
		mutable List<Type>* free_ = this;
		List<Type>* next_         = NULL;
};

/*constructors and destructor*/
/*{*/
template<typename Type>
List<Type>::List(std::shared_ptr<Type> t):
	t_(t),
	free_(this),
	next_(NULL)
{}

template<typename Type>
void List<Type>::set(){
	if(next_){ 
		delete next_;
		next_ = NULL;
	}
}
/*}*/

/*i/o methods*/
/*{*/
template<typename Type>
std::ostream& operator<<(std::ostream& flux, List<Type> const& l){
	flux<<"("<<l.size()<<")"; 
	while(l.go_to_next()){ flux<<" "<<l.get(); }
	return flux;
}

template<typename Type>
std::istream& operator>>(std::istream& flux, List<Type>& l){
	unsigned int size;
	flux>>size;
	Type tmp;
	while(size--){ 
		flux>>tmp;
		l.add_end( std::make_shared<Type>(tmp) ); 
	}
	return flux;
}

template<typename Type>
void List<Type>::header_rst(std::string const& s, RST& rst) const {
	rst.def(s,"List("+my::tostring(size())+")"); 
}

template<typename Type>
IOFiles& operator<<(IOFiles& w, List<Type> const& l){
	if(w.is_binary()){
		l.set_free();
		w<<l.size();
		while(l.go_to_next()){ w<<l.get(); }
	} else { w.stream()<<l; }
	return w;
}

template<typename Type>
IOFiles& operator>>(IOFiles& r, List<Type>& l){
	if(r.is_binary()){
		unsigned int size;
		r>>size;
		while(size--){ l.add_end(std::make_shared<Type>(r.read<Type>())); }
	} else { r.stream()>>l; }
	return r;
}
/*}*/

/*add an element to the list*/
/*{*/
template<typename Type>
void List<Type>::add_start(std::shared_ptr<Type> t){
	free_ = new List<Type>(t);
	if(next_){
		free_->free_ = next_->free_;
		next_->free_ = free_;
		free_->next_ = next_;
	} else {
		free_->free_ = free_;
	}
	next_ = free_;
	free_ = this;
}

template<typename Type>
void List<Type>::add_end(std::shared_ptr<Type> t){
	free_ = new List<Type>(t);
	if(next_){
		next_->free_->next_ = free_;
		free_->free_ = next_->free_;
		next_->free_ = free_;
	} else {
		free_->free_ = free_;
		next_ = free_;
	}
	free_ = this;
}

template<typename Type>
void List<Type>::add(std::shared_ptr<Type> t, unsigned int const& idx){
	if(idx == 0){ add_start(t); }
	else{
		free_ = next_;
		unsigned int i(0);
		while ( ++i<idx && free_ ) { free_ = free_->next_; }
		if(free_){
			if(free_->next_){
				free_->next_->free_ = new List<Type>(t);
				free_->next_->free_->next_ = free_->next_;
				free_->next_ = free_->next_->free_;
				free_->next_->free_ = free_;
				free_ = this;
			} else {
				add_end(t);
			}
		}
	}
}

template<typename Type>
void List<Type>::add_after_free(std::shared_ptr<Type> t){
	if(free_==this){ add_start(t); }
	else {
		if(free_->next_){
			free_->next_->free_ = new List<Type>(t);
			free_->next_->free_->next_ = free_->next_;
			free_->next_ = free_->next_->free_;
			free_->next_->free_ = free_;
			free_ = this;
		} else {
			add_end(t);
		}
	}
}

template<typename Type>
void List<Type>::add_sort(std::shared_ptr<Type> t, std::function<bool (Type const&, Type const&)> cmp){
	if(next_){
		if(next_->next_){
			if(cmp(*next_->free_->t_,*t)){ add_end(t); }
			else {
				free_ = this;
				while( cmp(*free_->next_->t_,*t) ){ free_ = free_->next_; }
				add_after_free(t);
			}
		} else {
			if( cmp(*next_->t_,*t) ){ add_end(t); } 
			else { add_start(t); }
		}
	} else { add_end(t); }
}
/*}*/

/*remove an element from the list*/
/*{*/
template<typename Type>
void List<Type>::pop_start(){
	if(next_){
		free_ = next_->next_;
		if(free_){
			next_->next_ = NULL;
			free_->free_ = next_->free_;
			delete next_;
			next_ = free_;
		} else {
			delete next_;
			next_ = NULL;
		}
	}
	free_ = this;
}

template<typename Type>
void List<Type>::pop_end(){
	if(next_){
		if(next_->next_){
			next_->free_ = next_->free_->free_;
			delete next_->free_->next_;
			next_->free_->next_ = NULL;
		} else {
			delete next_;
			next_ = NULL;
		}
	}
	free_ = this;
}

template<typename Type>
void List<Type>::pop(unsigned int const& idx){
	if(idx == 0){ pop_start(); }
	else {
		free_ = next_;
		unsigned int i(0);
		while ( i++<idx && free_ )
		{ free_ = free_->next_; }
		if(free_){
			if(free_->next_){
				free_->free_->next_ = free_->next_;
				free_->next_->free_ = free_->free_;
				free_->next_ = NULL;
				delete free_;
				free_ = this;
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
		std::swap(list_a->t_,list_b->t_); 
	}
}

template<typename Type>
void List<Type>::fuse_with_free(std::shared_ptr<Type> t, std::function<void (Type&, Type&)> fuse){
	fuse(*free_->t_,*t);
}
/*}*/

/*methods that return something*/
/*{*/
template<typename Type>
bool List<Type>::find_sorted(std::shared_ptr<Type> t, std::function<unsigned int (Type const&, Type const&)> cmp){
	if(next_){//check if there is at least one element
		if(next_->next_){//check if there is at least two elements
			switch(cmp(*next_->free_->t_,*t)){ //check cmp(last,t)
				case 0://the last element is "bigger" than t
					{ 
						switch(cmp(*next_->t_,*t)){
							case 0:
								{ 
									free_ = this;
									return false;
								}break;
							case 1:
								{
									free_ = this;
									unsigned int c(1);
									while(c==1){
										free_ = free_->next_;
										c = cmp(*free_->next_->t_,*t);
									}
									if(c==0){ return false; } 
									else { 
										free_ = free_->next_;
										return true; 
									}
								}break; 
							case 2:
								{
									free_ = next_;
									return true;
								}
						}
					} break;
				case 1://the last element is "smaller" than t
					{ 
						free_ = next_->free_;
						return false;
					} break;
				case 2://the last element is equal to t
					{ 
						free_ = next_->free_;
						return true;
					} break;
			}
		} else {
			switch(cmp(*next_->t_,*t)){//check the unique element of the list
				case 0://the unique element is "bigger"  than t
					{ 
						free_ = this;
						return false;
					} break;
				case 1://the unique element is "smaller" than t
					{ 
						free_ = next_;
						return false;
					} break;
				case 2://the unique element is equal to t 
					{ 
						free_ = next_;
						return true;
					} break;
			}
		}
		std::cerr<<"bool List<Type>::find_sorted(std::shared_ptr<Type> t, std::function<unsigned int (Type const&, Type const&)> cmp) : unexpected value returned from cmp"<<std::endl;
		return false;
	} else {
		free_ = this;
		return false;
	}
}

template<typename Type>
bool List<Type>::go_to_next() const { 
	free_ = free_->next_; 
	if(free_){ return true; }
	else { 
		free_ = const_cast<List<Type>* const>(this);
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
	free_ = next_;
	while(free_){ 
		N++;
		free_ = free_->next_;
	}
	free_ = const_cast<List<Type>* const>(this);
	return N;
}
/*}*/
#endif
