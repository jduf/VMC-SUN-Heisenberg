#ifndef DEF_CONTAINER
#define DEF_CONTAINER

#include "Miscellaneous.hpp"

/*{Variable*/
/*!Base class that will never be instantiate. It provides a way to create
 * GenericVariable with defferent type.*/
class Variable{
	public:
		/*!Constructor*/
		Variable(std::string const& name):name_(name){}
		/*!Destructor*/
		virtual ~Variable() = default;
		/*{Forbidden*/
		Variable() = delete;
		Variable(Variable&&) = delete;
		Variable& operator=(Variable const&) = delete;
		/*}*/

		/*!Method that allows a copy of the derived class*/
		virtual Variable* clone() const = 0;
		/*!Returns the name of the data*/
		std::string const& get_name() const { return name_;}

	protected:
		/*!Default copy constructor*/
		Variable(Variable const&) = default;
		
		std::string name_;//!< Name of the (Generic)Variable 
};
/*}*/

/*{GenericVariable*/
/*!Derived class that can contain any one value of any type and of name
 * Variable::name_*/
template<typename Type>
class GenericVariable : public Variable{
	public:
		/*!Constructor*/
		GenericVariable(std::string const& name, Type t):Variable(name),t_(t){}
		/*!Default destructor*/
		~GenericVariable() = default;
		/*{Forbidden*/
		GenericVariable() = delete;
		GenericVariable(GenericVariable&&) = delete;
		GenericVariable& operator=(GenericVariable) = delete;
		/*}*/
		/*!Returns a copy of*/
		GenericVariable<Type>* clone() const { return new GenericVariable<Type>(*this); }

		/*!Returns the value of the data*/
		Type const& get_val() const { return t_; }
		/*!Set the value*/
		void set(Type const& t){ t_=t; }

	private:
		/*!Default copy constructor accessible via clone()*/
		GenericVariable(GenericVariable const&) = default;

		Type t_;//!< Value of the GenericVariable
};
/*}*/

/*{Container*/
/*!Contains a vector of Variable*, can store diffrent GenericVariable*/
class Container{
	public:
		/*!Default constructor*/
		Container() = default;
		/*!Copy constructor*/
		Container(Container const& c){ 
			for(unsigned int i(0);i<c.data_.size();i++){
				data_.push_back((c.data_[i])->clone());
			}
		}
		/*!Destructor*/
		~Container(){
			for(unsigned int i(0);i<data_.size();i++){
				delete data_[i];
				data_[i] = NULL;
			}
		}
		/*{Forbidden*/
		Container(Container&&) = delete;
		Container& operator=(Container) = delete;
		/*}*/

		/*!Add one GenericVariable<Type> of value t and named name*/
		template<typename Type>
			void set(std::string const& name, Type const& t);

		/*!Returns the value GenericVariable<Type>::t_ named name*/
		template<typename Type>
			Type get(std::string const& name) const;
		/*!Set the GenericVariable<Type> named name to t=GenericVariable<Type>::t_ */
		template<typename Type>
			void get(std::string const& name, Type& t) const;
		/*!Sets t to data_[i].get_val()*/
		template<typename Type>
			void get(unsigned int i, Type& t) const;
		/*!Returns data_[i].get_val()*/
		template<typename Type>
			Type get(unsigned int i) const;

		/*!Returns the number of GenericVariable stored*/
		unsigned int size() const { return data_.size(); }

		virtual bool find(std::string const& pattern, unsigned int& i, bool iffail=true) const {
			(void)(iffail);
			i=0;
			while(i<data_.size()){
				if(data_[i]->get_name()==pattern){ return true; } 
				else { i++; }
			}
			return false;
		}

	protected:
		std::vector<Variable*> data_;
};

template<typename Type>
void Container::set(std::string const& name, Type const& t){
	data_.push_back(new GenericVariable<Type>(name,t));
}

template<typename Type>
Type Container::get(std::string const& name) const {
	unsigned int i(0);
	if(find(name,i)){ return static_cast< GenericVariable<Type> *>(data_[i])->get_val(); }
	else {
		std::cerr<<__PRETTY_FUNCTION__<<" : no data with name "<<name<<std::endl;
		return Type();
	}
}

template<typename Type>
void Container::get(std::string const& name, Type& t) const {
	unsigned int i(0);
	if(find(name,i)){ t = static_cast< GenericVariable<Type> *>(data_[i])->get_val(); }
	else { std::cerr<<__PRETTY_FUNCTION__<<" : no data with name '"<<name<<"' thus the value is unchanged : "<< t <<std::endl; }
}

template<typename Type>
void Container::get(unsigned int i, Type &t) const {
	if(i<data_.size()){ t = static_cast< GenericVariable<Type> *>(data_[i])->get_val(); }
	else { std::cerr<<__PRETTY_FUNCTION__<<" : "<<i<<"<"<< data_.size()<<"?" <<std::endl; }
}

template<typename Type>
Type Container::get(unsigned int i) const {
	if(i<data_.size()){
		return static_cast< GenericVariable<Type> *>(data_[i])->get_val();
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : "<<i<<"<"<< data_.size()<<"?" <<std::endl; 
		return Type();
	}
}
/*}*/
#endif
