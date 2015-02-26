#ifndef DEF_CONTAINER
#define DEF_CONTAINER

#include <vector>
#include <iostream>
#include <string>

/*{Variable*/
/*!Base class that will never be instantiate. It provides a way to create
 * GenericVariable with defferent type.*/
class Variable{
	public:
		/*!Constructor*/
		Variable(std::string const& name):name_(name){}
		/*!Destructor, must be virtual*/
		virtual ~Variable(){};

		/*!Method that allows a copy of the derived class*/
		virtual Variable* clone() const = 0;
		/*!Returns the name of the data*/
		std::string const& get_name() const { return name_;}

	protected:
		std::string name_;//!< Name of the (Generic)Variable 
};
/*}*/

/*{GenericVariable*/
/*!Derived class that can contain any one value of any type and of name
 * Variable::name_*/
template<typename Type>
class GenericVariable : public Variable {
	public:
		/*!Constructor*/
		GenericVariable(std::string const& name, Type t):Variable(name),t_(t){}
		/*!Destructor*/
		~GenericVariable(){};

		/*!Returns the value of the data*/
		Type const& get_val() const {return t_;}
		/*!Set the value*/
		void set(Type const& t) {t_=t;}

	private:
		Type t_;//!< Value of the GenericVariable

		/*!Method that implements*/
		GenericVariable<Type>* clone() const { return new GenericVariable<Type>(*this);}
};
/*}*/

/*{Container*/
/*!Contains a vector of Variable*, can store diffrent GenericVariable*/
class Container{
	public:
		Container(){};
		Container(Container const& c){ 
			for(unsigned int i(0);i<c.data_.size();i++){
				data_.push_back((c.data_[i])->clone());
			}
		}
		~Container(){
			for(unsigned int i(0);i<data_.size();i++){
				delete data_[i];
				data_[i] = NULL;
			}
		}

		/*!Add one GenericVariable<Type> of value t and named name*/
		template<typename Type>
			void set(std::string const& name, Type const& t);

		/*!Returns the value GenericVariable<Type>::t_ named name*/
		template<typename Type>
			Type get(std::string const& name);
		/*!Set the GenericVariable<Type> named name to t=GenericVariable<Type>::t_ */
		template<typename Type>
			void get(std::string const& name, Type& t);
		/*!Sets t to data_[i].get_val()*/
		template<typename Type>
			void get(unsigned int i, Type& t);
		/*!Returns data_[i].get_val()*/
		template<typename Type>
			Type get(unsigned int i);

		/*!Returns the number of GenericVariable stored*/
		unsigned int size() const { return data_.size(); }

		virtual bool find(std::string const& pattern, unsigned int& i, bool iffail=true){
			(void)(iffail);
			while(i<data_.size()){
				if(data_[i]->get_name()==pattern){ return true; } 
				else { i++; }
			}
			return false;
		}

	protected:
		std::vector<Variable*> data_;

		template<typename Type>
			Type const& do_cast(unsigned int const& i) const {
				return reinterpret_cast< GenericVariable<Type> *>(data_[i])->get_val();
				//dynamic_cast< GenericVariable<Type> *>(data_[i])->get_val();
			}
};

template<typename Type>
void Container::set(std::string const& name, Type const& t){
	data_.push_back(new GenericVariable<Type>(name,t));
}

template<typename Type>
Type Container::get(std::string const& name){
	unsigned int i(0);
	if(find(name,i)){ return do_cast<Type>(i); }
	else {
		std::cerr<<"Container : get(string name) : no data with name "<<name<<std::endl;
		return Type();
	}
}

template<typename Type>
void Container::get(std::string const& name, Type& t){
	unsigned int i(0);
	if(find(name,i)){ t = do_cast<Type>(i); }
	else {
		std::cerr<<"Container : get(string name) : no data with name '"<<name<<"' thus the value is unchanged : "<< t <<std::endl; 
	}
}

template<typename Type>
void Container::get(unsigned int i, Type &t){
	if(i<data_.size()){ t = do_cast<Type>(i); }
	else { std::cerr<<"Container : "<<i<<"<"<< data_.size()<<"?" <<std::endl; }
}

template<typename Type>
Type Container::get(unsigned int i){
	if(i<data_.size()){
		//return reinterpret_cast< GenericVariable<Type> *>(data_[i])->get_val();
		return do_cast<Type>(i);
	} else {
		std::cerr<<"Container : "<<i<<"<"<< data_.size()<<"?" <<std::endl; 
		return Type();
	}
}
/*}*/

///*{FileParser*/
//class FileParser{
	//public:
		//FileParser(std::string const& filename):r(filename,false){};
//
		//template<typename Type>
			//void extract(Type& t){ r>>t;}
//
		//template<typename Type>
			//void transfer_to_container(std::string const& name, Container& c){
				//Type t;
				//r>>t;
				//c.set(name,t);
			//}
//
	//private:
		//IOFiles r;
//};
///*}*/
#endif
