#ifndef DEF_CONTAINER
#define DEF_CONTAINER

#include <vector>
#include "IOFiles.hpp"

/*{ Data*/
/*!Base class that will never be instantiate. It provides a way to create
 * GenericData with defferent type. */
class Data{
	public:
		/*!Constructor*/
		Data(std::string name):name_(name){}
		/*!Destructor, must be virtual*/
		virtual ~Data(){};

		/*!Method that allows a copy of the derived class*/
		virtual Data* clone() const = 0;
		/*!Returns the name of the data*/
		std::string const& get_name() const { return name_;}

	protected:
		std::string name_;//!< Name of the (Generic)Data 
};
/*}*/

/*{GenericData*/
/*!Derived class that can contain any one value of any type and of name
 * Data::name_*/
template<typename Type>
class GenericData : public Data{
	public:
		/*!Constructor*/
		GenericData(std::string name, Type t): Data(name), t_(t){}
		/*!Destructor*/
		~GenericData(){};

		/*!Method that implements*/
		GenericData<Type>* clone() const { return new GenericData<Type>(*this);}

		/*!Returns the value of the data*/
		Type const& get_val() const {return t_;}
		/*!Set the value*/
		void set(Type const& t) {t_=t;}

	private:
		Type t_;//!< Value of the GenericData
};
/*}*/

/*{Container*/
/*!Contains a vector of Data*, can store diffrent GenericData*/
class Container{
	public:
		/*!if copyinstrng == true, the variables are saved as string*/
		//Container(bool copyinstring=false);
		Container();
		Container(Container const& c);
		~Container();

		/*!Add one GenericData<Type> of value t and named name*/
		template<typename Type>
			void set(std::string name, Type const& t);

		//template<typename Type>
			//void reset(std::string name, Type const& t);

		/*!Returns the value GenericData<Type>::t_ named name*/
		template<typename Type>
			Type get(std::string name) const;
		/*!Set the GenericData<Type> named name to t=GenericData<Type>::t_ */
		template<typename Type>
			void get(std::string name, Type& t) const;
		/*!Check if containers has a value named name*/
		bool check(std::string name) const;
		/*!Returns the name of the ith value*/
		std::string const& name(unsigned int const& i) const { return data_[i]->get_name(); }
		/*!Returns the number of GenericData stored*/
		unsigned int size() const { return data_.size(); }

	private:
		//bool copyinstring_;
		std::vector<Data*> data_;
		//std::vector<std::string> stringvalue_;
};

//std::ostream& operator<<(std::ostream& flux, Container const& input);

template<typename Type>
void Container::set(std::string name, Type const& t){
	data_.push_back(new GenericData<Type>(name,t));
	//if(copyinstring_){
		//std::ostringstream oss;
		//oss<<t;
		//stringvalue_.push_back(oss.str());
	//}
}

//template<typename Type>
//void Container::reset(std::string name, Type const& t){
	//unsigned int i(0);
	//do{
		//if(data_[i]->get_name()==name){
			//reinterpret_cast< GenericData<Type> *>(data_[i])->set(t);
			//if(copyinstring_){
				//std::ostringstream oss;
				//oss<<t;
				//stringvalue_[i]=oss.str();
			//}
		//}
	//} while ( ++i<data_.size());
	//if(i==data_.size()){
		//std::cerr<<"Container : reset(string name, Type const& t) : no data with name "<<name<<std::endl;
	//}
//}

template<typename Type>
Type Container::get(std::string name) const {
	for(unsigned int i(0);i<data_.size();i++){
		if(data_[i]->get_name()==name){
			return reinterpret_cast< GenericData<Type> *>(data_[i])->get_val();
		}
	}
	std::cerr<<"Container : get(string name) : no data with name "<<name<<std::endl;
	return Type();
}

template<typename Type>
void Container::get(std::string name, Type& t) const {
	unsigned int i(0);
	do{
		if(data_[i]->get_name()==name){
			t = dynamic_cast< GenericData<Type> *>(data_[i])->get_val();
			i = data_.size();
		}
	}while(++i<data_.size());
	if(i==data_.size()){
		std::cerr<<"Container : '"<<name<<"' wasn't found thus the value is unchanged : "<< t <<std::endl; 
	}
}
/*}*/

/*{FileParser*/
class FileParser {
	public:
		FileParser(std::string filename):r(filename,false){};

		template<typename Type>
			void extract(Type& t){ r>>t;}

		template<typename Type>
			void extract(std::string name, Container& c){
				Type t;
				r>>t;
				c.set(name,t);
			}

	private:
		IOFiles r;
};
/*}*/
#endif
