#ifndef DEF_CONTAINER
#define DEF_CONTAINER

#include <sstream>
#include <vector>
#include "Read.hpp"

/*{ Data*/
class Data{
	public:
		virtual ~Data(){}

		std::string const& get() const { return name_;}

	protected:
		Data(){}
		Data(std::string name):name_(name){}

		std::string name_;
};
/*}*/

/*{ GenericData*/
template<typename Type>
class GenericData : public Data{
	public:
		GenericData(std::string name, Type t): Data(name), t_(t) {}

		Type const& get() const {return t_;}

	private:
		Type t_;
};
/*}*/

/*{Container*/
class Container{
	public:
		Container(bool copyinstring=false):copyinstring_(copyinstring){ }
		~Container();

		template<typename Type>
			void set(std::string name, Type const& t);

		template<typename Type>
			Type get(std::string name) const;

		std::string value(unsigned int const& i) const{
			if(copyinstring_){ return stringvalue_[i]; }
			else {return ""; }
		}

		std::string const& name(unsigned int const& i) const { return d_[i]->get(); }

		unsigned int size() const { return d_.size(); }

	private:
		std::vector<Data*> d_;
		std::vector<std::string> stringvalue_;
		bool copyinstring_;
};

Container::~Container(){
	for(unsigned int i(0);i<d_.size();i++){
		delete d_[i];
		d_[i] = NULL;
	}
}

template<typename Type>
void Container::set(std::string name, Type const& t){
	d_.push_back(new GenericData<Type>(name,t));
	if(copyinstring_){
		std::ostringstream oss;
		oss<<t;
		stringvalue_.push_back(oss.str());
	}
}

template<typename Type>
Type Container::get(std::string name) const {
	for(unsigned int i(0);i<d_.size();i++){
		if(d_[i]->get()==name){
			return dynamic_cast< GenericData<Type> *>(d_[i])->get();
		}
	}
	std::cerr<<"Container : get(string name) : no data with name "<<name<<std::endl;
	return 0;
}

std::ostream& operator<<(std::ostream& flux, Container const& input){
	for(unsigned int i(0);i<input.size();i++){
		flux<<input.value(i)<<" ";
	}
	return flux;
}
/*}*/

/*{FileParser*/
class FileParser {
	public:
		FileParser(std::string filename):r(filename){};

		template<typename Type>
			void extract(Type& t){ r>>t;}

		template<typename Type>
			void extract(std::string name, Container& c){
				Type t;
				r>>t;
				c.set(name,t);
			}

	private:
		Read r;
};
/*}*/
#endif
