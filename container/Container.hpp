#ifndef DEF_CONTAINER
#define DEF_CONTAINER

#include <vector>
#include "Read.hpp"

/*{ Data*/
class Data{
	public:
		virtual ~Data() = 0;

		virtual Data* clone() const = 0;

		std::string const& get_name() const { return name_;}

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

		GenericData<Type>* clone() const { return new GenericData<Type>(static_cast<GenericData<Type> const&>(*this) );}

		Type const& get_val() const {return t_;}
		void reset(Type const& t) {t_=t;}

	private:
		Type t_;
};
/*}*/

/*{Container*/
class Container{
	public:
		/*!if copyinstrng == true, the variables are saved as string*/
		Container(bool copyinstring=false);
		Container(Container const& c);
		~Container();

		template<typename Type>
			void set(std::string name, Type const& t);

		template<typename Type>
			void reset(std::string name, Type const& t);

		template<typename Type>
			Type get(std::string name) const;

		template<typename Type>
			void get(std::string name, Type& t) const;

		bool check(std::string pattern) const;
		std::string value(unsigned int const& i) const;

		std::string const& name(unsigned int const& i) const { return d_[i]->get_name(); }

		unsigned int size() const { return d_.size(); }

	private:
		bool copyinstring_;
		std::vector<Data*> d_;
		std::vector<std::string> stringvalue_;
};

std::ostream& operator<<(std::ostream& flux, Container const& input);

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
void Container::reset(std::string name, Type const& t){
	unsigned int i(0);
	bool searching(true);
	do{
		if(d_[i]->get_name()==name){
			reinterpret_cast< GenericData<Type> *>(d_[i])->reset(t);
			if(copyinstring_){
				std::ostringstream oss;
				oss<<t;
				stringvalue_[i]=oss.str();
			}
			searching=false;
		}
	} while ( ++i<d_.size() && searching);
	if(searching){
		std::cerr<<"Container : reset(string name, Type const& t) : no data with name "<<name<<std::endl;
	}
}

template<typename Type>
Type Container::get(std::string name) const {
	for(unsigned int i(0);i<d_.size();i++){
		if(d_[i]->get_name()==name){
			return dynamic_cast< GenericData<Type> *>(d_[i])->get_val();
		}
	}
	std::cerr<<"Container : get(string name) : no data with name "<<name<<std::endl;
	return 0;
}

template<typename Type>
void Container::get(std::string name, Type& t) const {
	unsigned int i(0);
	do{
		if(d_[i]->get_name()==name){
			t = dynamic_cast< GenericData<Type> *>(d_[i])->get_val();
			i = d_.size();
		}
		i++;
	}while(i<d_.size());
	if(i==d_.size()){
		std::cerr<<"Container : -"<<name<<" wasn't found thus its value is unchanged : "<< t <<std::endl; 
	}
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
