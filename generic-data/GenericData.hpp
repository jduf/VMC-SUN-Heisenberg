#ifndef DEF_DATA
#define DEF_DATA

#include<string>
#include<vector>
#include "Read.hpp"

/*{ Data*/
class Data{
	public:
		Data(std::string name);
		virtual ~Data();
		virtual void print() = 0;
		virtual bool compare(std::string name);

	protected:
		std::string name_;
};

Data::Data(std::string name):
	name_(name)
{}

Data::~Data(){}
bool Data::compare(std::string name){
	if(name == name_) { return true; }
	else { return false; }
}
/*}*/

/*{ GenericData*/
template<typename Type>
class GenericData : public Data{
	public:
		GenericData(std::string name, Type t);
		void print();
		
		Type get();

	private:
		Type t_;
};

template<typename Type>
GenericData<Type>::GenericData(std::string name, Type t):
	Data(name),
	t_(t)
{}

template<typename Type>
Type GenericData<Type>::get(){
	return t_;
}

template<typename Type>
void GenericData<Type>::print(){
	std::cout<<name_<<" "<<t_<<std::endl;
}
/*}*/

/*{Input*/
class Input{
	public:
		Input();
		~Input();

		template<typename Type>
		void read(std::string name, Read& r);

		template<typename Type>
		void set(std::string name, Type const& t);

		template<typename Type>
		Type get(std::string name);
		
	private:
		std::vector<Data*> d_;
};

Input::Input(){ }

Input::~Input(){
	for(unsigned int i(0);i<d_.size();i++){
		delete d_[i];
		d_[i] = NULL;
	}
}

template<typename Type>
void Input::read(std::string name, Read& r){
	Type t;
	r>>t;
	d_.push_back(new GenericData<Type>(name,t));
}

template<typename Type>
void Input::set(std::string name, Type const& t){
	d_.push_back(new GenericData<Type>(name,t));
}

template<typename Type>
Type Input::get(std::string name){
	for(unsigned int i(0);i<d_.size();i++){
		if(d_[i]->compare(name)){
			return dynamic_cast< GenericData<Type> *>(d_[i])->get();
		}
	}
	std::cerr<<"Input : get(string name) : no data with name "<<name<<std::endl;
	return 0;
}
/*}*/
#endif
