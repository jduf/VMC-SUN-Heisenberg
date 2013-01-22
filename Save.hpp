#ifndef DEF_SAVE
#define DEF_SAVE

#include<iostream>
#include<fstream>
#include<string>


class Save{
	public:
		Save(std::string filename);
		~Save();

		template<typename T>
		Save& operator<<(T const& t);	
		
		static std::string endl;

	private:
		Save();
		Save(Save const& s);
		Save& operator=(Save const&);
		std::ofstream s;
		
};

std::string Save::endl="\n";

Save::Save(std::string filename){
	s.open(filename.c_str(),std::ios::out);
}

Save::~Save(){
	s.close();
}

template<typename T>
Save& Save::operator<<(T const& t){
	s << t;
	return (*this);
}
#endif

