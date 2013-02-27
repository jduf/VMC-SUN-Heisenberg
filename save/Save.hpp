#ifndef DEF_SAVE
#define DEF_SAVE

#include<iostream>
#include<fstream>
#include<string>


class Save{
	public:
		/*!
		 * \param filename name of the file where the datas will be saved
		 * \param overwrite if true will overwrite any exististing file*/
		Save(std::string filename, bool overwrite=true);
		/*!close the file*/
		~Save();

		/*!operator that writes the t value in the filename file */
		template<typename T>
		Save& operator<<(T const& t);	
		
		/*!provide a way to end a line (may be improved)*/
		static std::string endl;

	private:
		/*!forbids default constructor*/
		Save();
		/*!forbids copy constructor*/
		Save(Save const& s);
		/*!forbids assertion operator*/
		Save& operator=(Save const&);

		std::ofstream s; //!< filename
};

std::string Save::endl="\n";

Save::Save(std::string filename,bool overwrite){
	if(overwrite){ s.open(filename.c_str(),std::ios::out | std::ios::trunc); }
	else {s.open(filename.c_str(),std::ios::out | std::ios::app);}
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

