#ifndef DEF_READ
#define DEF_READ

#include"Matrice.hpp"
#include<iostream>
#include<fstream>
#include<string>


class Read{
	public:
		/*!
		 * \param filename name of the file where the datas will be readd
		 * \param overwrite if true will overwrite any exististing file*/
		Read(std::string filename, Matrice<double>& m);

	private:
		/*!forbids default constructor*/
		Read();
		/*!forbids copy constructor*/
		Read(Read const& R);
		/*!forbids assertion operator*/
		Read& operator=(Read const&);
};


#endif

