#ifndef DEF_CELL
#define DEF_CELL

#include <vector>
#include <iostream>
//#include "FonctionTest.hpp"
#include "Epsilon.hpp"

class LastCell {
	public:
		LastCell(double x, double y, double dx, double dy, Fonction* F);
		virtual ~LastCell();
		
		virtual void Affiche() const;

		bool virtual StepSharpen();

		virtual double Integrate();

		double const x_,y_,z_;
	protected:
		double const dx_,dy_,ds_;
		double const I_;
		Fonction* F_;
};

class Cell: public LastCell{
	public:
		Cell(double x, double y, double dx, double dy, Fonction* F);
		virtual ~Cell();	

		virtual void Affiche() const;

		void SharpenGrid(unsigned int max_step);
		bool virtual StepSharpen();

		virtual double Integrate();
	protected:
		std::vector<LastCell*> cells_;
		void Split(unsigned int i);
};


#endif
