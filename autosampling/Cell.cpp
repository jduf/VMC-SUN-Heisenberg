#include "Cell.hpp"

LastCell::LastCell(double x, double y, double dx, double dy, Fonction* F):x_(x),y_(y),z_((*F)(x,y)),dx_(dx),dy_(dy),ds_(dx*dy),I_(ds_*z_),F_(F){}

//je ne comprends pas pourquoi je dois définir ce destructeur...
LastCell::~LastCell(){}

void LastCell::Affiche() const {
	//std::cout<<"("<<x_<<","<<y_<<","<<z_<<") -> ("<<dx_<<","<<dy_<<","<<I_<<")"<<std::endl;
	std::cout<<x_<<" "<<y_<<" "<<z_<<std::endl;
}
 
void Cell::Affiche() const {
	for(unsigned int i(0); i<cells_.size(); i++){
		cells_[i]->Affiche();
	}
}

Cell::Cell(double x, double y, double dx, double dy, Fonction* F):LastCell(x,y,dx,dy,F){
	cells_.push_back(new LastCell(x-dx/4.0,y-dy/4.0,dx/2.0,dy/2.0,F));
	cells_.push_back(new LastCell(x-dx/4.0,y+dy/4.0,dx/2.0,dy/2.0,F));
	cells_.push_back(new LastCell(x+dx/4.0,y-dy/4.0,dx/2.0,dy/2.0,F));
	cells_.push_back(new LastCell(x+dx/4.0,y+dy/4.0,dx/2.0,dy/2.0,F));
}

Cell::~Cell(){
	for(unsigned int i(0); i<cells_.size();i++){
		delete cells_[i];
	}
	cells_.clear();
}

void Cell::Split(unsigned int i){
	double x(cells_[i]->x_);
	double y(cells_[i]->y_);
	delete cells_[i];
	cells_[i] = new Cell(x,y,dx_/2.0,dy_/2.0,F_);
}

double Cell::Integrate(){
	double I(0.0);
	for(unsigned int i(0); i<cells_.size(); i++){
		I += cells_[i]->Integrate(); 
	}
	return I;
}

double LastCell::Integrate(){ return I_;}

bool Cell::StepSharpen(){
	for(unsigned int i(0);i<cells_.size();i++){
		if(cells_[i]->StepSharpen()) {
			Split(i);
		}
	}
	return false;
}

bool LastCell::StepSharpen(){
	double I(0.0);
	I +=(*F_)( x_-dx_/2.0, y_-dy_/2.0 ) * ds_ / 4;
	I +=(*F_)( x_-dx_/2.0, y_+dy_/2.0 ) * ds_ / 4;
	I +=(*F_)( x_+dx_/2.0, y_-dy_/2.0 ) * ds_ / 4;
	I +=(*F_)( x_+dx_/2.0, y_+dy_/2.0 ) * ds_ / 4;

	if( fabs((I -I_)/I_) > 0.1){
		//std::cout<<"vrai "<<x_<<" "<<y_<< " "<<ds_<<" " <<I_<< " " <<I<<" "<<fabs((I -I_)/I_)<<std::endl;
		return true;
	} else {
		//std::cout<<"faux"<<std::endl;
		return false;
	}
}

void Cell::SharpenGrid(unsigned int max_step){
	for (unsigned int i(0);i<max_step;i++){
		StepSharpen();
	}
}
