#include "Save.hpp"
#include "Vecteur.hpp"

int main(){
	Save save_stream("filename",false);
	double x(10.1);
	save_stream<<"x="<<x<<Save::endl;

	save_stream<<"one can write any type in that file ";
	save_stream<<"oups, I forgot to end the line..."<<Save::endl;
	
	Vecteur<double> v(3,2.1);
	save_stream<<v<<Save::endl;
}

