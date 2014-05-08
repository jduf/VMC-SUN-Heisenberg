#include "Fit.hpp"

#include "Read.hpp"
#include "Gnuplot.hpp"

double f(double x, double *p);

int main(){
	Read r("data.dat");
	unsigned int s(2);
	Matrix<double> xy(49,2);
	Vector<double> x(49-s);
	Vector<double> y(49-s);
	r>>xy;
	for(unsigned int i(s);i<xy.row();i++){
		x(i-s) = xy(i,0);
		y(i-s) = xy(i,1);
	}

	Vector<double> p(4,1);
	Fit test(x,y,*f,p);
	Gnuplot gp("./","image");
	Vector<double> yfit(test.fx());
	Write w("fit.dat");
	for(unsigned int i(0);i<x.size();i++){
		w<<x(i)<<" "<<yfit(i)<<Write::endl;
	}
	gp="'fit.dat' t sprintf('$\\eta=%f$',"+tostring(p(2))+"), 'data.dat'";
	gp.save_file();
	gp.create_image(true);

}

double f(double x, double *p){
	return p[0]/(x*x)+p[1]*cos(2*M_PI*x/3.0)/pow(x,p[2])+p[3];
}

