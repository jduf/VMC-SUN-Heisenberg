/* @file test.cpp */

#include "PSpline.hpp"
#include "Rand.hpp"
#include "Gnuplot.hpp"

void test_given_function();
void test_211();
void test_212();

int main(){
	test_given_function();
	test_211();
	test_212();
}

void test_given_function(){
	PSpline s(2);

	Rand<double> rnd(-1,1);
	auto f = [](Vector<double> const& x){
		return x(0)*x(1)+cos(x(1))*exp(x(0)*x(0));
	};
	IOFiles data("surf.dat",true);
	Vector<double> c_tmp(2);
	double y_tmp;

	for(unsigned int i(0);i<100;i++){
		c_tmp(0) = rnd.get();
		c_tmp(1) = rnd.get();
		y_tmp = f(c_tmp);
		s.add_data(c_tmp,y_tmp);
		data<<c_tmp<<" "<<y_tmp<<IOFiles::endl;
	}

	s.compute_weights();
	IOFiles out("spline.dat",true);
	double min(-1.0);
	double dx(0.04);
	Vector<double> tmp(2);
	for(unsigned int i(0);i<50;i++){
		tmp(0) = min+i*dx;
		for(unsigned int j(0);j<50;j++){
			tmp(1) = min+j*dx;
			out<<tmp<<" "<<s.extrapolate(tmp)<<IOFiles::endl;
		}
		out<<tmp<<IOFiles::endl;
	}

	Gnuplot plot("./","plot");
	plot.range("x",-1,1);
	plot.range("y",-1,1);
	plot+="f(x,y) = x*y+cos(y)*exp(x*x)";
	plot+="splot 'surf.dat' u 1:2:3 notitle,\\";
	plot+="      'spline.dat' u 1:2:3 notitle,\\";
	plot+="      f(x,y)";
	plot.save_file();
}

void test_211(){
	IOFiles in("211.dat",false);
	Matrix<double> data(1055,8);//max = 1055
	in>>data;
	PSpline s(5);
	Vector<double> tmp(2);
	for(unsigned int i(0);i<data.row();i++){
		tmp(0) = data(i,1);
		tmp(1) = data(i,2);
		s.add_data(tmp,data(i,3));
	}
	s.compute_weights();

	IOFiles out("spline.dat",true);
	double min(-2.0);
	double dx(0.05);
	for(unsigned int i(0);i<80;i++){
		tmp(0) = min+i*dx;
		for(unsigned int j(0);j<80;j++){
			tmp(1) = min+j*dx;
			out<<tmp<<" "<<s.extrapolate(tmp)<<IOFiles::endl;
		}
		out<<tmp<<IOFiles::endl;
	}

	Gnuplot plot("./","plot");
	plot.range("x",-2,2);
	plot.range("y",-2,2);
	plot+="splot '211.dat' u 2:3:4 notitle,\\";
	plot+="      'spline.dat' u 1:2:3 notitle";
	plot.save_file();
}

void test_212(){
	IOFiles in("212.dat",false);
	Matrix<double> data(857,8);//max = 857
	in>>data;
	PSpline s(1);
	Vector<double> tmp(2);
	for(unsigned int i(0);i<data.row();i++){
		tmp(0) = data(i,1);
		tmp(1) = data(i,2);
		s.add_data(tmp,data(i,3));
	}
	s.compute_weights();

	IOFiles out("spline.dat",true);
	double min(-2.0);
	double dx(0.05);
	for(unsigned int i(0);i<80;i++){
		tmp(0) = min+i*dx;
		for(unsigned int j(0);j<80;j++){
			tmp(1) = min+j*dx;
			out<<tmp<<" "<<s.extrapolate(tmp)<<IOFiles::endl;
		}
		out<<tmp<<IOFiles::endl;
	}

	Gnuplot plot("./","plot");
	plot.range("x",-1,1);
	plot.range("y",-1,1);
	plot+="splot '212.dat' u 2:3:4 notitle,\\";
	plot+="      'spline.dat' u 1:2:3 notitle";
	plot.save_file();
}
