#include "Fit.hpp"
#include "Rand.hpp"
#include "IOFiles.hpp"
#include "Gnuplot.hpp"

class Myclass{
	public:
		Myclass():N_(3),m_(1),n_(60){};

		void fitquad(){
			const unsigned int m(20);
			Vector<double> y(m);
			Vector<double> yfit(m);
			Vector<double> x(m);
			Rand<double> rnd(0,5);
			for(unsigned int i(0);i<m;i++){
				x(i) = i;
				y(i) = 0.5*i*i+3*i+4+rnd.get();
			}

			Vector<double> p(3,1);
			Fit F(x,y,p,[this](double x, const double* p){return f(x,p);});

			IOFiles data("data.dat",true);
			for(unsigned int i(0);i<m;i++){
				 data<<i<<" "<<y(i)<<" "<<f(x(i),p.ptr())<<IOFiles::endl;
			}

			Gnuplot gp("./","plot-1");
			gp+="plot 'data.dat' u 1:2 notitle,\\";
			gp+="     'data.dat' u 1:3 notitle";
			gp.save_file();
			gp.create_image(true);
			std::cout<<p<<std::endl;
		}

		void fitcorr(){
			IOFiles r("chain-fermi-N3-m1-n60-M-20-20-20-P-long-range-corr.dat",false);

			unsigned int s(3);
			Matrix<double> xy(60,6);
			Vector<double> x(60-2*s);
			Vector<double> y(60-2*s);
			r>>xy;
			for(unsigned int i(s);i<xy.row()-s;i++){
				x(i-s) = xy(i,0);
				y(i-s) = xy(i,1);
			}

			Vector<double> p(3,1);

			auto func = [this](double x, const double* p){ 
				//return p[0]*cos(2*M_PI*x*m_/N_)*(pow(x,-p[1])+pow(n_-x,-p[1]))+p[2]*(pow(x,-p[3])+pow(n_-x,-p[3]));
				return p[0]*cos(2*M_PI*x*this->m_/this->N_)*(pow(x,-p[1])+pow(this->n_-x,-p[1]))+p[2]*(pow(x,-2)+pow(this->n_-x,-2));
			};
			Fit F(x,y,p,func);

			IOFiles data("data.dat",true);
			for(unsigned int i(0);i<x.size();i++){
				data<<x(i)<<" "<<y(i)<<" "<<func(x(i),p.ptr())<<IOFiles::endl;
			}

			Gnuplot gp("./","plot-2");
			gp+="plot 'data.dat' u 1:2 notitle,\\";
			gp+="     'data.dat' u 1:3 notitle";
			gp.save_file();
			gp.create_image(true);
			std::cout<<p<<std::endl;
		}
	private:
		double f(double x, const double* p);
		unsigned int N_;
		unsigned int m_;
		unsigned int n_;
};

int main() {
	Myclass mc;
	//mc.fitquad();
	mc.fitcorr();
}

double Myclass::f(double x, const double* p){
	return p[0]*x+p[1]*x*x+p[2];
}
