#include "Miscellaneous.hpp"
#include "Rand.hpp"

inline unsigned int to_fraction_old(double x, unsigned long long& num, unsigned long long& den, double& sign, double const& err = 1e-10){
	sign = my::sign(x);
	num  = std::floor(std::abs(x));
	x = std::abs(x)-num;
	unsigned long long lower_n(0);
	unsigned long long lower_d(1);
	unsigned long long upper_n(1);
	unsigned long long upper_d(1);
	unsigned long long middle_n;
	unsigned long long middle_d;
	unsigned int iter(0);
	while(iter++<1e12){
		middle_n = lower_n + upper_n;
		middle_d = lower_d + upper_d;
		if(middle_d*(x+err)<middle_n){
			upper_n = middle_n;
			upper_d = middle_d;
		} else if(middle_d*(x-err)>middle_n) {
			lower_n = middle_n;
			lower_d = middle_d;
		} else {
			num = num*middle_d+middle_n;
			den = middle_d;
			return iter;
		}
	}
	den = 1;
	std::cerr<<__PRETTY_FUNCTION__<<" : failed to find a fraction for "<<x+num<<std::endl;
	return 0;
}

void test(double a){
	unsigned long long num(1);
	unsigned long long den(1);
	double sign;
	unsigned int iter;
	std::cout<<"-------------------------------------------------------"<<std::endl;
	iter = to_fraction_old(a,num,den,sign);
	std::cout<<"btilly  "<<a<<" "<<sign*num/den<<"="<<num<<"/"<<den<<" in "<<iter<<" iterations"<<std::endl;
	iter = my::to_fraction(a,num,den,sign);                                               
	std::cout<<"richard "<<a<<" "<<sign*num/den<<"="<<num<<"/"<<den<<" in "<<iter<<" iterations"<<std::endl;
}

int main(){
	Rand<double> rnd(-3,3);
	double a(rnd.get());
	test(1./6.);
	test(1./3.);
	test(1./7.);
	test(5./7.);
	test(0.0000001);
	test(2./3.+3);
	test(sqrt(2.));
	test(M_PI);
	test(exp(1));
	test(a);
}
