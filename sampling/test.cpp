#include "Sampling.hpp"
#include "Gnuplot.hpp"
#include "Rand.hpp"
#include <stdlib.h>
#include <time.h>

void check_troyer();
void check_flip_coin();
void read_output_files();
void merge_two_sim();
void merge_two_sim_rnd();

int main(){
	check_troyer();
	//check_flip_coin();
	//read_output_files();
	//merge_two_sim();
	//merge_two_sim_rnd();
}

void check_troyer(){
	std::cout<<"tiré de Mathias Troyer"<<std::endl<<std::endl;
	std::cout<<"+ les barre d'erreur de plot.pdf doivent contenir P(x)"<<std::endl;
	std::cout<<"+ avec M=1e5 devrait reproduire fig 6 de son article"<<std::endl<<std::endl;

	unsigned int const M(1e5);
	unsigned int const N(50);
	unsigned int iter(0);
	unsigned int n(0);
	unsigned int r;
	double tol(1e-3);
	Rand<unsigned int> rnd(1,N);

	unsigned int B(50);
	unsigned int b(5);

	DataSet<double> H_correct;
	H_correct.set(N+1,B,b,true);
	do{
		r = rnd.get();
		if(r <= n){n--;}
		else{n++;}

		if(iter>M/5){
			H_correct[n].set_x(1.0);
			H_correct.add_sample();
			H_correct[n].set_x(0.0);
		}
	} while (iter++<M);
	H_correct.complete_analysis(tol);

	IOFiles hist_correct("hist_correct.dat",true,false);
	for(unsigned int i(0);i<N+1;i++){ 
		hist_correct<<i<<" "<<H_correct[i].get_x()<<" "<<H_correct[i].get_dx()<<" "<<H_correct[i].get_conv()<<IOFiles::endl;
	}

	Gnuplot gp("./","plot");
	unsigned int min(10);
	unsigned int max(40);
	gp.range("x",min,max);
	gp.range("y","0","");
	gp+="N = "+my::tostring(N);
	gp+="fac(x) = (int(x)==0) ? 1.0 : int(x) * fac(int(x)-1.0)";
	gp+="P(x) = fac(N)/(2**N*fac(x)*fac(N-x))";
	gp+="set sample "+my::tostring(max-min+1);
	gp+="plot P(x) w l t '$P(x)$',\\";
	gp+="     'hist_correct.dat' u 1:($4==1?$2:1/0):3 w errorbars t 'correct converged',\\";
	gp+="     'hist_correct.dat' u 1:($4==0?$2:1/0):3 w errorbars t 'correct not converged'";
	gp.save_file();
	gp.create_image(true,false);
}

void check_flip_coin(){
	std::cout<<"vérifier que ce soit bien correct car bcp de fois ou diff=1e-14... peut-être pas une erreur numérique "<<std::endl<<std::endl;
	std::cout<<"Il est normal si des erreurs s'affichent. "<<std::endl<<std::endl;
	std::cout<<"+ si diff<1e-10, c'est simplement que la tolérance est trop ";
	std::cout<<"basse et que par conséquent, le premier test de convergence";
	std::cout<<"est vrai. Alors forcément xtot!=H.get_mean()."<<std::endl;
	std::cout<<"+ si 1e-10<diff<1e-16, c'est une approximation numérique"<<std::endl<<std::endl;

	Rand<double> rnd(0.0,1.0);
	unsigned int N(10);
	unsigned int iter(0);

	unsigned int B(50);
	unsigned int b(5);

	double tol(5e-5);
	double sup(0.0);
	double diff(0.0);
	double conv(0.0);

	Data<double> H;
	H.set(B,b,conv);
	double xtot;
	double x;
	IOFiles w("coin.jdbin",true,false);
	w<<N;
	for(unsigned int i(0);i<N;i++){
		iter=0;
		xtot=0.0;
		do{ 
			iter++;
			x=rnd.get();
			H.set_x(x);
			H.add_sample(); 
			xtot += x;
			if(iter>1e4){ H.complete_analysis(tol); }
		} while (!H.get_conv());
		H.complete_analysis(tol);

		xtot/=iter;
		if(H.get_x()>0.5){ sup+=1.0;} 
		/*{if the following test is wrong it is because the xtot and
		 * H_get_mean() are not computed for the same number of samples. this
		 * is normal because CorrelatedSamples recompute x_ and dx_ only when one bin is
		 * added in the smallest CorrelatedSamples. So if one bin is partially "filled",
		 * the corresponding samples are not taken into account. They will be
		 * when the bin is filled because this will allow a new computation of
		 * x_ and dx_. To recover the equality at all time, the computation of
		 * x_ should be :
		 *
		 x_=0.0;
		 for(unsigned int i(0);i<Ml_(0);i++){ x_ += bin_[0](i); }
		 if(iter_!=0){ 
		 x_ *= DPL_;
		 x_ += bin_[0](Ml_(0));
		 x_ /= DPL_*Ml_(0)+iter_;
		 } else {
		 x_/=Ml_(0);
		 }

		 * where iter_ is the number of samples stored in the partially filled
		 * bin. But if one uses this x_, there won't be any corresponding dx_
		 *}*/
		if(std::abs(H.get_x()-xtot)<1e-14){ diff+=1.0; } 
		else { std::cout<<"diff "<<std::abs(H.get_x()-xtot)<<" "<<H.get_x()<<" "<<xtot<<std::endl; }
		if(H.get_conv()){ conv+=1.0;}
		w<<H;
		H.set(B,b,true);
	}
	std::cout<<std::endl<<"=> probabilité x in [0:0.5] = "<<sup/N*100<<"%"<<std::endl;
	std::cout<<"=> probabilité xtot==H.get_mean() = "<<diff/N*100<<"%"<<std::endl;
	std::cout<<"=> pourcentage de simulation ayant convergée = "<<conv/N*100<<"%"<<std::endl;
}

void read_output_files(){
	IOFiles flip("coin.jdbin",false,false);
	if(flip.is_open()){
		unsigned int N;
		flip>>N;
		Data<double> H;
		for(unsigned int i(0);i<N;i++){
			flip>>H;
			std::cout<<H<<std::endl;
		}
	}
	IOFiles troyer("coin.jdbin",false,false);
	if(troyer.is_open()){
	}
}

void merge_two_sim(){
	unsigned int iter(0);
	unsigned int N2(345);
	unsigned int N1(123);

	unsigned int B(5);
	unsigned int b(2);

	double tol(5e-5);
	double conv(0.0);

	Data<double> H1;
	H1.set(B,b,conv);
	Data<double> H2;
	H2.set(B,b,conv);
	Data<double> H3;
	H3.set(B,b,conv);
	do{ 
		iter++;
		H1.set_x(iter);
		H1.add_sample(); 

		H3.set_x(iter);
		H3.add_sample(); 
	} while (iter<N1);
	iter = 0;
	do{ 
		iter++;
		H2.set_x(iter);
		H2.add_sample(); 

		H3.set_x(iter);
		H3.add_sample(); 
	} while (iter<N2);

	std::cout<<H2<<std::endl;
	std::cout<<H1<<std::endl;

	H2.merge(H1);
	H2.complete_analysis(tol);
	H2.delete_binning();

	std::cout<<H2<<std::endl;

	H3.complete_analysis(tol);
	H3.delete_binning();
	std::cout<<H3<<std::endl;
}

void merge_two_sim_rnd(){
	unsigned int iter(0);
	Rand<double> rnd(0.,1.);
	double tmp;
	unsigned int N2(3456789);
	unsigned int N1(1234567);

	unsigned int B(5);
	unsigned int b(2);

	double tol(5e-5);
	double conv(0.0);

	Data<double> H1;
	H1.set(B,b,conv);
	Data<double> H2;
	H2.set(B,b,conv);
	Data<double> H3;
	H3.set(B,b,conv);
	do{ 
		tmp = rnd.get();
		iter++;
		H1.set_x(tmp);
		H1.add_sample(); 

		H3.set_x(tmp);
		H3.add_sample(); 
	} while (iter<N1);
	iter = 0;
	do{ 
		tmp = rnd.get();
		iter++;
		H2.set_x(tmp);
		H2.add_sample(); 

		H3.set_x(tmp);
		H3.add_sample(); 
	} while (iter<N2);

	std::cout<<H2<<std::endl;
	std::cout<<H1<<std::endl;

	H2.merge(H1);
	H2.complete_analysis(tol);
	H2.delete_binning();

	std::cout<<H2<<std::endl;

	H3.complete_analysis(tol);
	H3.delete_binning();
	std::cout<<H3<<std::endl;
}
