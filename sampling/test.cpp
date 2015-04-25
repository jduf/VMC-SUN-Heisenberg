#include "Sampling.hpp"
#include "Gnuplot.hpp"
#include "Rand.hpp"
#include<stdlib.h>
#include<time.h>

void check_troyer();
void check_flip_coin();
void read_output_files();
void merge_two_sim();

int main(){
	//check_troyer();
	//check_flip_coin();
	//read_output_files();
	merge_two_sim();
}

void check_troyer(){
	std::cout<<"tiré de Mathias Troyer"<<std::endl<<std::endl;
	std::cout<<"+ il faut que le graphe log.pdf montre une convergence (pente null)";
	std::cout<<"pour M=100000, mais pas pour M=10000 si p=24"<<std::endl;
	std::cout<<"+ les barre d'erreur de plot.pdf doivent contenir P(x)"<<std::endl;
	std::cout<<"+ avec M=1e4 devrait reproduire fig 4 de son article"<<std::endl<<std::endl;


	unsigned int const M(1e4);
	unsigned int const N(50);
	unsigned int iter(0);
	unsigned int n(0);
	unsigned int r;
	double tol(5e-4);
	Rand<double> rnd(0,N);
	do{
		r = rnd.get();
		if(r < n){n--;}
		if(r > n){n++;}
		if(r != n){ iter++; }
	} while (iter<M/5);

	unsigned int B(50);
	unsigned int b(5);

	Vector<double> H_wrong(N+1,0);
	DataSet<double> H_right;
	H_right.set(N+1,B,b,true);
	iter=0;

	do{
		r = rnd.get();
		if(r < n){n--;}
		if(r > n){n++;}
		if(r != n){
			iter++;
			H_wrong(n) +=1.0;

			H_right[n].set_x(1.0);
			H_right.add_sample();
			H_right[n].set_x(0.0);
			if(!(iter % (M/10))){ H_right.compute_convergence(tol); }
		}
	} while (iter<M);
	H_right.complete_analysis(tol);

	IOFiles hist_right("hist_right.dat",true);
	for(unsigned int i(0);i<N+1;i++){ 
		hist_right<<i<<" "<<H_right[i].get_x()<<" "<<H_right[i].get_dx()<<" "<<H_right[i].get_conv()<<IOFiles::endl;
	}
	IOFiles hist_wrong("hist_wrong.dat",true);
	for(unsigned int i(0);i<N+1;i++){ 
		hist_wrong<<i<<" "<<H_wrong(i)/M<<" "<<sqrt((H_wrong(i)/M)*(1-H_wrong(i)/M)/(iter-1))<<IOFiles::endl;
	}

	Gnuplot gp("./","plot");
	gp.range("x",10,40);
	gp.range("y",0,0.12);
	gp+="N = "+my::tostring(N);
	gp+="fac(x) = (int(x)==0) ? 1.0 : int(x) * fac(int(x)-1.0)";
	gp+="P(x) = fac(N)/(2**N*fac(x)*fac(N-x))";
	gp+="set sample 31";
	gp+="plot P(x) w l t '$P(x)$',\\";
	gp+="     'hist_right.dat' u 1:($4==1?$2:1/0):3 w errorbars t 'correct converged',\\";
	gp+="     'hist_right.dat' u 1:($4==0?$2:1/0):3 w errorbars t 'correct not converged',\\";
	gp+="     'hist_wrong.dat' w errorbars  t 'wrong'";
	gp.save_file();
	gp.create_image(true);
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
	IOFiles w("coin.jdbin",true);
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
			if(iter>1e4){ H.compute_convergence(tol); }
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
	IOFiles flip("coin.jdbin",false);
	if(flip.is_open()){
		unsigned int N;
		flip>>N;
		Data<double> H;
		for(unsigned int i(0);i<N;i++){
			flip>>H;
			std::cout<<H<<std::endl;
		}
	}
	IOFiles troyer("coin.jdbin",false);
	if(troyer.is_open()){
	}
}

void merge_two_sim(){
	Rand<double> rnd1(0,1);
	Rand<double> rnd2(0,0.5);
	unsigned int iter(0);

	unsigned int B(5);
	unsigned int b(2);

	double tol(5e-5);
	double conv(0.0);

	Data<double> H1;
	H1.set(B,b,conv);
	Data<double> H2;
	H2.set(B,b,conv);
	do{ 
		iter++;
		H1.set_x(rnd1.get());
		H1.add_sample(); 
	} while (iter<60);
	iter = 0;
	do{ 
		iter++;
		H2.set_x(rnd2.get());
		H2.add_sample(); 
	} while (iter<600);

	H1.complete_analysis(tol);
	H2.complete_analysis(tol);

	std::cout<<H2<<std::endl;
	std::cout<<H1<<std::endl;

	H2.merge(H1);
	H2.complete_analysis(tol);

	std::cout<<H2<<std::endl;
}
