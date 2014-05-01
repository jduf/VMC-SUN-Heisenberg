
#include "Binning.hpp"
#include "Rand.hpp"
#include <iostream>
#include<stdlib.h>
#include<time.h>

void check_troyer();
void check_flip_coin();

int main(){
	//check_troyer();
	check_flip_coin();
}

void check_troyer(){
	std::cout<<"tiré tiré de Mathias Troyer"<<std::endl<<std::endl;
	std::cout<<"+ il faut que le graphe log.pdf montre une convergence (pente null)";
	std::cout<<"pour M=100000, mais pas pour M=10000"<<std::endl;
	std::cout<<"+ les barre d'erreur de plot.pdf doivent contenir P(x)"<<std::endl<<std::endl;

	Rand r(10);
	unsigned int const M(1e6);
	unsigned int const N(50);
	unsigned int n(0);
	unsigned int iter(0);
	double tol(1e-3);
	do{
		iter++;
		if(r.get(N)+1 <= n ){ n--; }/*because 0<=r.get(N)<N */
		else {n++;}
	} while (iter<M/5);

	double* H_wrong(new double [N]);
	Binning* H_right(new Binning[N]);
	H_right[24].log("log");
	iter=0;
	do{
		iter++;
		if(r.get(N)+1 <= n ){ n--; }/*because 0<=r.get(N)<N */
		else {n++;}
		H_right[n].add_sample(1.0);
		H_wrong[n]+=1.0;
		for(unsigned int i(0);i<n;i++){ H_right[i].add_sample(0.0); }
		for(unsigned int i(n+1);i<N;i++){ H_right[i].add_sample(0.0); }
		if(!(iter % (M/10))){
			H_right[24].is_converged(tol);
		}
	} while (iter<M);

	for(unsigned int i(0);i<N;i++){ H_right[i].is_converged(tol); }
	Write hist_right("hist_right.dat");
	Write hist_wrong("hist_wrong.dat");
	for(unsigned int i(0);i<N;i++){ 
		hist_right<<i<<" "<<H_right[i].get_mean()<<" "<<H_right[i].get_variance()<<Write::endl; 
		hist_wrong<<i<<" "<<H_wrong[i]/M<<" "<<sqrt((H_wrong[i]/M)*(1-H_wrong[i]/M)/(iter-1))<<Write::endl;
	}
	H_right[24].plot();

	Gnuplot gp("./","plot","plot");
	gp.preplot("N = "+tostring(N));
	gp.preplot("fac(x) = (int(x)==0) ? 1.0 : int(x) * fac(int(x)-1.0)");
	gp.preplot("P(x) = fac(N)/(2**N*fac(x)*fac(N-x))");
	gp.preplot("set xrange [10:40]");
	gp.preplot("set sample 31");
	gp.add("P(x) w p t '$P(x)$', 'hist_wrong.dat' w errorbars  t 'wrong', 'hist_right.dat' w errorbars t 'correct'");
	gp.save_file();
	gp.create_image(true);

	std::cout<<"x="<<H_right[24].get_mean()<<" ("<<H_right[24].get_variance()<<")"<<std::endl;

	delete[] H_right;
	delete[] H_wrong;
}

void check_flip_coin(){
	std::cout<<"Il est normal si des erreurs s'affichent. "<<std::endl<<std::endl;
	std::cout<<"+ si diff<1e-10, c'est simplement que la tolérance est trop ";
	std::cout<<"basse et que par conséquent, le premier test de convergence";
	std::cout<<"est vrai. Alors forcément xtot!=H.get_mean()."<<std::endl;
	std::cout<<"+ si 1e-10<diff<1e-16, c'est une approximation numérique"<<std::endl<<std::endl;

	Rand r(10000);
	srand(time(NULL));
	unsigned int N(10000);
	unsigned int iter(0);
	double tol(1e-5);

	double sup(0.0);
	double diff(0.0);

	Binning H(100,5);
	double xtot;
	double x;
	bool converged;
	for(unsigned int i(0);i<N;i++){
		iter=0;
		xtot=0.0;
		converged = false;
		do{ 
			iter++;
			//x=r.get();
			x=(1.0*rand())/RAND_MAX;
			H.add_sample(x); 
			xtot += x;
			if(iter>1e4){ converged=H.is_converged(tol);}
		} while (!converged);

		xtot/=iter;
		if(H.get_mean()>0.5){ sup+=1.0;} 
		/*{if the following test is wrong it is because the xtot and
		 * H_get_mean() are not computed for the same number of samples. this
		 * is normal because Binning recompute x_ and dx_ only when one bin is
		 * added in the smallest binning. So if one bin is partially "filled",
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
		if(std::abs(H.get_mean()-xtot)<1e-14){ diff+=1.0; } 
		else { std::cout<<"diff "<<std::abs(H.get_mean()-xtot)<<std::endl; }
		H.set();
	}

	std::cout<<std::endl<<"=> probabilité x in [0:0.5] = "<<sup/N*100<<"%"<<std::endl;
	std::cout<<"=> probabilité xtot==H.get_mean() = "<<diff/N*100<<"%"<<std::endl;
}
