#ifndef DEF_PSO
#define DEF_PSO

#include"Rand.hpp"
#include"Read.hpp"
#include"Write.hpp"

class PSO {
	public:
		PSO(unsigned int Nbees, unsigned int Nfreedom, double cg, double cp, unsigned int maxiter);
		virtual ~PSO();

		void PSO_init();
		void PSO_set_limit(unsigned int param, double min, double max);
		void PSO_run(bool synchro=true);
		void PSO_print();
		void PSO_save(std::string filename);
		void PSO_load(std::string filename);

	protected:
		/*!This method must exist is the child class, it is the function that
		 * is minimized*/
		virtual double run(Vector<double> const& x)=0;

		unsigned int Nbees_; //!< numbre of bees
		unsigned int Nfreedom_; //!< number of degrees of freedom
		unsigned int maxiter_; //!< maximum number of iteration
		Vector<double>* pb_;//!< pb_[i] best position of the ith particle
		double* pfb_; //!< pfb_[i] value at the best position of the ith particle
		unsigned int bbee_; //!< best bee

	private:
		void launch(unsigned int i);
		void move(unsigned int i);
		void evaluate(unsigned int i);		

		Vector<double>* pv_;//!< pv_[i] velocity of the ith particle
		Vector<double>* px_;//!< px_[i] position of the ith particle
		Vector<double> min_;//!< min_(c) minimum value of the cth coordinate
		Vector<double> max_; //!< max(c) minimum value of the cth coordinate
		unsigned int forget_;
		bool* free_;//!< true if particle_[i] isn't running
		double cg_;
		double cp_;
		double chi_;
		Rand rnd_;
};
#endif
