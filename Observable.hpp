#ifndef DEF_OBSERVABLES
#define DEF_OBSERVABLES

#include "Sampling.hpp"
#include "Matrix.hpp"

class Observable{
	public:
		/*!Constructor*/
		Observable(unsigned int const& n, unsigned int const& B, unsigned int const& b, bool const& conv);
		/*!Constructor that reads from file*/
		Observable(IOFiles& r);
		/*!Copy constructor*/
		Observable(Observable const& s) = default;
		/*!Default destructor*/
		virtual ~Observable() = default;
		Observable(Observable&&) = default;
		Observable& operator=(Observable obs);

		void set(){ val_.set(); }
		void merge(Observable & obs){ val_.merge(obs.val_); }
		void delete_binning(){ val_.delete_binning(); }
		void complete_analysis(double const& convergence_criterion){ val_.complete_analysis(convergence_criterion); }
		
		unsigned int size() const { return links_.row(); }

		int const& operator()(unsigned int const& i, unsigned int const& j) const { return links_(i,j); }
		int& operator()(unsigned int const& i, unsigned int const& j){ return links_(i,j); }

		void set_x(unsigned int const& i, double const& val){ val_[i].set_x(val); }
		void add(unsigned int const& i, double const& val){ val_[i].add(val); }
		void add_sample(){ val_.add_sample(); }

		Data<double> const& operator[](unsigned int const& i) const { return val_[i]; }

	protected:
		/*!Default constructor*/
		Observable() = default;

		Matrix<int> links_;
		DataSet<double> val_;

		void swap_to_assign(Observable& ds1, Observable& ds2);
};

IOFiles& operator<<(IOFiles& w, Observable const& obs);
IOFiles& operator>>(IOFiles& r, Observable& obs);
#endif
