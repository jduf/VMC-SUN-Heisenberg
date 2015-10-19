#ifndef DEF_OBSERVABLES
#define DEF_OBSERVABLES

#include "Sampling.hpp"

class Observable{
	public:
		/*!Constructor*/
		Observable(unsigned int const& nlinks, unsigned int const& nval, unsigned int const& B, unsigned int const& b, bool const& conv);
		Observable(unsigned int const& nlinks);
		/*!Constructor that reads from file*/
		Observable(IOFiles& r);
		/*!Copy constructor*/
		Observable(Observable const&) = default;
		/*!Default destructor*/
		virtual ~Observable() = default;
		Observable(Observable&&) = default;
		Observable& operator=(Observable obs);
		Observable() = delete;

		void set(unsigned int const& nval, unsigned int const& B, unsigned int const& b, bool const& conv);

		void merge(Observable & obs){ val_.merge(obs.val_); }
		void delete_binning(){ val_.delete_binning(); }
		void complete_analysis(double const& convergence_criterion){ val_.complete_analysis(convergence_criterion); }

		unsigned int nlinks() const { return links_.row(); }
		unsigned int nval() const { return val_.size(); }

		void set_x(double const& val);
		void add(unsigned int const& i, double const& val);
		void add_sample();

		Data<double> const& operator[](unsigned int const& i) const { return val_[i]; }

		int const& operator()(unsigned int const& i, unsigned int const& j) const { return links_(i,j); }
		int& operator()(unsigned int const& i, unsigned int const& j){ return links_(i,j); }

		void write(IOFiles& w) const;
		void set_links(Matrix<int> const& links){ links_ = links; }
		Matrix<int> const& get_links() const { return links_; }

	protected:
		Matrix<int> links_;
		DataSet<double> val_;
		unsigned int modulo_;

	private:
		void swap_to_assign(Observable& ds1, Observable& ds2);
};

IOFiles& operator<<(IOFiles& w, Observable const& obs);
IOFiles& operator>>(IOFiles& r, Observable& obs);
IOFiles& operator>>(IOFiles& r, Observable& obs);
std::ostream& operator<<(std::ostream& flux, Observable const& obs);
#endif
