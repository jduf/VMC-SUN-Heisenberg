#ifndef DEF_OBSERVABLES
#define DEF_OBSERVABLES

#include "Sampling.hpp"

class Observable{
	public:
		/*!Constructor*/
		Observable(unsigned int const& nlinks, unsigned int const& nval, unsigned int const& B, unsigned int const& b, bool const& conv);
		/*Constructor that only sets links_*/
		Observable(Matrix<int> const& links);
		/*!Constructor that reads from file*/
		Observable(IOFiles& r);
		/*!Copy constructor*/
		Observable(Observable const&) = default;
		/*!Default move constructor*/
		Observable(Observable&&) = default;
		/*!Default destructor*/
		virtual ~Observable() = default;
		/*!Assignment operator that uses the copy-and-swap idiom*/
		Observable& operator=(Observable obs);
		/*{Forbidden*/
		Observable() = delete;
		/*}*/

		void set(unsigned int const& nval, unsigned int const& B, unsigned int const& b, bool const& conv);

		unsigned int nval() const { return val_.size(); }
		unsigned int nlinks() const { return links_.row(); }
		Matrix<int> const& get_links() const { return links_; }

		void merge(Observable & obs){ val_.merge(obs.val_); }
		void delete_binning(){ val_.delete_binning(); }
		void complete_analysis(double const& convergence_criterion){ val_.complete_analysis(convergence_criterion); }

		void set_x(double const& val){ for(unsigned int i(0);i<val_.size();i++){ val_[i].set_x(val); } }
		void add(unsigned int const& i, double const& val){ val_[links_(i,2)].add(val/modulo_); }
		void add_sample(){ val_.add_sample(); }

		Data<double> const& operator[](unsigned int const& i) const { return val_[i]; }

		int const& operator()(unsigned int const& i, unsigned int const& j) const { return links_(i,j); }
		int& operator()(unsigned int const& i, unsigned int const& j){ return links_(i,j); }

		void write(IOFiles& w) const;
		void print() const;

	protected:
		Matrix<int> links_;
		DataSet<double> val_;
		unsigned int modulo_ = 0;

	private:
		void swap_to_assign(Observable& ds1, Observable& ds2);
};

IOFiles& operator<<(IOFiles& w, Observable const& obs);
IOFiles& operator>>(IOFiles& r, Observable& obs);
IOFiles& operator>>(IOFiles& r, Observable& obs);
std::ostream& operator<<(std::ostream& flux, Observable const& obs);
#endif
