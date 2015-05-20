#ifndef DEF_MCSim
#define DEF_MCSim

#include "List.hpp"
#include "CreateSystem.hpp"
#include "SystemBosonic.hpp"
#include "SystemFermionic.hpp"
#include "MonteCarlo.hpp"

class MCSim {
	public:
		/*!Constructor that only sets param_*/
		MCSim(Vector<double> const& param);
		/*!Constructor that reads from file*/
		MCSim(IOFiles& r);
		/*!Default destructor*/
		~MCSim() = default;
		/*{Forbidden*/
		MCSim() = delete;
		MCSim(MCSim const&) = delete;
		MCSim(MCSim&&) = delete;
		MCSim& operator&=(MCSim) = delete;
		/*}*/

		static unsigned int cmp_for_fuse(MCSim const& list_elem, MCSim const& new_elem);
		static void fuse(MCSim& list_elem, MCSim& new_elem);

		Vector<double> const& get_param() const { return param_; }
		std::unique_ptr<MCSystem> const& get_S() const { return S_; }
		void create_S(Container* C, Vector<double> const* param);
		void copy_S(std::unique_ptr<MCSystem> const& S);

		void print() const { std::cout<<param_<<S_->get_energy(); }
		void write(IOFiles& w) const;
		bool is_created() const { return (S_.get()?!S_->get_status():false); }
		void save(Container* C) const;
		void run(double const& converged_criterion, unsigned int const& tmax);

	private:
		Vector<unsigned int> ref_;
		Vector<double> param_;
		std::unique_ptr<MCSystem> S_;
};
#endif
