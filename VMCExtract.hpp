#ifndef DEF_VMCEXTRACT
#define DEF_VMCEXTRACT

#include "VMCMinimization.hpp"

class VMCExtract : public VMCMinimization{
	public:
		VMCExtract(IOFiles& in, bool const& sort);
		/*!Default destructor*/
		virtual ~VMCExtract() = default;
		/*{Forbidden*/
		VMCExtract() = delete;
		VMCExtract(VMCExtract const&) = delete;
		VMCExtract(VMCExtract&&) = delete;
		VMCExtract& operator=(VMCExtract) = delete;
		/*}*/

		void refine(Vector<unsigned int> const& which_obs, double const& dEoE, unsigned int const& t, unsigned int maxiter = 0);
		void save(std::string const& filename) const;
		void print() const;
		List<MCSim>::Node* analyse(std::string const& path, std::string const& filename, List<MCSim>& kept_samples) const;

	private:
		class DiscardedSim{
			public:
				DiscardedSim(std::shared_ptr<MCSim> const& s):
					param_(s->get_param()),
					E_(s->get_energy()){
						E_.complete_analysis(1e-5);
						E_.delete_binning();
					}
				DiscardedSim(IOFiles& in):
					param_(in),
					E_(in){};

				void write(IOFiles& w) const { w<<param_<<E_; }
				static unsigned int sort(DiscardedSim const& list, DiscardedSim const& new_elem);
				Vector<double> const& get_param() const { return param_; }
				Data<double> const& get_energy() const { return E_; }

			private:
				Vector<double> param_;
				Data<double> E_;
		};

		List<VMCExtract::DiscardedSim> dis_sim_;
};
#endif
