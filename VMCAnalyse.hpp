#ifndef DEF_VMCANALYSE
#define DEF_VMCANALYSE

#include "VMCMinimization.hpp"

class VMCAnalyse : public VMCMinimization{
	public:
		VMCAnalyse(IOFiles& in);
		VMCAnalyse(List<IOFiles>& in);
		/*!Default destructor*/
		virtual ~VMCAnalyse() = default;
		/*{Forbidden*/
		VMCAnalyse() = delete;
		VMCAnalyse(VMCAnalyse const&) = delete;
		VMCAnalyse(VMCAnalyse&&) = delete;
		VMCAnalyse& operator=(VMCAnalyse) = delete;
		/*}*/

	private:
		void update_info(List<MCSim> const& merging);
};
#endif
