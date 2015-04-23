#ifndef DEF_PSOFERMIONIC
#define DEF_PSOFERMIONIC

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"
#include "PSO.hpp"
#include "List.hpp"

class MCSim {
	public:
		/*!Constructor that only sets param_*/
		MCSim(Vector<double> const& param): param_(param) {}
		/*!Constructor that reads from file*/
		MCSim(IOFiles& r);
		/*!Destructor*/
		~MCSim(){}

		static unsigned int cmp_for_fuse(MCSim const& list_elem, MCSim const& new_elem);
		static void fuse(MCSim& list_elem, MCSim& new_elem);

		Vector<double> const& get_param() const { return param_; }
		std::unique_ptr<MCSystem> const& get_S() const { return S_; }
		void create_S(Container* C);
		void copy_S(std::unique_ptr<MCSystem> const& S);

		void print() const { std::cout<<param_<<S_->get_energy(); }
		void write(IOFiles& w) const;

	private:
		MCSim();
		MCSim(MCSim const&);

		Vector<unsigned int> ref_;
		Vector<double> param_;
		std::unique_ptr<MCSystem> S_;
};

class MCParticle: public Particle{
	public:
		MCParticle(){}
		~MCParticle(){}

		void move(Vector<double> const& bx_all);
		bool update(std::shared_ptr<MCSim> const& new_elem);

		void print() const;

	private:
		MCParticle(MCParticle const&);

		List<MCSim> history_;

		static void fuse(MCSim& list_elem, MCSim& new_elem);
};

class PSOFermionic: public Swarm<MCParticle>{
	public:
		PSOFermionic(Parseur* P);
		virtual ~PSOFermionic(){}

		void complete_analysis(double tol);

		void plot() const;
		void print() const;

	private:
		PSOFermionic();
		PSOFermionic(PSOFermionic const&);

		bool is_better_x(unsigned int const& p);
		void create();

		unsigned int tmax_;
		Container system_;
		List<MCSim> all_results_;
};
#endif
