#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(hrex, HREX);
// clang-format on
#else

#ifndef LMP_HREX_H
#define LMP_HREX_H

#include "command.h"
#include "fix_hrex.h"

namespace LAMMPS_NS {

class HREX : public Command {
  friend class FixHREX; // requires access to MPI comms

public:
  HREX(class LAMMPS*);
  ~HREX();

  void command(int, char**) override;

  int get_ireplica() const {return ireplica;}
  int get_nreplicas() const {return nreplicas;}
  int get_iwalker() const {return iwalker;}
  int get_nwalkers() const {return nwalkers;}
  double get_lambda() const {return ireplica_to_lambda[ireplica];}

private:
  void attempt_swap(int, int, int);
  void root_swap(int);
  void build_comms();
  void s_to_wr(int isim, int& iwalker, int& ireplica) const
  { iwalker = isim / nreplicas; ireplica = isim - iwalker * nreplicas; }
  void wr_to_s(int iwalker, int ireplica, int& isim) const
  { isim = ireplica + iwalker * nreplicas; }
  void print_status();
  
  int iproc_world;
  int nprocs_world;
  int iproc_universe;
  int nprocs_universe;
  int iworld;
  int nworlds;
  int isim;
  int nsims;
  int iwalker;
  int nwalkers;
  int ireplica;
  int nreplicas;
  int tempfix;
  int swap_counter;
  int glob_swap_counter;
  std::vector<class FixHREX*> fixes_hrex;
  std::vector<class FixHREXPot*> fixes_pot;
  std::vector<class FixHREXMetad*> fixes_metad;
  std::vector<class FixHREXDump*> fixes_dump;
  std::vector<double> ireplica_to_lambda;
  std::unique_ptr<class RanPark> prng;
  MPI_Comm comm_sim_roots;
  MPI_Comm comm_walker_roots;
};

} // LAMMPS_NS

#endif
#endif
