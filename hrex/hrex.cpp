#include "hrex.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "finish.h"
#include "fix.h"
#include "fix_hrex_pot.h"
#include "fix_hrex_metad.h"
#include "fix_hrex_dump.h"
#include "force.h"
#include "integrate.h"
#include "modify.h"
#include "random_park.h"
#include "timer.h"
#include "universe.h"
#include "update.h"

#include <iostream>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

HREX::HREX(LAMMPS* lmp)
  : Command(lmp)
  , fixes_hrex()
  , fixes_pot()
  , fixes_metad()
  , fixes_dump()
  , ireplica_to_lambda()
  , comm_sim_roots(MPI_COMM_NULL)
  , comm_walker_roots(MPI_COMM_NULL)
  , swap_counter(0)
  , glob_swap_counter(0)
{}

HREX::~HREX()
{
  if(comm_sim_roots != MPI_COMM_NULL) MPI_Comm_free(&comm_sim_roots);
  if(comm_walker_roots != MPI_COMM_NULL) MPI_Comm_free(&comm_walker_roots);
}

void HREX::command(int narg, char** arg)
{
  if(timer->is_timeout()) return;
  
  // parse user input
  
  if(domain->box_exist == 0) {
    error->universe_all(FLERR,"hrex command before simulation box is defined");
  }

  int nsteps = utils::inumeric(FLERR,arg[0],false,lmp);
  int nevery = utils::inumeric(FLERR,arg[1],false,lmp);
  for(tempfix = 0; tempfix < modify->nfix; tempfix++) {
    if(strcmp(arg[2],modify->fix[tempfix]->id) == 0) break;
  }
  if(tempfix == modify->nfix) {
    error->universe_all(FLERR,"Tempering fix ID is not defined");
  }
  int seed = utils::inumeric(FLERR,arg[3],false,lmp);
  nreplicas = utils::inumeric(FLERR,arg[4],false,lmp);
  for(int i=0; i<nreplicas; i++)
    ireplica_to_lambda.push_back(utils::numeric(FLERR, arg[5+i], false, lmp));
  
  int nswaps = nsteps / nevery;
  if(nevery * nswaps != nsteps)
    error->universe_all(FLERR, "Number of steps must be a multiple of the swap interval");

  // MPI setup
  
  iproc_world = comm->me;
  nprocs_world = comm->nprocs;
  iproc_universe = universe->me;
  nprocs_universe = universe->nprocs;
  iworld = universe->iworld;
  nworlds = universe->nworlds;
  isim = iworld;
  nsims = nworlds;
  nwalkers = nsims / nreplicas;
  if(nwalkers * nreplicas != nsims) {
    error->universe_all(FLERR, "Number of partitions must be a multiple of the number of replicas");
  }
  s_to_wr(isim, iwalker, ireplica);

  build_comms();

  // find all hrex fixes (pot, metad, dump)

  for(int i=0; i<modify->nfix; i++) {
    if(utils::strmatch(modify->fix[i]->style, "^hrex/pot/")) {
      auto pot_fix = dynamic_cast<FixHREXPot*>(modify->fix[i]);
      fixes_pot.push_back(pot_fix);
      fixes_hrex.push_back(pot_fix);
    } else if(utils::strmatch(modify->fix[i]->style, "^hrex/metad/")) {
      auto metad_fix = dynamic_cast<FixHREXMetad*>(modify->fix[i]);
      fixes_metad.push_back(metad_fix);
      fixes_hrex.push_back(metad_fix);
    } else if(utils::strmatch(modify->fix[i]->style, "^hrex/dump/")) {
      auto dump_fix = dynamic_cast<FixHREXDump*>(modify->fix[i]);
      fixes_dump.push_back(dump_fix);
      fixes_hrex.push_back(dump_fix);
    }
  }
  for(auto it=fixes_hrex.begin(); it!=fixes_hrex.end(); ++it) {
    (*it)->set_hrex(this);
  }
  
  if(nreplicas > 1 && fixes_pot.size() == 0 && iproc_world == 0) {
    error->warning(FLERR, "All replicas have the same Hamiltonian");
  }

  if(nwalkers > 1 && fixes_metad.size() == 0 && iproc_world == 0) {
    error->universe_all(FLERR, "Multiple walkers were allocated but MetaD is not being used");
  }

  // initialise random number generator

  prng = std::unique_ptr<RanPark>(new RanPark(lmp, seed + iproc_universe));
  for(int i = 0; i < 100; i++) prng->uniform();

  // setup for hrex run

  update->whichflag = 1;
  timer->init_timeout();
  update->nsteps = nsteps;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + nsteps;
  if(update->laststep < 0) {
    error->universe_all(FLERR,"Too many timesteps");
  }
  lmp->init();

  // restore previous state from dump

  for(auto it=fixes_dump.begin(); it!=fixes_dump.end(); ++it) {
    int new_ireplica = (*it)->restore_ireplica();
    int new_isim = new_ireplica + iwalker * nreplicas;
    if(new_ireplica != ireplica && iproc_world == 0) {
      root_swap(new_isim);
    }
    isim = new_isim;
    s_to_wr(isim, iwalker, ireplica);
    build_comms();  
  }
  MPI_Barrier(universe->uworld);

  for(auto it=fixes_hrex.begin(); it!=fixes_hrex.end(); ++it) {
    (*it)->hrex_init();
  }
  update->integrate->setup(1);
  timer->init();
  timer->barrier_start();
  
  if(iproc_universe == 0) {
    if(universe->uscreen) {
      fprintf(universe->uscreen, "Step\tNum swaps\n");
    }
    if(universe->ulogfile) {
      fprintf(universe->ulogfile, "Step\tNum swaps\n");
    }
    print_status();
  }
  
  for(int iswap=0; iswap<nswaps; iswap++) {
    
    // integrate

    timer->init_timeout();
    update->integrate->run(nevery);
    
    int my_timeout=0;
    int any_timeout=0;
    if(timer->is_timeout()) {
      my_timeout=1;
    }
    MPI_Allreduce(&my_timeout, &any_timeout, 1, MPI_INT, MPI_SUM, universe->uworld);
    if(any_timeout) {
      timer->force_timeout();
      break;
    }

    // replica exchange

    int i0 = nreplicas > 2 ? ((1 + iswap) % 2) : 0;
    for(int j = 0; j < nwalkers; j++) {
      for(int i = i0; i < nreplicas - 1; i += 2) {
	attempt_swap(j, i, i+1);
      }
    }

    MPI_Allreduce(&swap_counter, &glob_swap_counter, 1, MPI_INT, MPI_SUM, comm_sim_roots);
    if(iproc_universe == 0) print_status();
  }

  // cleanup
  
  timer->barrier_stop();
  for(auto it=fixes_hrex.begin(); it!=fixes_hrex.end(); ++it) {
    (*it)->hrex_end();
    (*it)->set_hrex(nullptr);
  }
  update->integrate->cleanup();
  Finish finish(lmp);
  finish.end(1);
  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}

void HREX::attempt_swap(int walkerab, int repa, int repb)
{
  const int prev_isim = isim;
  int sima, simb;
  wr_to_s(walkerab, repa, sima);
  wr_to_s(walkerab, repb, simb);
  bool involved = isim == sima || isim == simb; // is my proc involved in the swap?

  if(involved && iproc_world == 0) {
    double my_lambda = ireplica_to_lambda[ireplica];
    int partner_isim, partner_ireplica;
    double partner_lambda;
    if(isim == sima) {
      partner_isim = simb;
      partner_ireplica = repb;
      partner_lambda = ireplica_to_lambda[repb];
    } else {
      partner_isim = sima;
      partner_ireplica = repa;
      partner_lambda = ireplica_to_lambda[repa];
    }

    // delta pe contribution from this replica if we were to swap
    
    double delta_pe = 0;
    for(auto it=fixes_pot.begin(); it!=fixes_pot.end(); ++it)
      delta_pe += (*it)->root_unlambed_pe() * (partner_lambda - my_lambda);

    for(auto it=fixes_metad.begin(); it!=fixes_metad.end(); ++it) {
      std::vector<double> my_cv = (*it)->root_collective_variables();
      std::vector<double> partner_cv(my_cv.size());
	
      // swap cv with partner

      if(isim < partner_isim) {
	MPI_Recv(partner_cv.data(), my_cv.size(), MPI_DOUBLE, partner_isim,
		 0, comm_sim_roots, MPI_STATUS_IGNORE);
	MPI_Send(my_cv.data(), my_cv.size(), MPI_DOUBLE, partner_isim,
		 0, comm_sim_roots);
      } else {
	MPI_Send(my_cv.data(), my_cv.size(), MPI_DOUBLE, partner_isim,
		 0, comm_sim_roots);
	MPI_Recv(partner_cv.data(), my_cv.size(), MPI_DOUBLE, partner_isim,
		 0, comm_sim_roots, MPI_STATUS_IGNORE);
      }
      
      // delta bias potential contribution from this replica if we were to swap
      
      delta_pe += (*it)->root_bias_potential(partner_cv) - (*it)->root_bias_potential(my_cv);
    }
    
    // should we swap?
    
    int swap_flag;
    if(isim < partner_isim) {
      double partner_delta_pe;
      MPI_Recv(&partner_delta_pe, 1, MPI_DOUBLE, partner_isim,
	       0, comm_sim_roots, MPI_STATUS_IGNORE);
      delta_pe += partner_delta_pe;

      // Metropolis criterion
      
      int dim;
      double temperature = * (double*) modify->fix[tempfix]->extract("t_target", dim);
      swap_flag = delta_pe < 0 || prng->uniform() < exp(-delta_pe / (force->boltz * temperature));
      MPI_Send(&swap_flag, 1, MPI_INT, partner_isim, 0, comm_sim_roots);
    } else {
      MPI_Send(&delta_pe, 1, MPI_DOUBLE, partner_isim, 0, comm_sim_roots);
      MPI_Recv(&swap_flag, 1, MPI_INT, partner_isim, 0, comm_sim_roots, MPI_STATUS_IGNORE);
    }

    // swap
    
    if(swap_flag) {
      if(isim < partner_isim) {
	++swap_counter;
      }
      root_swap(partner_isim);
    }
  }
  
  // inform non-root procs of their (possibly new) isim
  
  if(involved) {
    MPI_Bcast(&isim, 1, MPI_INT, 0, world);
    s_to_wr(isim, iwalker, ireplica);
  }
  build_comms();  

  // notify hrex fixes if swap occurred
  
  if(involved && isim != prev_isim) {
    for(auto it = fixes_hrex.begin(); it != fixes_hrex.end(); ++it) {
      (*it)->hrex_swap();
    }
  }
}

void HREX::root_swap(int partner_isim)
{
  // swap the metad bias potentials
  
  for(auto it=fixes_metad.begin(); it!=fixes_metad.end(); ++it) {
    std::vector<double> packed_pot = (*it)->root_pack_potential();
    std::vector<double> partner_packed_pot(packed_pot.size());
      
    if(isim < partner_isim) {
      MPI_Send(packed_pot.data(), packed_pot.size(), MPI_DOUBLE,
	       partner_isim, 0, comm_sim_roots);
      MPI_Recv(partner_packed_pot.data(), packed_pot.size(), MPI_DOUBLE,
	       partner_isim, 0, comm_sim_roots, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(partner_packed_pot.data(), packed_pot.size(), MPI_DOUBLE,
	       partner_isim, 0, comm_sim_roots, MPI_STATUS_IGNORE);
      MPI_Send(packed_pot.data(), packed_pot.size(), MPI_DOUBLE,
	       partner_isim, 0, comm_sim_roots);
    }

    (*it)->root_unpack_potential(partner_packed_pot);
  }
  
  // swap isim
  
  isim = partner_isim;
}

void HREX::build_comms()
{
  if(comm_sim_roots != MPI_COMM_NULL) MPI_Comm_free(&comm_sim_roots);
  if(comm_walker_roots != MPI_COMM_NULL) MPI_Comm_free(&comm_walker_roots);
  
  MPI_Comm_split(universe->uworld,
		 std::min(iproc_world, 1),
		 isim, &comm_sim_roots);
  MPI_Comm_split(universe->uworld,
		 iproc_world == 0 ? ireplica : nreplicas,
		 isim, &comm_walker_roots);
}

void HREX::print_status()
{
  if(universe->uscreen) {
    fprintf(universe->uscreen, "%ld\t%d\n", update->ntimestep, glob_swap_counter);
  }
  if(universe->ulogfile) {
    fprintf(universe->ulogfile, "%ld\t%d\n", update->ntimestep, glob_swap_counter);
  }  
}
