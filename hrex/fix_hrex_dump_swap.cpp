#include "fix_hrex_dump_swap.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "hrex.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair.h"
#include "update.h"
#include "universe.h"

#include <iostream>

#define MAXLINE 1024

using namespace LAMMPS_NS;
using namespace FixConst;

FixHREXDumpSwap::FixHREXDumpSwap(LAMMPS* lmp, int nargs, char** arg)
  : FixHREXDump(lmp, nargs, arg)
{
  filename = new char[MAXLINE];
  if(nargs < 4) {
    error->all(FLERR, "Illegal fix hrex/dump/swap command");
  }
  sprintf(filename, "%s.%d", arg[3], universe->iworld);
}

FixHREXDumpSwap::~FixHREXDumpSwap()
{
  delete filename;
}

int FixHREXDumpSwap::restore_ireplica()
{
  // restore last ireplica if file already exists

  int ireplica = get_ireplica();
  if(fp = fopen(filename, "r")) {
    bigint ntimestep;
    int iwalker;
    double lambda;
    while(fscanf(fp, "%ld %d %d %lf\n", &ntimestep, &ireplica, &iwalker, &lambda) == 4) {
      if(comm->me==0 && universe->me==0)
	std::cerr<<ntimestep<<"\t"<<ireplica<<"\t"<<iwalker<<"\t"<<lambda<<std::endl;
    }
    fclose(fp);
  }

  return ireplica;
}

void FixHREXDumpSwap::hrex_init()
{
  fp = fopen(filename, "a");
  hrex_swap();
}

void FixHREXDumpSwap::hrex_end()
{
  if(fp) fclose(fp);
}

void FixHREXDumpSwap::hrex_swap()
{
  if(comm->me == 0) {
    fprintf(fp, "%ld %d %d %f\n", update->ntimestep, get_ireplica(), get_iwalker(), get_lambda());
    fflush(fp);
  }
}
