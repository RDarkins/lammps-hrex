// clang-format off

#include "fix_hrex_pot_srnb.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "math_const.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

FixHREXPotSRNB::FixHREXPotSRNB(LAMMPS* lmp, int narg, char** arg)
  : FixHREXPot(lmp, narg, arg)
  , unlambed_pe(0.)
{
  scalar_flag = 1;
  global_freq = 1;
  energy_global_flag = 1;

  // parse user input
  
  jgroup = group->find(utils::strdup(arg[3]));
  if (jgroup == -1) error->all(FLERR, "fix hrex/pot/srnb group ID does not exist");
  jgroupbit = group->bitmask[jgroup];
  scale_factor = utils::numeric(FLERR, arg[4], false, lmp);
}

double FixHREXPotSRNB::root_unlambed_pe()
{
  return unlambed_pe;
}

int FixHREXPotSRNB::setmask()
{
  return PRE_FORCE;
}

void FixHREXPotSRNB::init()
{
  if(!utils::strmatch(update->integrate_style,"^verlet")) {
    error->all(FLERR, "Fix hrex/pot/srnb does not support RESPA");
  }
  if (force->pair == nullptr) {
    error->all(FLERR, "No pair style defined for fix hrex/pot/srnb");
  }
  if (force->pair_match("^hybrid", 0) == nullptr && force->pair->single_enable == 0) {
    error->all(FLERR, "Pair style does not support fix hrex/pot/srnb");
  }
  neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
}

void FixHREXPotSRNB::init_list(int /*id*/, NeighList* ptr)
{
  list = ptr;
}

void FixHREXPotSRNB::setup_pre_force(int vflag)
{
  pre_force(vflag);
}

/**
 * Reduce the forces between two groups by a fraction
 * (scale_factor * lambda).
 * 
 * Add forces using pre_force so that contributions to
 * the ghost atoms are automatically communicated
 * by the subsequent pairwise force evaluations.
 */
void FixHREXPotSRNB::pre_force(int /*vflag*/)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, unlambed_eng, fpair, factor_coul, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  Pair* pair = force->pair;
  double** cutsq = pair->cutsq;
  double lambda = get_lambda();

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double my_unlambed_pe = 0;
  for(ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if(!(mask[i] & groupbit || mask[i] & jgroupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for(jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[j >> SBBITS & 3];
      factor_coul = special_coul[j >> SBBITS & 3];
      j &= NEIGHMASK;
      
      if(!(mask[j] & groupbit || mask[j] & jgroupbit)) continue;
      if(mask[i] & groupbit  && mask[j] & groupbit) continue;
      if(mask[i] & jgroupbit && mask[j] & jgroupbit) continue;
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if(rsq < cutsq[itype][jtype]) {
        unlambed_eng = - pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair)
	  * scale_factor;
	fpair *= - lambda * scale_factor;

	f[i][0] += delx * fpair;
	f[i][1] += dely * fpair;
	f[i][2] += delz * fpair;

        if(newton_pair || j < nlocal) {
	  my_unlambed_pe += unlambed_eng;
	  f[j][0] -= delx * fpair;
	  f[j][1] -= dely * fpair;
	  f[j][2] -= delz * fpair;
        } else {
          my_unlambed_pe += 0.5 * unlambed_eng;
        }
      }
    }
  }
  
  MPI_Allreduce(&my_unlambed_pe, &unlambed_pe, 1, MPI_DOUBLE, MPI_SUM, world);
}

double FixHREXPotSRNB::compute_scalar()
{
  return get_lambda() * unlambed_pe;
}
