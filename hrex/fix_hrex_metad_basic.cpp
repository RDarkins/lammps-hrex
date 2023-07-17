#include "fix_hrex_metad_basic.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "comm.h"
#include "atom.h"
#include "domain.h"
#include "universe.h"
#include "update.h"

#include <iostream>

#define MAXLINE 1024

using namespace LAMMPS_NS;
using namespace FixConst;

FixHREXMetadBasic::FixHREXMetadBasic(LAMMPS* lmp, int narg, char** arg)
  : FixHREXMetad(lmp, narg, arg)
  , root_grid(), welltempering(0)
  , xflag(0)
  , yflag(0)
  , zflag(0)
  , rflag(0)
  , currentcv(3, 0.)
  , record_filename(nullptr)
  , nequ(0)
  , ntimestep_begin(0)
  , biasfactora(1.)
  , biasfactorb(1.)
{
  // parse user input

  if(narg < 9) {
    error->all(FLERR, "Illegal fix hrex/metad/basic command");
  }

  nequ = utils::inumeric(FLERR, arg[3], false, lmp);
  nevery = utils::inumeric(FLERR, arg[4], false, lmp);
  gheight = utils::numeric(FLERR, arg[5], false, lmp);
  stdev = utils::numeric(FLERR, arg[6], false, lmp);
  biastemp = utils::numeric(FLERR, arg[7], false, lmp);
  record_filename = strdup(arg[8]);
  int iarg = 9;
  while(iarg < narg) {
    if(strcmp(arg[iarg], "x") == 0) {
      xflag = 1;
      double x0 = utils::numeric(FLERR, arg[iarg+1], false, lmp);
      double x1 = utils::numeric(FLERR, arg[iarg+2], false, lmp);
      double skinx = utils::numeric(FLERR, arg[iarg+3], false, lmp);
      kwallx = utils::numeric(FLERR, arg[iarg+4], false, lmp);
      xuw = x1 - skinx;
      xlw = x0 + skinx;
      iarg += 5;
      if(comm->me == 0) {
	root_grid.add_x_dim(x0, x1);
      }
    } else if(strcmp(arg[iarg], "y") == 0) {
      yflag = 1;
      double y0 = utils::numeric(FLERR, arg[iarg+1], false, lmp);
      double y1 = utils::numeric(FLERR, arg[iarg+2], false, lmp);
      double skiny = utils::numeric(FLERR, arg[iarg+3], false, lmp);
      kwally = utils::numeric(FLERR, arg[iarg+4], false, lmp);
      yuw = y1 - skiny;
      ylw = y0 + skiny;
      iarg += 5;
      if(comm->me == 0) {
	root_grid.add_y_dim(y0, y1);
      }
    } else if(strcmp(arg[iarg], "z") == 0) {
      zflag = 1;
      double z0 = utils::numeric(FLERR, arg[iarg+1], false, lmp);
      double z1 = utils::numeric(FLERR, arg[iarg+2], false, lmp);
      double skinz = utils::numeric(FLERR, arg[iarg+3], false, lmp);
      kwallz = utils::numeric(FLERR, arg[iarg+4], false, lmp);
      zuw = z1 - skinz;
      zlw = z0 + skinz;
      iarg += 5;
      if(comm->me == 0) {
	root_grid.add_z_dim(z0, z1);
      }
    } else if(strcmp(arg[iarg], "r") == 0) {
      rflag = 1;
      jgroup = group->find(utils::strdup(arg[iarg+1]));
      if (jgroup == -1) error->all(FLERR, "Group ID does not exist");
      jgroupbit = group->bitmask[jgroup];
      double r0 = utils::numeric(FLERR, arg[iarg+2], false, lmp);
      double r1 = utils::numeric(FLERR, arg[iarg+3], false, lmp);
      double skinr = utils::numeric(FLERR, arg[iarg+4], false, lmp);
      kwallr = utils::numeric(FLERR, arg[iarg+5], false, lmp);
      ruw = r1 - skinr;
      rlw = r0 + skinr;
      iarg += 6;
      if(comm->me == 0) {
	// use the x dim of the grid to store r gaussians
	root_grid.add_x_dim(r0, r1);
      }
    } else if(strcmp(arg[iarg], "bf") == 0) {
      biasfactora = utils::numeric(FLERR, arg[iarg+1], false, lmp);
      biasfactorb = utils::numeric(FLERR, arg[iarg+2], false, lmp);
      iarg += 3;
    }
  }

  if(!rflag && !xflag && !yflag && !zflag) {
    error->all(FLERR, "Must enable at least one CV");
  }
  
  if(rflag && (xflag || yflag || zflag)) {
    error->all(FLERR, "Cannot use the r CV with the x/y/z CVs");
  }

  if(comm->me == 0) {
    root_grid.init(stdev);
  }
}

FixHREXMetadBasic::~FixHREXMetadBasic()
{
  free(record_filename);
}

void FixHREXMetadBasic::init()
{
  reset_biasfactor();
  load_record();
  
  ntimestep_begin = update->ntimestep + nequ;
  masstotal = group->mass(igroup);
  inv_masstotal = 1. / masstotal;
  if(rflag) {
    masstotal2 = group->mass(jgroup);
    inv_masstotal2 = 1. /  masstotal2;
  }
  
  record_header();
}

void FixHREXMetadBasic::reset_biasfactor()
{
  // linearly interpolate bias factors across replicas
  
  double biasfactor_i = biasfactora +
    (biasfactorb - biasfactora) * (double)get_ireplica() / (double)(get_nreplicas() - 1);
  
  if(biasfactor_i > 1.001) {
    welltempering = 1;
    wtscale = 1. / (biasfactor_i - 1.) / (force->boltz * biastemp);
    dumpscale = biasfactor_i / (biasfactor_i - 1.);
  } else {
    welltempering = 0;
    wtscale = 0.;
    dumpscale = 1.;
  }
}

void FixHREXMetadBasic::hrex_swap()
{
  reset_biasfactor();
}

void FixHREXMetadBasic::post_force(int vflag)
{
  if(rflag)
    post_force_r(vflag);
  else
    post_force_xyz(vflag);
}

/**
 * Force due to the distance (r) CV
 */
void FixHREXMetadBasic::post_force_r(int /*vflag*/)
{
  // distance between COM of groups
  
  double xcm[3],xcm2[3];
  group->xcm(igroup,masstotal,xcm);
  group->xcm(jgroup,masstotal2,xcm2);

  double dx = xcm2[0] - xcm[0];
  double dy = xcm2[1] - xcm[1];
  double dz = xcm2[2] - xcm[2];
  double r = MAX(sqrt(dx*dx + dy*dy + dz*dz), 1e-10);
  currentcv[0] = r;

  // metadynamics force

  double grad_V;
  if(comm->me == 0)
    grad_V = root_grid.dV_dx(currentcv[0], currentcv[1], currentcv[2]);
  MPI_Bcast(&grad_V, 1, MPI_DOUBLE, 0, world);
  
  double fmag = grad_V / r;

  if(r > ruw) {
    double del = r - ruw;
    fmag += kwallr * del / r;
  }
  
  double fx_per_m = fmag * inv_masstotal * dx;
  double fy_per_m = fmag * inv_masstotal * dy;
  double fz_per_m = fmag * inv_masstotal * dz;
  
  double fx2_per_m = fmag * inv_masstotal2 * dx;
  double fy2_per_m = fmag * inv_masstotal2 * dy;
  double fz2_per_m = fmag * inv_masstotal2 * dz;

  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double massone;
  
  if(rmass) {
    for(int i = 0; i < nlocal; i++) {
      if(mask[i] & groupbit) {
        massone = rmass[i];
        f[i][0] += fx_per_m*massone;
        f[i][1] += fy_per_m*massone;
        f[i][2] += fz_per_m*massone;
      }
      if(mask[i] & jgroupbit) {
        massone = rmass[i];
        f[i][0] -= fx2_per_m*massone;
        f[i][1] -= fy2_per_m*massone;
        f[i][2] -= fz2_per_m*massone;
      }
    }
  } else {
    for(int i = 0; i < nlocal; i++) {
      if(mask[i] & groupbit) {
        massone = mass[type[i]];
        f[i][0] += fx_per_m*massone;
        f[i][1] += fy_per_m*massone;
        f[i][2] += fz_per_m*massone;
      }
      if(mask[i] & jgroupbit) {
        massone = mass[type[i]];
        f[i][0] -= fx2_per_m*massone;
        f[i][1] -= fy2_per_m*massone;
        f[i][2] -= fz2_per_m*massone;
      }
    }
  }
}

/**
 * Force due to x/y/z CVs
 */
void FixHREXMetadBasic::post_force_xyz(int /*vflag*/)
{
  // current CV

  group->xcm(igroup,masstotal,currentcv.data());
  domain->remap(currentcv.data());
  
  // metadynamics force
  
  double grad_V[3];
  if(comm->me == 0) {
    grad_V[0] = xflag ? root_grid.dV_dx(currentcv[0], currentcv[1], currentcv[2]) : 0.;
    grad_V[1] = yflag ? root_grid.dV_dy(currentcv[0], currentcv[1], currentcv[2]) : 0.;
    grad_V[2] = zflag ? root_grid.dV_dz(currentcv[0], currentcv[1], currentcv[2]) : 0.;
  }
  MPI_Bcast(&grad_V[0], 3, MPI_DOUBLE, 0, world);
  
  double fx_per_m = - grad_V[0] * inv_masstotal;
  double fy_per_m = - grad_V[1] * inv_masstotal;
  double fz_per_m = - grad_V[2] * inv_masstotal;

  // force from confining walls

  double wfx = 0;
  double wfy = 0;
  double wfz = 0;

  if(xflag) {
    if(currentcv[0] > xuw) {
      double del = currentcv[0] - xuw;
      wfx -= kwallx * del;
    }
    if(currentcv[0] < xlw) {
      double del = currentcv[0] - xlw;
      wfx -= kwallx * del;
    }
  }
  
  if(yflag) {
    if(currentcv[1] > yuw) {
      double del = currentcv[1] - yuw;
      wfy -= kwally * del;
    }
    if(currentcv[1] < ylw) {
      double del = currentcv[1] - ylw;
      wfy -= kwally * del;
    }
  }

  if(zflag) {
    if(currentcv[2] > zuw) {
      double del = currentcv[2] - zuw;
      wfz -= kwallz * del;
    }
    if(currentcv[2] < zlw) {
      double del = currentcv[2] - zlw;
      wfz -= kwallz * del;
    }
  }

  fx_per_m += wfx * inv_masstotal;
  fy_per_m += wfy * inv_masstotal;
  fz_per_m += wfz * inv_masstotal;

  // apply the forces to the atoms

  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double massone;

  if(rmass) {
    for(int i = 0; i < nlocal; i++)
      if(mask[i] & groupbit) {
	massone = rmass[i];
	f[i][0] += fx_per_m * massone;
	f[i][1] += fy_per_m * massone;
	f[i][2] += fz_per_m * massone;
      }
  } else {
    for(int i = 0; i < nlocal; i++)
      if(mask[i] & groupbit) {
	massone = mass[type[i]];
	f[i][0] += fx_per_m * massone;
	f[i][1] += fy_per_m * massone;
	f[i][2] += fz_per_m * massone;
      }
  }
}

double FixHREXMetadBasic::root_bias_potential(CV const& cv)
{
  return root_grid.V(cv[0], cv[1], cv[2]);
}

CV FixHREXMetadBasic::root_collective_variables()
{
  return currentcv;
}

Gaussian FixHREXMetadBasic::root_build_gaussian()
{
  double h = gheight;
  if(welltempering) {
    double V = root_grid.V(currentcv[0], currentcv[1], currentcv[2]);
    h *= exp( - V * wtscale);
  }
  return { currentcv[0], currentcv[1], currentcv[2], h };
}

void FixHREXMetadBasic::root_deposit_gaussian(Gaussian const& g)
{
  if(update->ntimestep < ntimestep_begin) return;
  root_grid.deposit_gaussian(g);
  record_gaussian(g);
}

std::vector<double> FixHREXMetadBasic::root_pack_potential()
{
  return root_grid.potential;  
}

void FixHREXMetadBasic::root_unpack_potential(std::vector<double>& new_potential)
{
  root_grid.potential.swap(new_potential);
}

FILE* FixHREXMetadBasic::fopen_log(const char* mode) const
{
  char fn[MAXLINE];
  sprintf(fn, "%s.%d", record_filename, get_ireplica());
  return fopen(fn, mode);
}

/**
 * Reload the gaussian files if they already exist
 */
void FixHREXMetadBasic::load_record()
{
  static bool already_loaded = false;
  if(already_loaded) return;
  already_loaded = true;
  
  if(comm->me) return;
  
  FILE* fp;
  if(fp = fopen_log("r")) {

    // determine file format
    
    enum {X = 1 << 0, Y = 1 << 1, Z = 1 << 2};
    int format = 0;
    if(root_grid.nx) format |= X;
    if(root_grid.ny) format |= Y;
    if(root_grid.nz) format |= Z;
    
    // process each line
    
    char line[MAXLINE];
    while(std::fgets(line, MAXLINE, fp)) {
      if(line[0] == '#') continue;
      
      double x = 0, y = 0, z = 0, h;
      switch(format) {
	break;case X:
	sscanf(line, "%lf %lf", &x, &h);
	break;case X|Y:
	sscanf(line, "%lf %lf %lf", &x, &y, &h);
	break;case X|Z:
	sscanf(line, "%lf %lf %lf", &x, &z, &h);
	break;case X|Y|Z:
	sscanf(line, "%lf %lf %lf %lf", &x, &y, &z, &h);
	break;case Y:
	sscanf(line, "%lf %lf", &y, &h);
	break;case Y|Z:
	sscanf(line, "%lf %lf %lf", &y, &z, &h);
	break;case Z:
	sscanf(line, "%lf %lf", &z, &h);
      }
      
      root_grid.deposit_gaussian( {x, y, z, h / dumpscale} );
    }
    fclose(fp);
  }
}

void FixHREXMetadBasic::record_header() const
{
  if(comm->me || get_iwalker()) return;
  
  FILE* fp = fopen_log("a");
  fprintf(fp, "# %f", stdev);
  if(rflag) {
    fprintf(fp, " r %f %f", root_grid.x0, root_grid.x1);
  } else {
    if(root_grid.nx) fprintf(fp, " x %f %f", root_grid.x0, root_grid.x1);
    if(root_grid.ny) fprintf(fp, " y %f %f", root_grid.y0, root_grid.y1);
    if(root_grid.nz) fprintf(fp, " z %f %f", root_grid.z0, root_grid.z1);
  }
  fprintf(fp, "\n");
  fclose(fp);
}

void FixHREXMetadBasic::record_gaussian(const Gaussian& g) const
{
  if(comm->me || get_iwalker()) return;
  
  FILE* fp = fopen_log("a");
  if(root_grid.nx) fprintf(fp, "%f ", g[0]);
  if(root_grid.ny) fprintf(fp, "%f ", g[1]);
  if(root_grid.nz) fprintf(fp, "%f ", g[2]);
  fprintf(fp, "%f\n", g[3] * dumpscale);
  fclose(fp);
}
