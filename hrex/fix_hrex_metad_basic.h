#ifdef FIX_CLASS
// clang-format off
FixStyle(hrex/metad/basic, FixHREXMetadBasic);
// clang-format on
#else

#ifndef LMP_FIX_HREX_METAD_BASIC_H
#define LMP_FIX_HREX_METAD_BASIC_H

#include "fix_hrex_metad.h"
#include <vector>

namespace LAMMPS_NS {

class FixHREXMetadBasic : public FixHREXMetad {
public:
  FixHREXMetadBasic(class LAMMPS*, int, char**);
  ~FixHREXMetadBasic() override;
  
  void init() override;
  void post_force(int) override;
  void post_force_r(int);
  void post_force_xyz(int);
  
  CV root_collective_variables() override;
  Gaussian root_build_gaussian() override;
  void root_deposit_gaussian(Gaussian const&) override;
  double root_bias_potential(CV const&) override;
  std::vector<double> root_pack_potential() override;
  void root_unpack_potential(std::vector<double>&) override;
  
  void hrex_swap() override;
  
private:

  /**
   * Grid for storing the MetaD bias potential in 1, 2 or 3 dimensions.
   */
  struct Grid {
    enum { DISABLED, ENABLED }; 
    
    Grid()
      : nx(DISABLED)
      , ny(DISABLED)
      , nz(DISABLED)
      , potential()
      , neigh_list()
    {}

    void add_x_dim(double x0_, double x1_)
    {
      x0 = x0_;
      x1 = x1_;
      nx = ENABLED;
    }
    
    void add_y_dim(double y0_, double y1_)
    {
      y0 = y0_;
      y1 = y1_;
      ny = ENABLED;
    }

    void add_z_dim(double z0_, double z1_)
    {
      z0 = z0_;
      z1 = z1_;
      nz = ENABLED;
    }
    
    void init(double stdev)
    {
      neigh_list.clear();

      double approx_del = stdev / 10.; // rule of thumb
      rcut = 4. * stdev;
      rcut2 = rcut * rcut;
      scale = 1. / (2. * stdev * stdev);
      
      if(nx) {
	lx = x1 - x0;
	inv_lx = 1./lx;
	nx = (int)(lx / approx_del);
	dx = lx / nx;
	inv_dx = 1./dx;
      }
      if(ny) {
	ly = y1 - y0;
	inv_ly = 1./ly;
	ny = (int)(ly / approx_del);
	dy = ly / ny;
	inv_dy = 1./dy;
      }
      if(nz) {
	lz = z1 - z0;
	inv_lz = 1./lz;
	nz = (int)(lz / approx_del);
	dz = lz / nz;
	inv_dz = 1./dz;
      }
      
      int ir = nx ? (int)(rcut * inv_dx) : 0;
      int jr = ny ? (int)(rcut * inv_dy) : 0;
      int kr = nz ? (int)(rcut * inv_dz) : 0;
      for(int i=-ir; i<=ir; i++)
	for(int j=-jr; j<=jr; j++)
	  for(int k=-kr; k<=kr; k++)
	    neigh_list.push_back(std::make_tuple(i, j, k));

      int n=1;
      if(nx) n *= nx;
      if(ny) n *= ny;
      if(nz) n *= nz;
      potential=std::vector<double>(n, 0.);
    }

    void deposit_gaussian(Gaussian g)
    {
      pbc(g[0], g[1], g[2]);
      int i, j, k;
      xyz_to_ijk(g[0], g[1], g[2], i, j, k);
      for(auto it=neigh_list.cbegin(); it!=neigh_list.cend(); ++it) {
	double x, y, z;
	int ii = i + std::get<0>(*it);
	int jj = j + std::get<1>(*it);
	int kk = k + std::get<2>(*it);
	ijk_to_xyz(ii, jj, kk, x, y, z);
	double dx = g[0] - x;
	double dy = g[1] - y;
	double dz = g[2] - z;
	double dr2 = dx*dx + dy*dy + dz*dz;
	if(dr2 < rcut2)
	  {
	    pbc_ijk(ii, jj, kk);
	    potential[ijk_to_index(ii, jj, kk)] += g[3] * exp(-dr2*scale);
	  }
      }
    }

    double V(double x, double y, double z) const
    {
      // trilinear interpolation

      pbc(x, y, z);
      
      int i, j, k;
      xyz_to_ijk(x, y, z, i, j, k);

      double V000, V100, V110, V010,
	V001, V101, V111, V011;
      
      V000 = V_ijk(i, j, k);
      V100 = V_ijk(i+1, j, k);
      V110 = V_ijk(i+1, j+1, k);
      V010 = V_ijk(i, j+1, k);
      V001 = V_ijk(i, j, k+1);
      V101 = V_ijk(i+1, j, k+1);
      V111 = V_ijk(i+1, j+1, k+1);
      V011 = V_ijk(i, j+1, k+1);

      double xref, yref, zref;
      ijk_to_xyz(i, j, k, xref, yref, zref);

      double xd = (x - xref) * inv_dx;
      double yd = (y - yref) * inv_dy;
      double zd = (z - zref) * inv_dz;
      double V00 = V000 * (1. - xd) + V100 * xd;
      double V01 = V001 * (1. - xd) + V101 * xd;
      double V10 = V010 * (1. - xd) + V110 * xd;
      double V11 = V011 * (1. - xd) + V111 * xd;
      double V0 = V00 * (1. - yd) + V10 * yd;
      double V1 = V01 * (1. - yd) + V11 * yd;
      return V0 * (1. - zd) + V1 * zd;
    }
    
    double dV_dx(double x, double y, double z) const
    {
      if(!nx) return 0.;
      double V0 = V(x - dx, y, z);
      double V1 = V(x + dx, y, z);
      return (V1 - V0) * 0.5 * inv_dx;
    }
    
    double dV_dy(double x, double y, double z) const
    {
      if(!ny) return 0.;
      double V0 = V(x, y - dy, z);
      double V1 = V(x, y + dy, z);
      return (V1 - V0) * 0.5 * inv_dy;
    }

    double dV_dz(double x, double y, double z) const
    {
      if(!nz) return 0.;
      double V0 = V(x, y, z - dz);
      double V1 = V(x, y, z + dz);
      return (V1 - V0) * 0.5 * inv_dz;
    }

    double V_ijk(int i, int j, int k) const
    {
      pbc_ijk(i, j, k);
      return potential[ijk_to_index(i, j, k)];
    }

    void pbc(double& x, double& y, double& z) const
    {
      if(!nx) {
	x = 0;
      } else if(x < x0 || x >= x1) {
	int offset = (int)((x - x0) * inv_lx);
	if(x < x0) --offset;
	x -= offset * lx;
      }

      if(!ny) {
	y = 0;
      } else if(y < y0 || y >= y1) {
	int offset = (int)((y - y0) * inv_ly);
	if(y < y0) --offset;
	y -= offset * ly;
      }

      if(!nz) {
	z = 0;
      } else if(z < z0 || z >= z1) {
	int offset = (int)((z - z0) * inv_lz);
	if(z < z0) --offset;
	z -= offset * lz;
      }
    }

    void xyz_to_ijk(double x, double y, double z, int& i, int& j, int& k) const
    {
      i = nx ? (int)((x - x0) * inv_dx) : 0;
      j = ny ? (int)((y - y0) * inv_dy) : 0;
      k = nz ? (int)((z - z0) * inv_dz) : 0;
    }
    
    void ijk_to_xyz(int i, int j, int k, double& x, double& y, double& z) const
    {
      x = nx ? (x0 + i * dx) : 0.;
      y = ny ? (y0 + j * dy) : 0.;
      z = nz ? (z0 + k * dz) : 0.;
    }

    int ijk_to_index(int i, int j, int k) const
    {
      return i + j * (nx ? nx : 1) + k * (nx ? nx : 1) * (ny ? ny : 1);
    }

    void pbc_ijk(int& i, int& j, int& k) const
    {
      if(!nx) {
	i = 0;
      } else
	{
	  while(i < 0) i += nx;
	  while(i >= nx) i -= nx;
	}
      
      if(!ny) {
	j = 0;
      } else
	{
	  while(j < 0) j += ny;
	  while(j >= ny) j -= ny;
	}

      if(!nz) {
	k = 0;
      } else
	{
	  while(k < 0) k += nz;
	  while(k >= nz) k -= nz;
	}
    }

    std::vector<double> potential;
    double x0, x1;
    double y0, y1;
    double z0, z1;
    double lx, ly, lz;
    double inv_lx, inv_ly, inv_lz;
    int nx, ny, nz;
    double stdev;
    double approx_del;
    double dx, dy, dz;
    double inv_dx, inv_dy, inv_dz;
    double rcut, rcut2;
    double scale;
    std::vector<std::tuple<int,int,int>> neigh_list;
  };

  void reset_biasfactor();
  FILE* fopen_log(const char*) const;
  void load_record();
  void record_header() const;
  void record_gaussian(const Gaussian&) const;

  int nequ;
  bigint ntimestep_begin;
  int welltempering;
  double biasfactora, biasfactorb;
  double biastemp;
  double wtscale;
  double dumpscale;
  double gheight;
  double stdev;
  double kwallx, kwally, kwallz, kwallr;
  int jgroup, jgroupbit;
  double xlw, xuw;
  double ylw, yuw;
  double zlw, zuw;
  double rlw, ruw;
  int xflag, yflag, zflag, rflag;
  double masstotal, inv_masstotal;
  double masstotal2, inv_masstotal2;
  CV currentcv;
  char* record_filename;
  Grid root_grid;
};

} // LAMMPS_NS

#endif
#endif
