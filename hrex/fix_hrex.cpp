#include "fix_hrex.h"
#include "hrex.h"

using namespace LAMMPS_NS;

FixHREX::FixHREX(class LAMMPS* lmp, int narg, char** arg)
  : Fix(lmp, narg, arg)
  , hrex(nullptr)
{}

int FixHREX::get_ireplica() const
{
  return hrex ? hrex->get_ireplica() : 0;
}

int FixHREX::get_nreplicas() const
{
  return hrex ? hrex->get_nreplicas() : 1;
}

int FixHREX::get_iwalker() const
{
  return hrex ? hrex->get_iwalker() : 0;
}

int FixHREX::get_nwalkers() const
{
  return hrex ? hrex->get_nwalkers() : 1;
}

double FixHREX::get_lambda() const
{
  return hrex ? hrex->get_lambda() : 0.0;
}

/**
 * MPI wrapper for intra-walker gather. Must be called by rank 0 processes only
 */
int FixHREX::root_MPI_Allgather_walkers(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
				       void *recvbuf, int recvcount, MPI_Datatype recvtype)
{
  if(!hrex || hrex->iproc_world) return MPI_ERR_COMM;
  return MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, hrex->comm_walker_roots);
}

void FixHREX::set_hrex(HREX* t_hrex)
{
  hrex = t_hrex;
}
