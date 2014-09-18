#ifndef ARGWEAVER_MCMCMC_H
#define ARGWEAVER_MCMCMC_H

#include "logging.h"

#ifdef ARGWEAVER_MPI
//extern class MPI::Intracomm;
namespace MPI {
    class Intracomm;
}
#endif

namespace argweaver {

class Mc3Config
{
 public:

    Mc3Config() {
        group=0;
        max_group=0;
        heat_interval=0.05;
        heat=1.0;
    }

    Mc3Config(int group, double heat_interval);

    int group;
    int max_group;
    double heat_interval;
    double heat;
#ifdef ARGWEAVER_MPI
    MPI::Intracomm *group_comm;
#endif
};

} //namespace argweaver

#endif
