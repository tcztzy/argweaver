#ifdef ARGWEAVER_MPI
#include "mpi.h"
#endif

#include "logging.h"

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

 Mc3Config(int group, double heat_interval) :
   group(group), heat_interval(heat_interval) {
#ifdef ARGWEAVER_MPI
    int numthread=MPI::COMM_WORLD.Get_size();
    int *groups = (int*)malloc(numthread*sizeof(int));
    MPI::COMM_WORLD.Allgather(&group, 1, MPI::INT, groups, 1, MPI::INT);
    heat = 1.0 - heat_interval * group;
    group_comm = MPI::COMM_WORLD.Split(group, 0);

    //check that configuration makes sense.
    // Groups should be numbered 0.. max_group, there should be an equal
    // number of each
    int min_group = groups[0];
    max_group = groups[0];
    for (int i=0; i < numthread; i++) {
        if (groups[i] > max_group) max_group = groups[i];
        if (groups[i] < min_group) min_group = groups[i];
    }
    if (min_group != 0) {
        printError("Error: Should have a zero group for --mcmcmc-group");
        exit(1);
    }
    int counts[max_group+1];
    for (int i=0; i <=max_group; i++) counts[i]=0;
    for (int i=0; i < numthread; i++)
        counts[groups[i]]++;
    for (int i=1; i <= max_group; i++) {
        if (counts[i] != counts[0]) {
            printError("Error: Not all groups have same size for --mcmcmc-group");
            for (i=0; i < numthread; i++) printf("groups[%i]=%i\n", i, groups[i]); fflush(stdout);
            for (i=0; i <= max_group; i++) printf("counts[%i]=%i\n", i, counts[i]); fflush(stdout);
            exit(1);
        }
    }
    free(groups);
#endif
    }

    int group;
    int max_group;
    double heat_interval;
    double heat;
#ifdef ARGWEAVER_MPI
    MPI::Intracomm group_comm;
#endif
};

} //namespace argweaver
