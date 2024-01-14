#include "local_tree.h"
using namespace argweaver;
namespace argweaver {
extern "C" {
LocalTrees *read_local_trees(const char *filename, const double *times,
                             int ntimes);
void get_local_trees_blocks(const LocalTrees *trees, int *starts, int *ends);
void get_treelens(const LocalTrees *trees, const double *times, int ntimes,
                  double *treelens);
}
} // namespace argweaver
