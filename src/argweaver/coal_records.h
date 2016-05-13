//=============================================================================
// ArgHmm input/output of coalescence records


#ifndef ARGWEAVER_COAL_RECORDS_H
#define ARGWEAVER_COAL_RECORDS_H

// c/c++ includes
#include <math.h>
#include <set>
#include <map>
#include <list>
#include <algorithm>
#include <vector>
#include <string>
#include <queue>

// arghmm includes
#include "common.h"
#include "model.h"
#include "pop_model.h"
#include "local_tree.h"
#include "sequences.h"

namespace argweaver {

class CoalRecord {
 public:
    CoalRecord(int start,
               int recomb_node, int recomb_time,
               int coal_node, int coal_time,
               int pop_path=-1);

    void write(const char *chrom, FILE *file=stdout);

    int start;
    int end;
    int recomb_node;
    int recomb_time;
    int coal_node;
    int coal_time;
    int pop_path;
};

class CoalRecords {
 public:
    vector<CoalRecord> records;
    CoalRecords(const ArgModel *model, const LocalTree *first_tree,
                int start_pos, const char *chrom);
    ~CoalRecords();

    /* adds record to list. If outfile is not NULL, then write AND pop
       records as they are completed (i.e., once )
     */
    void addRecord(const Spr &spr, int pos, const LocalTree *last_tree,
                   int *mapping, FILE *outfile=NULL);

    void lastRecord(int end_pos, FILE *outfile=NULL);
    void writeAndPopCompleteRecords(FILE *file);

 private:
    queue<int> empty_records;  //indices in records vector which are free
    queue<int> order;       // order of indexes in records vector
    int *branch_created;    // which element of record each branch was created
    int *node_created;      // which element of record each node was created
    const ArgModel *model;
    const char *chrom;
    int nnodes;
};


void write_coal_record_header(FILE *file, const ArgModel *model,
                              const LocalTrees *trees, const Sequences &seqs);

void write_finished_coal_records(FILE *file, CoalRecords coal_records);

void write_coal_records(FILE *file, const ArgModel *model,
                        const LocalTrees *local_trees,
                        const Sequences *sequences,
                        bool header=true);

bool read_coal_records(FILE *file, const ArgModel *model,
                       LocalTrees *trees, vector<string> &seqnames);

/// TODO
//LocalTree *read_coal_records(FILE *file, int region_start, int region_end);

} //namespace argweaver

#endif  // ARGWEAVER_COAL_RECORDS
