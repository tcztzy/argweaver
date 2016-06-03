#include <utility>
#include "coal_records.h"
#include "local_tree.h"

namespace argweaver {


CoalRecord::CoalRecord(int start, int recomb_node, int recomb_time,
                       int coal_node, int coal_time, int pop_path) :
    start(start), recomb_node(recomb_node), recomb_time(recomb_time),
    coal_node(coal_node), coal_time(coal_time), pop_path(pop_path) {
    end = -1;
}

void CoalRecord::write(const char *chrom, FILE *file) {
    if (pop_path == -1) {  // just for debugging right now
        assert(0);
    }
    fprintf(file, "%s\t%i\t%i\t%i\t%i\t%i\t%i",
            chrom, start, end, recomb_node, recomb_time, coal_node,
            coal_time);
    if (pop_path != -1) fprintf(file, "\t%i", pop_path);
    fputc('\n', file);
}


/*
 Coal records written as BED files with columns
 chr, start, end, recomb node, recomb age, coal node, coal time, pop path

  pop path optional depending on whether model->pop_trees is NULL

  NEED TO THINK ABOUT
  Header at beginning of file defining times, pop tree, pop paths
  Do we want to recycle node numbers or have them constantly increase?
  constantly increasing would be consistent with McVean paper but not
  necessary with SMC assumption. Maybe use ARGweaver node numbers for now

 */

void write_coal_record_header(FILE *file, const ArgModel *model,
                              const LocalTrees *trees,
                              const Sequences &seqs) {
    /* first write header */
    int numleaf = trees->get_num_leaves();
    fprintf(file, "# ARGweaver coal_record\n");
    fprintf(file, "# ntimes %i\n", model->ntimes);
    fprintf(file, "# times");
    for (int i=0; i < model->ntimes; i++)
        fprintf(file, " %f", model->times[i]);
    fprintf(file, "\n");
    fprintf(file, "# num_sample %i\n", numleaf);
    fprintf(file, "# samples");
    for (int i=0; i < numleaf; i++) {
        if (i < (int)seqs.names.size())
            fprintf(file, " %s", seqs.names[trees->seqids[i]].c_str());
        else fprintf(file, " %d", trees->seqids[i]);
    }
    fprintf(file, "\n");
    if (model->pop_tree != NULL) {
        fprintf(file, "# npop %i\n", model->num_pops());
        fprintf(file, "# npath %i\n", model->num_pop_paths());
        for (int i=0; i < model->num_pop_paths(); i++) {
            fprintf(file, "# path%i", i);
            for (int j=0; j < model->ntimes; j++)
                fprintf(file, " %i", model->get_pop(i,j));
            fprintf(file, "\n");
        }
    }
}

typedef std::pair<int,int> intpair;
bool comparator (const intpair &l, const intpair &r) {
    return l.first < r.first;
}

void CoalRecords::writeAndPopCompleteRecords(FILE *file) {
    queue<int> tmp=order;
    vector<intpair> end_pos;
    while (! tmp.empty()) {
        int pos=0;
        int idx = tmp.front();
        int start_pos = records[idx].start;
        if (records[idx].end < 0) break;
        end_pos.clear();
        end_pos.push_back(intpair(records[idx].end, pos));
        tmp.pop();
        pos++;
        bool done_all=true;
        while (! tmp.empty()) {
            idx = tmp.front();
            if (records[idx].start != start_pos)
                break;
            if (records[idx].end < 0) {
                done_all = false;
                break;
            }
            end_pos.push_back(intpair(records[idx].end, pos));
            pos++;
            tmp.pop();
        }
        if (done_all) {
            sort(end_pos.begin(), end_pos.end(), comparator);
            vector<int> indexes(end_pos.size());
            for (unsigned int i=0; i < end_pos.size(); i++) {
                indexes[i] = order.front();
                order.pop();
            }
            for (unsigned int i=0; i < end_pos.size(); i++) {
                assert(end_pos[i].second >=0 && end_pos[i].second < (int)end_pos.size());
                idx = indexes[end_pos[i].second];
                records[idx].write(chrom, file);
                // NOTE: inefficient. could fix...
                for (int j=0; j < nnodes; j++) {
                    if (branch_created[j] == idx) branch_created[j]=-1;
                    if (node_created[j] == idx) node_created[j]=-1;
                }
                empty_records.push(idx);
            }
        } else break;
    }
}


void CoalRecords::addRecord(const Spr &spr, int pos, const LocalTree *last_tree,
                            int *mapping, FILE *outfile) {
    int recomb_node = spr.recomb_node;
    int coal_node = spr.coal_node;
    int broken_node = last_tree->nodes[recomb_node].parent;

    if (branch_created[recomb_node] != -1)
        records[branch_created[recomb_node]].end = pos;
    if (node_created[broken_node] != -1)
        records[node_created[broken_node]].end = pos;
    int pop_path = 0;
    if (model->pop_tree != NULL) {
        int path1 = model->consistent_path(last_tree->nodes[recomb_node].pop_path,
                                      spr.pop_path,
                                      last_tree->nodes[recomb_node].age,
                                      spr.recomb_time,
                                      spr.coal_time);
        int path2 = model->pop_tree->path_to_root(last_tree->nodes, coal_node,
                                                  spr.coal_time);
	pop_path = model->consistent_path(path1, path2, last_tree->nodes[recomb_node].age,
                                          spr.coal_time, model->ntimes - 1);
    }
    CoalRecord cr = CoalRecord(pos, mapping[recomb_node], spr.recomb_time,
                               mapping[coal_node], spr.coal_time, pop_path);
    int idx=-1;
    if (empty_records.size() >= 1) {
      idx = empty_records.front();
      empty_records.pop();
      records[idx] = cr;
    } else {
      idx = records.size();
      records.push_back(cr);
    }
    node_created[broken_node] = idx;
    order.push(idx);
    if (outfile != NULL)
        writeAndPopCompleteRecords(outfile);
}


 void CoalRecords::lastRecord(int end_pos, FILE *outfile) {
    int nbranch = (nnodes+1)/2;
    for (int i=0; i < nbranch; i++) {
        int idx = branch_created[i];
        if (idx >= 0) {
            records[idx].end = end_pos;
        }
    }
    for (int i=0; i < nnodes; i++) {
        int idx = node_created[i];
        if (idx >= 0) {
            records[idx].end = end_pos;
        }
    }
    if (outfile != NULL)
        writeAndPopCompleteRecords(outfile);
    assert(order.empty());
    for (int i=0; i < nbranch; i++)
        assert(branch_created[i] < 0);
    for (int i=0; i < nnodes; i++)
        assert(node_created[i] < 0);
}

CoalRecords::CoalRecords(const ArgModel *model, const LocalTree *first_tree,
                         int start_pos, const char *chrom) :
    model(model), chrom(chrom)
{
    nnodes = first_tree->nnodes;
    int nbranch = (nnodes+1)/2;
    branch_created = new int[nnodes];
    node_created = new int[nnodes];
    int i=0;
    for ( ; i < nbranch; i++) {
        branch_created[i] = i;
        node_created[i] = -1;
        int pop_path = ( model->pop_tree == NULL ? 0 :
                         model->pop_tree->path_to_root(first_tree->nodes, i) );
        CoalRecord cr(start_pos, i, -1, first_tree->nodes[i].parent,
                      first_tree->nodes[i].age, pop_path);
        records.push_back(cr);
        order.push(i);
    }
    for ( ; i < nnodes; i++) {
        node_created[i] = i;
        branch_created[i] = -1;
        int pop_path = ( model->pop_tree == NULL ? 0 :
                         model->pop_tree->path_to_root(first_tree->nodes, i) );
        CoalRecord cr(start_pos, i, -1, first_tree->nodes[i].parent,
                      first_tree->nodes[i].age, pop_path);
        records.push_back(cr);
        order.push(i);
    }
}


CoalRecords::~CoalRecords() {
    if (model != NULL) {
        delete [] branch_created;
        delete [] node_created;
    }
}


void write_coal_records(FILE *file, const ArgModel *model,
                        const LocalTrees *local_trees,
                        const Sequences *sequences, bool header) {

    LocalTrees::const_iterator it=local_trees->begin();
    int pos=local_trees->start_coord;
    CoalRecords cr(model, it->tree, pos, local_trees->chrom.c_str());
    int nnodes = it->tree->nnodes;
    int *total_mapping = new int[nnodes];
    int *tmp_mapping = new int[nnodes];
    for (int i=0; i < nnodes; i++)
        total_mapping[i]=i;
    pos += it->blocklen;
    const LocalTree *last_tree = it->tree;
    it++;
    if (header) write_coal_record_header(file, model, local_trees,
                                         sequences);
    while (it != local_trees->end()) {
        cr.addRecord(it->spr, pos, last_tree, total_mapping, file);
        pos += it->blocklen;
        last_tree = it->tree;
        int *mapping = it->mapping;
        for (int i=0; i < nnodes; i++)
            tmp_mapping[i] = total_mapping[i];
        for (int i=0; i < nnodes; i++) {
            if (mapping[i] != -1) {
                total_mapping[mapping[i]] = tmp_mapping[i];
            } else {
                int recoal = get_recoal_node(last_tree, it->spr, mapping);
                total_mapping[recoal] = tmp_mapping[i];
            }
        }
        it++;
    }
    cr.lastRecord(pos, file);
    delete [] total_mapping;
    delete [] tmp_mapping;
}


// read coal records and fill in trees and seqnames
// check that times match the ones in the model
bool read_coal_records(FILE *file, const ArgModel *model,
                        LocalTrees *trees, vector<string> &seqnames) {
    char *line = NULL;
    trees->clear();
    LocalTree *last_tree = NULL;
    Spr spr;
    spr.set_null();
    int numleaf=0;
    int nnodes = 0;
    int lineno = 0;
    int biggest_pos = -1;
    bool error=false;
    //    bool times_verified=false;
    bool ntimes_verified=false;
    int pos = -1;
    while ((line = fgetline(file))) {
        lineno++;
        chomp(line);
        if (strncmp(line, "# ntimes ", 9)==0) {
            int ntimes0;
            if (sscanf(&line[9], "%d", &ntimes0) != 1) {
                error=true;
                break;
            }
            if (ntimes0 != model->ntimes) {
                exitError("number of times in ARG file (%i) does not match ntimes in the model (%i)\n",
                          ntimes0, model->ntimes);
            }
            ntimes_verified=true;
        } else if (strncmp(line, "# times ", 8)==0) {
            assert(ntimes_verified);
            vector<string> times0;
            split(&line[8], " ", times0);
            if ((int)times0.size() != model->ntimes) {
                error=true;
                break;
            }
            for (int i=0; i < model->ntimes; i++) {
                double tmp = atof(times0[i].c_str());
                if (fabs(tmp - model->times[i]) > 0.001) {
                    exitError("time %i in ARG file (%f) does not match time in model (%f)\n",
                              i, tmp, model->times[i]);
                }
            }
            //            times_verified=true;
        } else if (strncmp(line, "# num_sample ", 13)==0) {
            if (sscanf(&line[13], "%d", &numleaf) != 1) {
                error=true;
                break;
            }
            nnodes = numleaf*2-1;
        } else if (strncmp(line, "# samples ", 10)==0) {
            assert(numleaf > 0);
            split(&line[10], " ", seqnames);
            assert((int)seqnames.size() == numleaf);
        } else if (strncmp(line, "# npop ", 7)==0) {
            int npop0;
            if (sscanf(&line[7], "%d", &npop0) != 1) {
                error=true;
                break;
            }
            if (npop0 != model->num_pops()) {
                exitError("Number of populations in ARG file (%i) does not match number in model (%i)\n",
                          npop0, model->num_pops());
            }
        } else if (strncmp(line, "# npath ", 8)==0) {
            int npath0;
            if (sscanf(&line[8], "%d", &npath0) != 1) {
                error=true;
                break;
            }
            if (npath0 != model->num_pop_paths()) {
                exitError("Number of population paths in ARG file (%i) does not match number in model (%i)\n",
                          npath0, model->num_pop_paths());
            }
            for (int i=0; i < npath0; i++) {
                assert(line=fgetline(file));
                lineno++;
                char tmp[100];
                sprintf(tmp, "# path%i", i);
                if (strncmp(line, tmp, strlen(tmp))!=0) {
                    error=true;
                }
                vector<string> path0;
                split(&line[strlen(tmp)+1], " ", path0);
                assert((int)path0.size() == model->ntimes);
                for (int j=0; j < model->ntimes; j++)
                    if (atoi(path0[j].c_str()) != model->get_pop(i, j)) {
                        exitError("Population paths in ARG file do not match paths in model\n");
                    }
            }
        } else if (strcmp(line, "# ARGweaver coal_record")==0) {
        } else if (strncmp(line, "#", 1)==0) {
            fprintf(stderr, "Unrecognized header line in ARG file '%s'\n", line);
        } else if (last_tree == NULL) {
            // first step is to read the first nnodes entries to form the first tree
            vector<string> fields;
            split(line, "\t", fields);
            assert((fields.size()==7 && model->pop_tree==NULL) ||
                   fields.size()==8);
            trees->start_coord = atoi(fields[1].c_str());
            biggest_pos = atoi(fields[2].c_str());
            int *ptree = new int[nnodes];
            int *ages = new int[nnodes];
            int *paths=NULL;
            for (int i=0; i < nnodes; i++) ptree[i]=-1;
            if (fields.size() == 8) paths = new int[nnodes];
            trees->chrom = fields[0].c_str();
            int node = atoi(fields[3].c_str());
            assert(node >= 0 && node < nnodes);
            assert(atoi(fields[4].c_str()) == -1);
            ptree[node] = atoi(fields[5].c_str());
            ages[node] = atoi(fields[6].c_str());
            if (fields.size()==8)
                paths[node] = atoi(fields[7].c_str());
            for (int i=1; i < nnodes; i++) {
                assert(line=fgetline(file));
                lineno++;
                split(line, "\t", fields);
                assert(strcmp(fields[0].c_str(), trees->chrom.c_str())==0);
                assert(atoi(fields[1].c_str())==trees->start_coord);
                if (atoi(fields[2].c_str()) > biggest_pos)
                    biggest_pos = atoi(fields[2].c_str());
                node = atoi(fields[3].c_str());
                assert(atoi(fields[4].c_str())==-1);
                ptree[node] = atoi(fields[5].c_str());
                ages[node] = atoi(fields[6].c_str());
                if (fields.size() == 8)
                    paths[node] = atoi(fields[7].c_str());
            }
            last_tree = new LocalTree(ptree, nnodes, ages, paths);
            spr.set_null();
            pos = trees->start_coord;
        } else {
            vector<string> fields;
            split(line, "\t", fields);
            assert((fields.size()==7 && model->pop_tree==NULL) ||
                   fields.size()==8);
            assert(strcmp(fields[0].c_str(), trees->chrom.c_str()) == 0);
            int new_pos = atoi(fields[1].c_str());
            int *mapping = NULL;
            if (!spr.is_null()) {
                mapping = new int [nnodes];
                for (int i=0; i < nnodes; i++)
                    mapping[i] = i;
                mapping[last_tree->nodes[spr.recomb_node].parent] = -1;
            }
            trees->trees.push_back(LocalTreeSpr(last_tree, spr, new_pos - pos,
                                                mapping));
            pos = new_pos;
            if (atoi(fields[2].c_str()) > biggest_pos)
                biggest_pos = atoi(fields[2].c_str());
            spr.recomb_node = atoi(fields[3].c_str());
            spr.recomb_time = atoi(fields[4].c_str());
            spr.coal_node = atoi(fields[5].c_str());
            spr.coal_time = atoi(fields[6].c_str());
            if (fields.size() == 8) {
                spr.pop_path = atoi(fields[7].c_str());
            } else spr.pop_path = 0;
            LocalTree *tree = new LocalTree(nnodes);
            tree->copy(*last_tree);
            apply_spr(tree, spr, model->pop_tree);
            last_tree = tree;
        }
        if (error) {
            exitError("Error reading ARG file on line %i\n", lineno);
        }
    }
    if (last_tree) {
        int *mapping = NULL;
        mapping = new int[nnodes];
        for (int i=0; i < nnodes; i++)
            mapping[i] = i;
        mapping[last_tree->nodes[spr.recomb_node].parent] = -1;
        trees->trees.push_back(LocalTreeSpr(last_tree, spr, biggest_pos - pos,
                                            mapping));
    }
    trees->end_coord = biggest_pos;
    //coords are 0-based in structure and file; no need to adjust

    trees->nnodes = nnodes;
    trees->set_default_seqids();
    assert_trees(trees);
    return true;
}

}  //namespace argweaver
