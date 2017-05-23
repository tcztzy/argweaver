/*=============================================================================

  Matt Rasmussen
  Copyright 2010-2011

  Tree datastructure

  =============================================================================*/

// c++ headers
#include <assert.h>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <stdio.h>

// spidir headers
#include "Tree.h"
#include "common.h"
#include "parsing.h"
#include "logging.h"
#include "model.h"

namespace spidir {

using namespace argweaver;

bool isNewickChar(char c) {
    char vals[9] = "(),:#![]";
    return (c==vals[0] || c==vals[1] || c==vals[2] || c==vals[3] ||
            c==vals[4] || c==vals[5] || c==vals[6] || c==vals[7]);
}

//create a tree from a newick string
Tree::Tree(string newick, const ArgModel *model)
{
    int len = newick.length();
    Node *node = NULL;
    vector <int> stack;
    nnodes=0;
    int nbracket=0;
    for (int i=0; i < len; i++) {
        if (newick[i]=='[') nbracket++;
        else if (newick[i]==']') nbracket--;
        else if (newick[i]=='(' && nbracket==0) nnodes++;
    }
    nnodes += (nnodes+1);  //add in leaves
    nodes.setCapacity(nnodes);
    for (int i=0; i < nnodes; i++) {
        nodes[i] = new Node();
    }
    root = nodes[0];
    root->name = 0;
    stack.push_back(0);
    nnodes = 1;

    for (int i=0; i < len; i++) {
        switch (newick[i]) {
        case ',':
            stack.pop_back();
        case '(':
            node = nodes[nnodes];
            if (stack.size()==0) {
                printError("bad newick: error parsing tree");
                abort();
            } else {
                node->parent = nodes[stack.back()];
            }
            stack.push_back(nnodes);
            node->name = nnodes++;
            break;
        case ')': {
            stack.pop_back();
            node = nodes[stack.back()];
            break;
        }
        case ':':  { //optional dist next
            int j=i+1;
            while (j < len && !isNewickChar(newick[j]))
                j++;
            if (sscanf(&newick[i+1], "%lf", &node->dist) != 1) {
                printError("bad newick: error reading distance");
                abort();
            }
            i=j-1;
            break;
        }
        case '[': { // comment next; parse pop_path
            int count=1;
            int j=i+1;
            while (count != 0) {
                if (j==len) {
                    printError("bad newick: no closing bracket in NHX comment");
                    abort();
                }
                if (newick[j]==']') count--;
                else if (newick[j]=='[') count++;
                j++;
            }
            string tmpstr(&newick[i], j-i);
            int path_pos = tmpstr.find("pop_path=");
            if (path_pos >= 0) {
                // careful not to confuse pop_path with spr_pop_path
                if (path_pos < 4 || tmpstr.substr(path_pos-4, 4).compare("spr_") != 0) {
                    sscanf(tmpstr.substr(path_pos).c_str(), "pop_path=%i",
                           &(node->pop_path));
                }
            }
            i=j-1;
            break;
        }
        case ';':
            break;
        default:
            int j=i+1;
            while (j < len && !isNewickChar(newick[j]))
                j++;
            if (node->longname.length() > 0) {
                printError("bad newick format; got multiple names for a node");
                abort();
                break;
            }
            node->longname = newick.substr(i, j-i);
            trim(node->longname);
            i=j-1;
            break;
        }
    }
    if (node != root) {
        printError("bad newick format: did not end with root");
        abort();
    }
    //All done, now fill in children
    for (int i=0; i < nnodes; i++) {
        if (nodes[i]->parent != NULL) {
            nodes[i]->parent->addChild(nodes[i]);
        }
    }
    // do not assume a molecular clock! But assume that some samples have age 0, set
    // the age of the root to max distance from leaf to root, and all other ages accordingly
    ExtendArray<Node*> prenodes;
    getTreePreOrder(this, &prenodes);
    double maxage = 0.0;
    for (int i=0; i < prenodes.size(); i++) {
	if (prenodes[i]->nchildren==0) {
	    double tmpage=0.0;
	    Node *tmpnode = prenodes[i];
	    while (tmpnode->parent != NULL) {
		tmpage += tmpnode->dist;
		tmpnode = tmpnode->parent;
	    }
	    if (tmpage > maxage) maxage=tmpage;
	}
    }
    assert(prenodes[0] == root);
    root->age = maxage;
    for (int i=1; i < prenodes.size(); i++) {
	prenodes[i]->age = prenodes[i]->parent->age - prenodes[i]->dist;
    }
    if (model != NULL)
        this->correct_times(model, 1);

    for (int i=0; i < nnodes; i++) {
        if (nodes[i]->longname.length() > 0)
            nodename_map[nodes[i]->longname] = i;
    }
 }

double Tree::age_diff(double age1, double age2) {
    double diff = age1 - age2;
    if (diff < 0) {
        if (diff < -2) {
            fprintf(stderr, "got age diff=%.8f (age1=%.8f, age2=%.8f)\n",
                    diff, age1, age2);
            fflush(stderr);
            assert(0);
        }
        return 0.0;
    }
    return diff;
}

//void Tree::correct_times(map<string,double> times) {
void Tree::correct_times(const ArgModel *model, double tol) {
    unsigned int lasttime=0, j;
    ExtendArray<Node*> postnodes;
    getTreePostOrder(this, &postnodes);
    for (int i=0; i < postnodes.size(); i++) {
        if (postnodes[i]->nchildren == 0)
            lasttime = 0;
        double newage = postnodes[i]->age + postnodes[i]->dist;
        j = model->discretize_time(newage, lasttime, tol);
        postnodes[i]->dist = age_diff(model->times[j], postnodes[i]->age);
        if (postnodes[i]->parent != NULL)
            postnodes[i]->parent->age = model->times[j];
        lasttime = j;
    }
}


string Tree::format_newick_recur(Node *node, bool internal_names,
                                 char *branch_format_str,
                                 const NodeSpr *spr, bool oneline) {
    string rv;
    char tmp[1000];
    if (node->nchildren > 0) {
        int first=0, second=1;
        if (node->children[0]->longname.size() > 0 &&
            node->children[1]->longname.size() > 0 &&
            node->children[0]->longname.compare(node->children[1]->longname) > 0)
            {
                first=1;
                second=0;
            }
        rv.append("(");
        rv.append(format_newick_recur(node->children[first],
                                      internal_names,
                                      branch_format_str,
                                      spr, oneline));
        for (int i=1; i < node->nchildren; i++) {
            rv.append(",");
            rv.append(format_newick_recur(node->children[second],
                                          internal_names,
                                          branch_format_str,
                                          spr, oneline));
        }
        rv.append(")");
        if (internal_names) rv.append(node->longname);
    } else {
        rv.append(node->longname);
    }
    //    fprintf(f, "(%i)", node->name);
    if (branch_format_str != NULL && node->parent != NULL) {
        rv.append(":");
        sprintf(tmp, branch_format_str, node->dist);
        rv.append(tmp);
    }
    if (spr != NULL && ( node == spr->recomb_node || node == spr->coal_node)) {
        rv.append("[&&NHX");
        /*        if (have_pops) {
            sprintf(tmp, ":pop_path=%i", node->pop_path);
            rv.append(tmp);
            }*/
        if (spr != NULL) {
            if (node == spr->recomb_node) {
                sprintf(tmp, ":recomb_time=%.1f", spr->recomb_time);
                rv.append(tmp);
            }
            if (node == spr->coal_node) {
                sprintf(tmp, ":coal_time=%.1f", spr->coal_time);
                rv.append(tmp);
            }
        }
        rv.append("]");
    }
    if (!oneline && node->nchildren > 0) rv.append("\n");
    return rv;
}

string Tree::format_newick(bool internal_names, bool branchlen,
                           int num_decimal, const NodeSpr *spr,
                           bool oneline) {
    char *format_str=NULL;
    if (branchlen) {
        format_str = new char[100];
        snprintf(format_str, 100, "%%.%if", num_decimal);
    }
    string rv = format_newick_recur(root, internal_names,
                                    format_str, spr, oneline);
    rv.append(";");
    if (!oneline) rv.append("\n");
    if (format_str != NULL) {
        delete [] format_str;
    }
    return rv;
}


void Tree::write_newick(FILE *f, bool internal_name, bool branchlen,
                        int num_decimal, const NodeSpr *spr,
                        bool oneline) {
    string str = format_newick(internal_name, branchlen,
                               num_decimal, spr,
                               oneline);
    fprintf(f, "%s", str.c_str());
}

// return a copy of the tree
Tree *Tree::copy()
{
    Tree *tree2 = new Tree(nnodes);
    Node **nodes2 = tree2->nodes;

    for (int i=0; i<nnodes; i++) {
        nodes2[i]->setChildren(nodes[i]->nchildren);
        nodes2[i]->name = i;
        nodes2[i]->dist = nodes[i]->dist;
        nodes2[i]->longname = nodes[i]->longname;
        nodes2[i]->age = nodes[i]->age;
        nodes2[i]->pop_path = nodes[i]->pop_path;
    }

    for (int i=0; i<nnodes; i++) {
        for (int j=0; j<nodes[i]->nchildren; j++) {
            Node *child = nodes[i]->children[j];
            if (child)
                nodes2[i]->children[j] = nodes2[child->name];
            else
                nodes2[i]->children[j] = NULL;
        }
        Node *parent = nodes[i]->parent;
        if (parent)
            nodes2[i]->parent = nodes2[parent->name];
        else
            nodes2[i]->parent = NULL;
    }

    tree2->root = nodes2[root->name];
    tree2->nodename_map = nodename_map;
    return tree2;
}



//private function called by update_spr
//given a nhx tag within a newick string, returns the index of the node
//in the tree that the tag refers to.
// This is done assuming that only leafs have names, so first it finds
// the leaf below the node with the NHX tag, then counts the close
// parenthesis to figure out how many nodes up from the leaf to go.
int Tree::get_node_from_newick(char *newick, char *nhx) {
    int num_paren=0;
    char nhx_template[10]="[&&NHX";
    while (!str_starts_with(nhx, nhx_template)) {
        assert(nhx != newick);
        nhx--;
    }
    nhx--;
    while (1) {
        while (':' != nhx[0] && ')' != nhx[0]) {
            assert(nhx != newick);
            if (nhx[0]==']') {
                while (nhx[0]!='[') nhx--;
            }
            nhx--;
        }
        if (nhx[0]==':') nhx--;
        if (nhx[0]==')') {
            num_paren++;
            nhx--;
        } else {
            char *tmp=&nhx[1];
            assert(nhx[1]==':');
            nhx[1]='\0';
            while (!isNewickChar(nhx[0])) {
                nhx--;
            }
            nhx++;
            map<string,int>::iterator it = nodename_map.find(string(nhx));
            if (it == nodename_map.end()) { //error
                printf("nodename_map size=%i\n", (int)nodename_map.size());
                printf("nhx=%s\n", nhx);
                assert(it != nodename_map.end());
            }
            int n = it->second;
            assert(nodes[n]->nchildren==0);  // should be leaf
            tmp[0]=':';
            for (int i=0; i < num_paren; i++) {
                if (nodes[n]->parent == NULL) {
                    printf("assert failed newick=%s nhx=%s node=%i "
                           "num_paren=%i i=%i\n", newick, nhx, nodes[n]->name,
                           num_paren, i);
                }
                assert(nodes[n]->parent != NULL);
                n = nodes[n]->parent->name;
            }
            //printf("returning %i\n", n);
            return n;
        }
    }
}


void NodeMap::remap_node(Node *n, int id, int *deleted_branch) {
    int old_id = nm[n->name];
    //  printf("remap_node %i %i->%i\n", n->name, old_id, id);
    if (old_id == id) return; // no change
    inv_nm[old_id].erase(n->name);
    if (old_id >= 0 && inv_nm[old_id].size() == 0) {
        assert((*deleted_branch)==-1 || (*deleted_branch)==old_id);
        *deleted_branch = old_id;
        //    printf("deleted_branch=%i\n", *deleted_branch); fflush(stdout);
    }
    nm[n->name] = id;
    inv_nm[id].insert(n->name);
}


void NodeMap::propogate_map(Node *n, int *deleted_branch, int count,
                            int count_since_change, int maxcount,
                            int maxcount_since_change) {
    Node *c0, *c1;

    if (count==maxcount) return;
    if (count_since_change == maxcount_since_change) return;
    if (n->nchildren == 0)
        return propogate_map(n->parent, deleted_branch, count+1,
                             count_since_change+1, maxcount,
                             maxcount_since_change);
    c0 = n->children[0];
    c1 = n->children[1];
    /*  printf("propogate_map %i (%i) %i (%i) %i (%i) count=%i\n",
                n->name, nm[n->name], c0->name, nm[c0->name],
                c1->name, nm[c1->name], count); fflush(stdout); */
    int change=0;
    if (nm[c0->name] == -1 && nm[c1->name] == -1) {
        if (nm[n->name] != -1) {
            remap_node(n, -1, deleted_branch);
            change=1;
        }
    } else if (nm[c0->name] == -1 ||
               nm[c1->name] == -1) {
        Node *c = nm[c0->name] == -1 ? c1 : c0;
        if (nm[n->name] != nm[c->name]) {
            remap_node(n, nm[c->name], deleted_branch);
            change=1;
        }
    } else {  // neither are -1
        if (nm[c0->name] == nm[c1->name]) {
            printf("here n=%i c0=%i c1=%i %i %i %i %i\n",
                   n->name, c0->name, c1->name, nm[n->name],
                   nm[c0->name], nm[c1->name],
                   *deleted_branch);
            fflush(stdout);
            assert(0);
        }
        if (nm[n->name] == -1 || nm[n->name]==-3 ||
            nm[n->name] == nm[c0->name] ||
            nm[n->name] == nm[c1->name]) {
            change=1;
            remap_node(n, -2, deleted_branch);
        }
    }
    if (n->parent == NULL) return;
    return propogate_map(n->parent, deleted_branch, count+1,
                         change==0 ? count+1 : 0, maxcount,
                         maxcount_since_change);
}


//apply the spr operation to the tree.
//if node_map is not NULL, update it so that it maps to branches of
//pruned tree after SPR opreration on both trees
void Tree::apply_spr(NodeSpr *spr, NodeMap *node_map, const ArgModel *model) {
    Node *recomb_parent, *recomb_grandparent, *recomb_sibling,
        *coal_parent=NULL;
    int x;
    Node *recomb_node = spr->recomb_node;
    Node *coal_node = spr->coal_node;
    double coal_time = spr->coal_time;

    if (recomb_node == NULL) return;
    if (recomb_node == root)  assert(coal_node == root);
    if (recomb_node == coal_node) {
        recomb_node->pop_path = model->consistent_path(recomb_node->pop_path,
                                                       spr->pop_path,
                                                       recomb_node->age, spr->recomb_time,
                                                       ( recomb_node == root ? -1.0 :
                                                         recomb_node->parent->age));
        return;
    }

    recomb_parent = recomb_node->parent;
    assert(recomb_parent != NULL);  //ie, recomb_node should not be root
    assert(recomb_parent->nchildren == 2);
    x = (recomb_parent->children[0] == recomb_node ? 0 : 1);
    recomb_sibling = recomb_parent->children[!x];

    //recomb_grandparent might be NULL
    recomb_grandparent = recomb_parent->parent;

    //coal_parent might be NULL too
    coal_parent = coal_node->parent;

    //special case; topology doesn't change; just adjust branch lengths/ages
    if (coal_parent == recomb_parent) {
        assert(coal_node == recomb_sibling);
        if (model->pop_tree != NULL) {
            recomb_parent->pop_path =
                model->consistent_path(recomb_sibling->pop_path,
                                       coal_parent->pop_path,
                                       spr->coal_time,
                                       recomb_parent->age,
                                       recomb_parent == root ? -1 :
                                       recomb_parent->parent->age);
            recomb_node->pop_path =
                model->consistent_path(recomb_node->pop_path,
                                       spr->pop_path,
                                       recomb_node->age,
                                       spr->recomb_time, spr->coal_time);
        }
        coal_parent->age = coal_time;
        coal_node->dist = age_diff(coal_time, coal_node->age);
        recomb_node->dist = age_diff(coal_time, recomb_node->age);
        if (recomb_grandparent != NULL)
            recomb_parent->dist = age_diff(recomb_grandparent->age, coal_time);
        //fprintf(stderr, "done trivial update SPR\n");
        return;
    }
    // similar other special case
    if (coal_node == recomb_parent) {
        if (model->pop_tree != NULL) {
            recomb_node->pop_path =
                model->consistent_path(recomb_node->pop_path,
                                       spr->pop_path,
                                       recomb_node->age,
                                       spr->recomb_time, spr->coal_time);
            recomb_sibling->pop_path =
                model->consistent_path(recomb_sibling->pop_path,
                                       recomb_parent->pop_path,
                                       recomb_sibling->age,
                                       recomb_parent->age,
                                       spr->coal_time);
        }
        coal_node->age = coal_time;
        recomb_node->dist = age_diff(coal_time, recomb_node->age);
        recomb_sibling->dist = age_diff(coal_time, recomb_sibling->age);
        if (coal_parent != NULL)
            coal_node->dist = age_diff(coal_parent->age, coal_time);
        //fprintf(stderr, "done trivial update SPR2\n");
        return;
    }

    //now apply SPR
    if (model->pop_tree != NULL) {
        recomb_node->pop_path =
            model->consistent_path(recomb_node->pop_path,
                                   spr->pop_path,
                                   recomb_node->age, spr->recomb_time,
                                   spr->coal_time);
        recomb_sibling->pop_path =
            model->consistent_path(recomb_sibling->pop_path,
                                   recomb_parent->pop_path,
                                   recomb_sibling->age,
                                   recomb_parent->age,
                                   ( recomb_parent == root ? -1 :
                                     recomb_parent->parent->age ));
    }
    recomb_sibling->parent = recomb_grandparent;
    if (recomb_grandparent != NULL) {
        int x1 = (recomb_grandparent->children[0]==recomb_parent ? 0 : 1);
        recomb_grandparent->children[x1] = recomb_sibling;
        recomb_sibling->dist += recomb_parent->dist;
    } else {
        root = recomb_sibling;
        recomb_sibling->parent = NULL;
    }

    //recomb_parent is extracted; re-use as new_node. one child is still
    //recomb_node
    recomb_parent->children[!x] = coal_node;
    coal_node->dist = age_diff(coal_time, coal_node->age);
    recomb_node->dist = age_diff(coal_time, recomb_node->age);
    coal_node->parent = recomb_parent;
    recomb_parent->age = coal_time;
    recomb_parent->pop_path = coal_node->pop_path;
    if (coal_parent != NULL) {
        recomb_parent->parent = coal_parent;
        recomb_parent->dist = age_diff(coal_parent->age, coal_time);
        coal_parent->children[coal_parent->children[0]==coal_node ? 0 : 1] =
            recomb_parent;
    } else {
        root = recomb_parent;
        recomb_parent->parent = NULL;
    }
    //  this->write_newick(stdout, 1, NULL, 0); printf("\n"); fflush(stdout);
    if (node_map != NULL) {
        int deleted_branch=-1;
        //set recomb_node and recomb_parent maps to -3 = unknown
        node_map->remap_node(recomb_parent, -3, &deleted_branch);
        node_map->propogate_map(coal_node, &deleted_branch, 0, 0, 1, 1);
        node_map->propogate_map(recomb_node, &deleted_branch, 0, 0, 1, 1);
        node_map->propogate_map(recomb_sibling, &deleted_branch, 0, 0, -1, 4);
        node_map->propogate_map(recomb_parent, &deleted_branch, 0, 0, -1, 4);

        map<int,set<int> >::iterator it = node_map->inv_nm.find(-2);
        if (it != node_map->inv_nm.end() && it->second.size() > 0) {
            if (deleted_branch == -1) {
                assert(deleted_branch != -1);
            }
            set<int> rename_nodes = it->second;
            for (set<int>::iterator it2 = rename_nodes.begin();
                 it2 != rename_nodes.end(); ++it2) {
                node_map->remap_node(nodes[*it2], deleted_branch,
                                     &deleted_branch);
            }
        }
        for (int i=0; i < nnodes; i++) {
            if (node_map->nm[i] == -2) {
                assert(0);
            }
        }
        //    node_map->print();
    }
}




void NodeSpr::correct_recomb_times(const double *times, int ntimes) {
    int recomb_idx = time_index(recomb_time, times, ntimes);
    int coal_idx = time_index(coal_time, times, ntimes, recomb_idx);
    recomb_time = times[recomb_idx];
    coal_time = times[coal_idx];
}

/* get new SPR from newick string */
void NodeSpr::update_spr_from_newick(Tree *tree, char *newick,
                                     const ArgModel *model) {
    char search1[100]="recomb_time=";
    char search2[100]="coal_time=";
    char search3[100]="spr_pop_path=";
    char *x = strstr(newick, search1);
    if (x == NULL) {
        recomb_node = NULL;
        coal_node = NULL;
        return;
    }
    assert(1==sscanf(x, "recomb_time=%lg", &recomb_time));
    recomb_node = tree->nodes[tree->get_node_from_newick(newick, x)];

    x = strstr(newick, search2);
    assert(x != NULL);
    assert(1 == sscanf(x, "coal_time=%lg", &coal_time));
    coal_node = tree->nodes[tree->get_node_from_newick(newick, x)];

    x = strstr(newick, search3);
    if (x != NULL) {
        assert(1 == sscanf(x, "spr_pop_path=%i", &pop_path));
    } else pop_path=0;

    if (model != NULL) correct_recomb_times(model->times, model->ntimes);
}

//update the SPR on pruned tree based on node_map in big tree
void SprPruned::update_spr_pruned(const ArgModel *model) {
    if (orig_spr.recomb_node == NULL) {
        pruned_spr.recomb_node = pruned_spr.coal_node = NULL;
        return;
    }
    int num=node_map.nm[orig_spr.recomb_node->name];
    if (num == -1 || pruned_tree->nodes[num] == pruned_tree->root) {
        pruned_spr.recomb_node = pruned_spr.coal_node = NULL;
    } else {
        assert(num>=0);
        if (model != NULL && model->pop_tree != NULL) {
            int path1 = model->path_to_root(orig_spr.coal_node,
                                            orig_spr.coal_time);
            pruned_spr.pop_path = model->consistent_path(orig_spr.pop_path,
                                                         path1,
                                                         (double)orig_spr.recomb_time,
                                                         (double)orig_spr.coal_time,
                                                         -1.0);
        }
        pruned_spr.recomb_node = pruned_tree->nodes[num];
        pruned_spr.recomb_time = orig_spr.recomb_time;
        num = node_map.nm[orig_spr.coal_node->name];
        if (num == -1) {
            // coal node does not map; need to trace back until it does
            Node *n = orig_spr.coal_node;
            while (node_map.nm[n->name] == -1) {
                //root should always map to pruned tree
                assert(n->parent != NULL);
                n = n->parent;
            }
            assert(orig_spr.coal_time-1 <= n->age);
            pruned_spr.coal_time = n->age;
            pruned_spr.coal_node = pruned_tree->nodes[node_map.nm[n->name]];
        } else {
            assert(num >= 0);
            pruned_spr.coal_node = pruned_tree->nodes[num];
            pruned_spr.coal_time = orig_spr.coal_time;
        }
        if (pruned_spr.recomb_node == pruned_spr.coal_node) {
            pruned_spr.recomb_node = pruned_spr.coal_node = NULL;
        }
    }
    if (pruned_spr.recomb_node != NULL) assert(pruned_spr.coal_node != NULL);
}

bool NodeSpr::is_invisible() const {
    if (recomb_node == NULL || coal_node == NULL) return false;
    if (recomb_node == coal_node) return true;
    if (recomb_node->parent == NULL) return false;
    const Node *parent = recomb_node->parent;
    if (coal_node == parent &&
        fabs(coal_time - parent->age) < 0.001) {
        return true;
    }
    for (int i=0; i < parent->nchildren; i++) {
         if (parent->children[i] != recomb_node &&
            coal_node == parent->children[i] &&
            fabs(coal_time - parent->children[i]->age) < 0.001) {
            return true;
        }
    }
          return false;
}


void SprPruned::update(char *newick, const ArgModel *model) {
    //in this first case need to parse newick tree again
    if (orig_spr.recomb_node == NULL) {
        // this should only happen at the end of regions analyzed
        // by arg-sample, there will be no SPR there, so need to
        // parse the newick. The pruned tree can have NULL recomb_node
        // though because the SPR operation may have been pruned.
        update_slow(newick, model);
    } else {
        //otherwise, apply the SPR and node map, and get next SPR
        orig_tree->apply_spr(&orig_spr, inds.size() > 0 ? &node_map : NULL,
                             model);
        orig_spr.update_spr_from_newick(orig_tree, newick, model);
        if (pruned_tree != NULL) {
            if (pruned_spr.recomb_node != NULL)
                pruned_tree->apply_spr(&pruned_spr, NULL, model);
            update_spr_pruned(model);
        }
    }
}


NodeMap Tree::prune(set<string> leafs, bool allBut) {
    ExtendArray<Node*> newnodes = ExtendArray<Node*>(0);
    map<int,int> node_map;  //maps original nodes to new nodes
    ExtendArray<Node*> postnodes;
    getTreePostOrder(this, &postnodes);
    vector<bool> is_leaf(postnodes.size());
    node_map.clear();
    for (int i=0; i < postnodes.size(); i++) {
        is_leaf[i] = (postnodes[i]->nchildren == 0);
    }
    for (int i=0; i < postnodes.size(); i++) {
        if (postnodes[i]->nchildren == 0) {
            int prune;
            if (!is_leaf[i]) {
                //in this case, node was not originally a leaf, but now
                // has no children, so should be pruned
                prune = true;
            } else {
                prune = (leafs.find(postnodes[i]->longname) != leafs.end());
                if (allBut) prune=!prune;
            }
            if (prune) {
                Node *parent = postnodes[i]->parent;
                node_map[postnodes[i]->name] = -1;
                if (parent != NULL) {
                    int j, maxj=parent->nchildren;
                    for (j=0; j < maxj; j++) {
                        if (parent->children[j] == postnodes[i]) {
                            parent->children[j] =
                                parent->children[parent->nchildren-1];
                            parent->nchildren--;
                            delete postnodes[i];
                            break;
                        }
                    }
                    if (j == maxj) {
                        fprintf(stderr,
                                "error in tree.prune(): didn't find child\n");
                        exit(-1);
                    }
                } else {  //entire tree has been pruned!
                    for (i=0; i < nnodes; i++)
                        delete nodes[i];
                    nnodes=0;
                    root=NULL;
                }
            } else {
                node_map[postnodes[i]->name] = newnodes.size();
                newnodes.append(postnodes[i]);
            }
        } else if (postnodes[i]->nchildren == 1) {
            if (postnodes[i] == root) {
                node_map[postnodes[i]->name] =
                    node_map[postnodes[i]->children[0]->name];
                root = postnodes[i]->children[0];
                root->parent = NULL;
                delete postnodes[i];
            } else {
                Node *parent = postnodes[i]->parent;
                int j, maxj=parent->nchildren;
                for (j=0; j < maxj; j++) {
                    if (parent->children[j] == postnodes[i]) {
                        parent->children[j] = postnodes[i]->children[0];
                        postnodes[i]->children[0]->dist += postnodes[i]->dist;
                        postnodes[i]->children[0]->parent = parent;
                        node_map[postnodes[i]->name] =
                            node_map[postnodes[i]->children[0]->name];
                        delete postnodes[i];
                        break;
                    }
                }
                if (j == maxj) {
                    fprintf(stderr,
                            "error in tree.prune(): didn't find child2\n");
                    exit(-1);
                }
            }
        } else {
            node_map[postnodes[i]->name] = newnodes.size();
            newnodes.append(postnodes[i]);
        }
    }
    nodes.clear();
    for (int i=0; i < newnodes.size(); i++) {
        nodes.append(newnodes[i]);
        nodes[i]->name = i;
    }
    nnodes = nodes.size();
    nodename_map.clear();
    for (int i=0; i < nnodes; i++) {
        if (nodes[i]->longname.length() > 0)
            nodename_map[nodes[i]->longname] = i;
    }
    return NodeMap(node_map);
}

void SprPruned::update_slow(char *newick, const ArgModel *model) {
    if (orig_tree  != NULL) delete(orig_tree);
    if (pruned_tree != NULL) delete(pruned_tree);
    orig_tree = new Tree(newick, model);
    orig_spr = NodeSpr(orig_tree, newick, model);
    if ((int)inds.size() >= (orig_tree->nnodes+1)/2) {
        bool have_all_leafs = true;
        for (int i=0; i < orig_tree->nnodes; i++)
            if (orig_tree->nodes[i]->children == NULL &&
                inds.find(orig_tree->nodes[i]->longname) == inds.end()) {
                have_all_leafs = false;
                break;
            }
        if (have_all_leafs) inds.clear();
    }
    if (inds.size() > 0) {
        pruned_tree = orig_tree->copy();
        pruned_spr = orig_spr;
        node_map = pruned_tree->prune(inds, true);
        update_spr_pruned(model);
    } else pruned_tree = NULL;
}

// assumes both trees have same number of nodes
// and have same leaves
void Tree::setTopology(Tree *other)
{
    assert(nnodes == other->nnodes);
    Node **onodes = other->nodes;

    for (int i=0; i<nnodes; i++) {
        Node *node = nodes[i];
        Node *onode = onodes[i];

        if (onode->parent)
            node->parent = nodes[onode->parent->name];
        else
            node->parent = NULL;


        if (node->isLeaf()) {
            assert(onode->isLeaf());
        } else {
            // copy child structure
            nodes[i]->setChildren(onodes[i]->nchildren);
            for (int j=0; j<onodes[i]->nchildren; j++) {
                node->children[j] = nodes[onode->children[j]->name];
            }
        }
    }
}


// root tree by a new branch/node
void Tree::reroot(Node *newroot, bool onBranch)
{
    // handle trivial case, newroot is root
    if (root == newroot ||
        (onBranch &&
         root->nchildren == 2 &&
         (root->children[0] == newroot ||
          root->children[1] == newroot)))
        return;

    // determine where to stop ascending
    Node *oldroot = root;
    Node *stop1=NULL, *stop2=NULL;

    if (isRooted()) {
        stop1 = root->children[0];
        stop2 = root->children[1];
    } else {
        stop1 = root;
    }

    // start the reversal
    Node *ptr1 = NULL, *ptr2 = NULL;
    double nextDist = 0;
    double rootdist;

    if (onBranch) {
        if (isRooted()) {
            // just need to stick current root somewhere else
            Node *other = newroot->parent;
            rootdist = stop1->dist + stop2->dist;

            oldroot->children[0] = newroot;
            oldroot->children[1] = other;
            newroot->parent = oldroot;
            newroot->dist /= 2.0;

            ptr1 = other;

            int oldchild = find_array(ptr1->children, ptr1->nchildren, newroot);
            assert(oldchild != -1);

            // prepare for reversing loop
            ptr1->children[oldchild] = oldroot;
            ptr2 = oldroot;
            nextDist = newroot->dist;
        } else {
            // need to add a new node to be root
            // TODO: not implemented
            assert(0);
        }
    } else {
        if (isRooted()) {
            // need to remove the root node, and make tribranch
            // TODO: not implemented
            assert(0);
        } else {
            // just need to swap node positions
            // TODO: not implemented
            assert(0);
        }
    }


    // reverse parent child relationships
    while (ptr1 != stop1 && ptr1 != stop2) {
        int oldchild = find_array(ptr1->children, ptr1->nchildren, ptr2);
        assert(oldchild != -1);

        Node *next = ptr1->parent;

        // ptr1 is now fixed
        ptr1->children[oldchild] = next;
        ptr1->parent = ptr2;

        // swap distances
        double tmpdist = ptr1->dist;
        ptr1->dist = nextDist;
        nextDist = tmpdist;

        // move pointers
        ptr2 = ptr1;
        ptr1 = next;
    }


    // handle last two nodes
    if (stop2 != NULL) {
        // make stop1 parent of stop2
        if (stop2 == ptr1) {
            Node *tmp = stop1;
            stop1 = ptr1;
            stop2 = tmp;
        }
        assert(ptr1 == stop1);

        int oldchild = find_array(stop1->children, stop1->nchildren, ptr2);
        stop1->children[oldchild] = stop2;
        stop1->parent = ptr2;
        stop1->dist = nextDist;
        stop2->parent = stop1;
        stop2->dist = rootdist;
    } else {
        assert(0);
    }


    // renumber nodes
    // - all leaves don't change numbers
    assert(root->name = nnodes-1);
}


void Tree::reroot(Node *node1, Node *node2)
{
    // determine new root
    Node *newroot;
    if (node1->parent == node2)
        newroot = node1;
    else if (node2->parent == node1)
        newroot = node2;
    else if (node1->parent == root ||
             node2->parent == root)
        // do nothing
        return;
    else
        // not a valid branch
        assert(0);

    reroot(newroot);
}





// store a hash key representing the topology into the key array
// key is a parent tree representation where the internal nodes are
// given a consistent numbering
void Tree::hashkey(int *key)
{
    // get post order of nodes
    ExtendArray<Node*> postnodes;
    getTreePostOrder(this, &postnodes);

    // order children
    ExtendArray<int> ordering(nnodes);
    for (int i=0; i<postnodes.size(); i++)
        {
            Node *node=postnodes[i];

            if (node->isLeaf()) {
                ordering[node->name] = node->name;
            } else {
                // propogate the min order to the parent
                int minorder = ordering[node->children[0]->name];
                for (int j=1; j<node->nchildren; j++) {
                    int order = ordering[node->children[j]->name];
                    if (order < minorder)
                        minorder = order;
                }
                ordering[node->name] = minorder;
            }
        }

    // get a sorted post ordering of nodes
    ExtendArray<Node*> sortpostnodes;
    getTreeSortedPostOrder(this, &sortpostnodes, ordering);

    // generate a unique key for this topology
    // postfix notation for a tree
    // ((A,B),C) is represented as
    // A, B, -1, C, -1
    for (int i=0; i<sortpostnodes.size(); i++) {
        Node *node = sortpostnodes[i];

        if (node->isLeaf())
            key[i] = node->name;
        else
            key[i] = -1;
    }
}


bool Tree::sameTopology(Tree *other)
{
    if (other->nnodes != nnodes)
        return false;

    typedef ExtendArray<int> TopologyKey;
    TopologyKey key1(nnodes);
    TopologyKey key2(other->nnodes);

    hashkey(key1);
    other->hashkey(key2);

    for (int i=0; i<nnodes; i++) {
        if (key1[i] != key2[i])
            return false;
    }
    return true;
}


void Tree::reorderLeaves(string *order)
{
    // count the leaves in the tree
    int nleaves = 0;
    for (int i=0; i<nnodes; i++)
        if (nodes[i]->isLeaf())
            nleaves++;

    ExtendArray<Node*> tmp(nleaves);

    // rename leaves
    for (int i=0; i<nleaves; i++) {
        bool found = false;
        for (int j=0; j<nleaves; j++) {
            if (nodes[i]->longname == order[j]) {
                found = true;
                nodes[i]->name = j;
                tmp[j] = nodes[i];
                break;
            }
        }
        assert(found);
    }

    // reorder leaves by name
    for (int i=0; i<nleaves; i++)
        nodes[i] = tmp[i];
}


// assert that the tree datastructure is self-consistent
bool Tree::assertTree()
{
    if (root == NULL) {
        fprintf(stderr, "root == NULL\n");
        return false;
    }
    if (nnodes != nodes.size()) {
        fprintf(stderr, "nnodes != nodes.size()\n");
        return false;
    }
    if (root->parent != NULL) {
        fprintf(stderr, "root->parent != NULL\n");
        return false;
    }
    /*if (root->name != nnodes - 1) {
      fprintf(stderr, "root->name != nnodes - 1\n");
      return false;
      }*/

    bool leaves = true;
    for (int i=0; i<nnodes; i++) {
        //printf("assert %d\n", i);
        if (nodes[i] == NULL) {
            fprintf(stderr, "nodes[i] == NULL\n");
            return false;
        }

        // names are correct
        if (nodes[i]->name != i) {
            fprintf(stderr, "nodes[i]->name != i\n");
            return false;
        }

        // do leaves come first
        if (nodes[i]->isLeaf()) {
            if (!leaves) {
                fprintf(stderr, "!leaves\n");
                return false;
            }
        } else
            leaves = false;

        // check parent child pointers
        for (int j=0; j<nodes[i]->nchildren; j++) {
            //printf("assert %d %d\n", i, j);
            if (nodes[i]->children[j] == NULL) {
                fprintf(stderr, "nodes[i]->children[j] == NULL\n");
                return false;
            }
            //printf("assert %d %d parent\n", i, j);
            if (nodes[i]->children[j]->parent != nodes[i]) {
                fprintf(stderr, "nodes[i]->children[j]->parent != nodes[i]\n");
                return false;
            }
        }
    }

    //printf("done\n");

    return true;
}



void getTreeSortedPostOrder(Tree *tree, ExtendArray<Node*> *nodes,
                            int *ordering, Node *node)
{
    if (!node)
        node = tree->root;

    // make a child index array
    int childperm[node->nchildren];
    int childorder[node->nchildren];
    for (int i=0; i<node->nchildren; i++) {
        childperm[i] = i;
        childorder[i] = ordering[node->children[i]->name];
    }

    // sort index array by order
    ranksort(childperm, childorder, node->nchildren);

    // recurse
    for (int i=0; i<node->nchildren; i++)
        getTreeSortedPostOrder(tree, nodes, ordering,
                               node->children[childperm[i]]);

    // record post-process
    nodes->append(node);
}

void getTreePostOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node)
{
    if (!node)
        node = tree->root;

    // recurse
    for (int i=0; i<node->nchildren; i++)
        getTreePostOrder(tree, nodes, node->children[i]);

    // record post-process
    nodes->append(node);
}

void getTreePreOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node)
{
    if (!node)
        node = tree->root;

    // record pre-process
    nodes->append(node);

    // recurse
    for (int i=0; i<node->nchildren; i++)
        getTreePreOrder(tree, nodes, node->children[i]);
}





//=============================================================================
// primitive input/output


void printFtree(int nnodes, int **ftree)
{
    for (int i=0; i<nnodes; i++) {
        printf("%2d: %2d %2d\n", i, ftree[i][0], ftree[i][1]);
    }
}


// write out the names of internal nodes
void printTree(Tree *tree, Node *node, int depth)
{
    if (node == NULL) {
        if (tree->root != NULL) {
            printTree(tree, tree->root, 0);
            printf(";\n");
        }
    } else {
        if (node->nchildren == 0) {
            for (int i=0; i<depth; i++) printf("  ");
            printf("%d=%s:%lf", node->name, node->longname.c_str(), node->dist);
        } else {
            // indent
            for (int i=0; i<depth; i++) printf("  ");
            printf("%d=(\n", node->name);

            for (int i=0; i<node->nchildren - 1; i++) {
                printTree(tree, node->children[i], depth+1);
                printf(",\n");
            }

            printTree(tree, node->children[node->nchildren-1], depth+1);
            printf("\n");

            for (int i=0; i<depth; i++) printf("  ");
            printf(")");

            if (depth > 0)
                printf(":%lf", node->dist);
        }
    }
}


//=============================================================================
// Tree statistics

double Tree::total_branchlength() {
    double len=0.0;
    ExtendArray<Node*> postnodes;
    getTreePostOrder(this, &postnodes);
    for (int i=0; i < postnodes.size(); i++) {
        Node *node = postnodes[i];
        if (node != root) len += node->dist;
    }
    return len;
}



//Note: assumes all leaf nodes have same distance to root!
double Tree::tmrca() {
    return root->age;
}

// Returns an estimate of population size based on coalescence times
// in local tree
double Tree::avg_pairwise_distance() {
    ExtendArray<Node*> postnodes;
    int numDec[nnodes];
    int num_leaf = (nnodes+1)/2;
    double pi=0.0;
    getTreePostOrder(this, &postnodes);
    for (int i=0; i < postnodes.size(); i++) {
        if (postnodes[i]->parent == NULL) continue;
        int id = postnodes[i]->name;
        if (postnodes[i]->nchildren == 0)
            numDec[id] = 1;
        else {
            numDec[id] = 0;
            for (int j=0; j < postnodes[i]->nchildren; j++)
                numDec[id] += numDec[postnodes[i]->children[j]->name];
        }
        pi += postnodes[i]->dist * (double)(num_leaf - numDec[id])*numDec[id];
    }
    return pi*2.0/(num_leaf * (num_leaf-1));
}


double Tree::popsize() {
    int numleaf = (nnodes+1)/2;
    vector<double>ages;
    double lasttime=0, popsize=0;
    int k=numleaf;
    for (int i=0; i < nnodes; i++)
        if (nodes[i]->nchildren > 0)
            ages.push_back(nodes[i]->age);
    std::sort(ages.begin(), ages.end());
    for (unsigned int i=0; i < ages.size(); i++) {
        popsize += (double)k*(k-1)*(ages[i]-lasttime);
        lasttime = ages[i];
        k--;
    }
    return popsize/(4.0*numleaf-4);
}

//assume that times is sorted!
vector<double> Tree::coalCounts(const double *times, int ntimes) {
    vector<double> counts(ntimes, 0.0);
    vector<double> ages;
    unsigned int total=0;
    for (int i=0; i < nnodes; i++) {
        if (nodes[i]->nchildren > 0)
            ages.push_back(nodes[i]->age);
    }
    std::sort(ages.begin(), ages.end());
    int idx=0;
    for (int i=0; i < (int)ages.size(); i++) {
        while (1) {
            if (fabs(ages[i]-times[idx]) < 0.00001) {
                counts[idx]++;
                total++;
                break;
            }
            idx++;
            assert(idx < ntimes);
        }
    }
    assert(total == ages.size());
    return counts;
}


double Tree::num_zero_branches() {
    int count=0;
    for (int i=0; i < nnodes; i++) {
        if (nodes[i] != root && fabs(nodes[i]->dist) < 0.0001)
            count++;
    }
    return count;
}



double tmrca_half_rec(Node *node, int numnode, const vector<int> &numnodes) {
    if (node->nchildren != 2) {
        fprintf(stderr, "Error: tmrca_half only works for bifurcating trees\n");
    }
    if (numnodes[node->name] == numnode) return node->age;
    if (numnodes[node->children[0]->name] == numnode &&
        numnodes[node->children[1]->name] == numnode) {
        return min(node->children[0]->age, node->children[1]->age);
    }
    if (numnodes[node->children[0]->name] >= numnode) {
        assert(numnodes[node->children[1]->name] < numnode);
        return tmrca_half_rec(node->children[0], numnode, numnodes);
    } else if (numnodes[node->children[1]->name] >= numnode) {
        assert(numnodes[node->children[0]->name] < numnode);
        return tmrca_half_rec(node->children[1], numnode, numnodes);
    }
    return node->age;
}


double Tree::tmrca_half() {
    vector<int> numnodes(nnodes);
    ExtendArray<Node*> postnodes;
    getTreePostOrder(this, &postnodes);
    for (int i=0; i < postnodes.size(); i++) {
        numnodes[postnodes[i]->name] = 1;
        for (int j=0; j < postnodes[i]->nchildren; j++) {
            numnodes[postnodes[i]->name] +=
                numnodes[postnodes[i]->children[j]->name];
        }
    }
    assert(nnodes == numnodes[root->name]);
    return tmrca_half_rec(root, (nnodes-1)/2, numnodes);
}


double Tree::rth() {
    return this->tmrca_half()/this->tmrca();
}

double Tree::distBetweenLeaves(Node *n1, Node *n2) {
    if (n1 == n2) return 0.0;
    ExtendArray<Node*> postnodes;
    getTreePostOrder(this, &postnodes);
    vector<int> count(postnodes.size());
    int s=0;
    double rv=0.0;
    for (int i=0; i < postnodes.size(); i++) {
        if (postnodes[i] == n1 || postnodes[i] == n2) {
            count[postnodes[i]->name] = 1;
            s++;
        }
        if (postnodes[i]->nchildren == 2) {
            count[postnodes[i]->name] = count[postnodes[i]->children[0]->name] +
                count[postnodes[i]->children[1]->name];
            if (count[postnodes[i]->name] == 2) break;
        }
        if (count[postnodes[i]->name]) rv += postnodes[i]->dist;
    }
    assert(s == 2);
    return rv;
}


//look at descendants of parent node and return group number if all have
// same group, otherwise ngroup
int Tree::getDescGroups(Node *parent, map<string,int> groups,
                        int ngroup, int currgroup) {
    if (parent->nchildren==0) {
        map<string,int>::iterator it=groups.find(parent->longname);
        if (it == groups.end())
            return currgroup;
        int group = it->second;
        if (currgroup == group || currgroup == -1)
            return group;
        return ngroup;
    }
    for (int i=0; i < parent->nchildren; i++) {
        currgroup = getDescGroups(parent->children[i], groups,
                                  ngroup, currgroup);
        if (currgroup == ngroup) return ngroup;
    }
    return currgroup;
}

void Tree::countDescGroups(Node *node, string hap, map<string,int> groups,
                           int ngroups, double addval, vector<double> &rv) {
    Node *parent = node->parent;
    assert(parent != NULL);
    Node *sib = node->parent->children[0];
    if (sib == node) sib = node->parent->children[1];
    Node *grandparent = parent->parent;
    Node *aunt=NULL;
    if (grandparent != NULL) {
        aunt = grandparent->children[0];
        if (aunt == parent) aunt = grandparent->children[1];
    }

    int coalgroup_sib = getDescGroups(sib, groups, ngroups);
    int coalgroup_aunt = aunt != NULL ?
        getDescGroups(aunt, groups, ngroups) : -1;
    int this_subgroup = groups[hap];
    if (coalgroup_aunt < 0) {
        if (this_subgroup >= 0 && this_subgroup == coalgroup_sib)
            rv[this_subgroup] += addval;
        else rv[ngroups] += addval;
    } else {
        if (coalgroup_aunt == this_subgroup)
            rv[this_subgroup] += addval;
        else if (coalgroup_aunt == coalgroup_sib)
            rv[coalgroup_sib] += addval;
        else rv[ngroups] += addval;
    }
}

vector<double> Tree::coalGroup(string hap1, string hap2,
                               map<string,int> groups, int numgroup) {
    int node1, node2;
    bool nodesTogether=true;
    vector<double> rv(numgroup+1, 0);
    double val=1.0;

    map <string,int>::iterator it = nodename_map.find(hap1);
    if (it == nodename_map.end()) {
        printError("No leaf named %s", hap1.c_str());
        abort();
    }
    node1 = it->second;
    it = nodename_map.find(hap2);
    if (it == nodename_map.end())
        node2=-1;
    else node2 = it->second;

    Node *parent = nodes[node1]->parent;
    Node *sib = parent->children[0];
    if (sib == nodes[node1])
        sib = parent->children[1];
    else assert(parent->children[1] == nodes[node1]);

    Node *node=NULL;
    if (node2 != -1 && sib == nodes[node2]) {
        node = parent;
        parent = node->parent;
        nodesTogether=true;
        val = 1.0;
    } else if (node2 != -1) {
        val = 0.5;
        nodesTogether=false;
        node = nodes[node1];
    }
    if (node != NULL)
	countDescGroups(node, hap1, groups, numgroup, val, rv);
    if (node2 != -1 && !nodesTogether)
        countDescGroups(nodes[node2], hap2, groups, numgroup, val, rv);
    return rv;
}


bool Tree::isGroup(set<string> group) {
    if (group.size() <= 1) return true;
    ExtendArray<Node*> postnodes;
    getTreePostOrder(this, &postnodes);
    vector<int> ingroup(postnodes.size());
    vector<int> outgroup(postnodes.size());
    for (int i=0; i < postnodes.size(); i++) {
        if (postnodes[i]->nchildren == 0) {
            ingroup[postnodes[i]->name] =
                ( group.find(postnodes[i]->longname) != group.end() );
            outgroup[postnodes[i]->name] = !ingroup[postnodes[i]->name];
        } else {
            int curr = postnodes[i]->name;
            int child[2];
            int isroot = (postnodes[i] == root);
            child[0] = postnodes[i]->children[0]->name;
            child[1] = postnodes[i]->children[1]->name;

            ingroup[curr] = ingroup[child[0]] + ingroup[child[1]];
            outgroup[curr] = outgroup[child[0]] + outgroup[child[1]];
            if (ingroup[curr] == (int)group.size()) {
                if (outgroup[curr] == 0) return true;
                if (!isroot) return false;
            }
            if (ingroup[curr] > 0) {
                if (outgroup[child[0]] > 0 && outgroup[child[1]] > 0) return false;
            }
        }
    }
    //got to root
    assert(ingroup[root->name] == (int)group.size());
    if (outgroup[root->children[0]->name] == 0 ||
        outgroup[root->children[1]->name] == 0) return true;
    return false;
}

// want to return set of nodes above which mutations happened under infinite
// sites to cause site pattern.
// assumes tree has been pruned to remove non-informative leafs
// this is not particularly efficient! may want to think of something faster.
set<Node*> Tree::lca(set<Node*> derived) {
    set<Node*>::iterator it;
    set<Node*> rv;

    if (derived.size() == 1) return derived;
    ExtendArray<Node*> postnodes;
    getTreePostOrder(this, &postnodes);
    for (int i=0; i < postnodes.size(); i++) {
        if (postnodes[i]->nchildren == 0) continue;
        if (postnodes[i] == root) {
            assert(derived.size() == 1);
            rv.insert(*(derived.begin()));
            return rv;
        }
        int count=0;
        for (int j=0; j < postnodes[i]->nchildren; j++) {
            if (derived.find(postnodes[i]->children[j]) != derived.end())
                count++;
        }
        if (count == postnodes[i]->nchildren) { // all children are derived
            for (int j=0; j < postnodes[i]->nchildren; j++)
                derived.erase(postnodes[i]->children[j]);
            derived.insert(postnodes[i]);
        }  else if (count != 0) {
            for (int j=0; j < postnodes[i]->nchildren; j++) {
                if (derived.find(postnodes[i]->children[j]) != derived.end()) {
                    rv.insert(postnodes[i]->children[j]);
                    derived.erase(postnodes[i]->children[j]);
                }
            }
        }
        if (derived.size() == 0) return rv;
    }
    fprintf(stderr, "got to end of LCA\n"); fflush(stderr);
    return rv;
}

bool Tree::haveMig(int p[2], int t[2], const ArgModel *model) {
    assert(t[0] < t[1]);
    double dt[2];
    dt[0] = model->times[t[0]];
    dt[1] = model->times[t[1]];
    for (int i=0; i < nnodes; i++) {
        if (nodes[i]->age-2 <= dt[0] &&
            (nodes[i]->parent == NULL ||
             nodes[i]->parent->age+2 >= dt[1])) {
            if (model->get_pop(nodes[i]->pop_path, t[0]) == p[0] &&
                model->get_pop(nodes[i]->pop_path, t[1]) == p[1])
                return true;
        }
    }
    return false;
}

//=============================================================================
// primitive tree format conversion functions

extern "C" {

/*
// creates a forward tree from a parent tree
// Note: assumes binary tree
void makeFtree(int nnodes, int *ptree, int ***ftree)
{
    *ftree = new int* [nnodes];
    int **ftree2 = *ftree;

    // initialize
    for (int i=0; i<nnodes; i++) {
        ftree2[i] = new int [2];
        ftree2[i][0] = -1;
        ftree2[i][1] = -1;
    }

    // populate
    for (int i=0; i<nnodes; i++) {
        int parent = ptree[i];

        if (parent != -1) {
            if (ftree2[parent][0] == -1)
                ftree2[parent][0] = i;
            else
                ftree2[parent][1] = i;
        }
    }
}


void freeFtree(int nnodes, int **ftree)
{
    for (int i=0; i<nnodes; i++)
    delete [] ftree[i];
    delete [] ftree;
}
*/


// create a tree object from a parent tree array
void ptree2tree(int nnodes, int *ptree, Tree *tree)
{
    Node **nodes = tree->nodes;

    // allocate children
    for (int i=0; i<nnodes; i++) {
        nodes[i]->allocChildren(2);
        nodes[i]->name = i;
        nodes[i]->nchildren = 0;
    }

    // store parent and child pointers
    for (int i=0; i<nnodes; i++) {
        int parent = ptree[i];

        if (parent != -1) {
            Node *parentnode = nodes[parent];
            parentnode->children[parentnode->nchildren++] = nodes[i];
            nodes[i]->parent = parentnode;
        } else {
            nodes[i]->parent = NULL;
        }
    }

    // set root
    tree->root = nodes[nnodes - 1];
    assert(tree->assertTree());
}


// create a parent tree from a tree object array
void tree2ptree(Tree *tree, int *ptree)
{
    Node **nodes = tree->nodes;
    int nnodes = tree->nnodes;

    for (int i=0; i<nnodes; i++) {
        if (nodes[i]->parent)
            ptree[i] = nodes[i]->parent->name;
        else
            ptree[i] = -1;
    }
}


Tree *makeTree(int nnodes, int *ptree)
{
    Tree *tree = new Tree(nnodes);
    ptree2tree(nnodes, ptree, tree);
    return tree;
}


void deleteTree(Tree *tree)
{
    delete tree;
}

void setTreeDists(Tree *tree, double *dists)
{
    tree->setDists(dists);
}




} // extern C


} // namespace spidir
