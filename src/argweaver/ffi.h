#ifndef ARGWEAVER_FFI_H
#define ARGWEAVER_FFI_H
#include "argweaver/compress.h"
#include "local_tree.h"
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using namespace std;

inline vector<string> empty_seqnames() { return vector<string>(); }
inline void vector_push_string(std::vector<std::string> &vec,
                               std::string &element) {
  vec.emplace_back(std::move(element));
}
inline vector<double> empty_times() { return vector<double>(); }
inline void vector_push_double(std::vector<double> &vec, double element) {
  vec.emplace_back(element);
}

inline unique_ptr<argweaver::LocalTrees>
read_local_trees(string filename, vector<double> &times) {
  argweaver::CompressStream instream(filename.c_str(), "r");
  if (instream.stream == NULL) {
    argweaver::printError("cannot read file '%s'\n", filename.c_str());
    return NULL;
  }

  vector<string> seqnames;
  unique_ptr<argweaver::LocalTrees> trees(new argweaver::LocalTrees());
  bool result = argweaver::read_local_trees(
      instream.stream, times.data(), times.size(), trees.get(), seqnames);

  instream.close();
  return trees;
}
inline void write_local_trees_as_bed(string filename,
                                     unique_ptr<argweaver::LocalTrees> trees,
                                     vector<string> seqnames,
                                     unique_ptr<argweaver::ArgModel> model,
                                     int sample) {
  FILE *outfile = fopen(filename.c_str(), "w");
  write_local_trees_as_bed(outfile, trees.get(), seqnames, model.get(), sample);
  fclose(outfile);
}
#endif // ARGWEAVER_FFI_H
