#ifndef ARGWEAVER_MULTIARRAY_H
#define ARGWEAVER_MULTIARRAY_H

#include <stdarg.h>

namespace argweaver {


/* matrix of any dimension
 */
 struct MultiArray {

 public:
     MultiArray(int ndim, ...) : ndim(ndim), defaultVal(0) {
         va_list ap;
         va_start(ap, ndim);
         dimSize = new int[ndim];
         multipliers = new int[ndim];
         multipliers[ndim-1] = 1;
         matSize = 1;
         for (int i=0; i < ndim; i++) {
             dimSize[i] = va_arg(ap, int);
             matSize *= dimSize[i];
         }
         va_end(ap);
         if (ndim >= 2)
             multipliers[ndim-2] = dimSize[ndim-1];
         for (int i=ndim-3; i >= 0; i--)
             multipliers[i] = multipliers[i+1]*dimSize[i+1];
         mat = new double[matSize];
     }
     ~MultiArray() {
         delete [] dimSize;
         delete [] multipliers;
         delete [] mat;
     }
     void setDefault(double val) {
         defaultVal = val;
     }

     void set(int *idxArr, double val) {
         int idx=0;
         for (int i=0; i < ndim; i++)
             idx += idxArr[i] * multipliers[i];
         mat[idx] = val;
     }
     void set(double val, ...) {
         int idx=0;
         va_list ap;
         va_start(ap, val);
         for (int i=0; i < ndim; i++)
             idx += va_arg(ap, int) * multipliers[i];
         va_end(ap);
         mat[idx] = val;
     }
     void set_all(double val) {
         for (int i=0; i < matSize; i++)
             mat[i] = val;
     }

     void set(double val, int pos0, int pos1) {
#ifdef DEBUG
         assert(ndim == 2);
#endif
         int idx = pos0 * multipliers[0] + pos1;
#ifdef DEBUG
         assert(idx >= 0 && idx < matSize);
#endif
         mat[idx] = val;
     }

     void set(double val, int pos0, int pos1, int pos2) {
#ifdef DEBUG
         assert(ndim == 3);
#endif
         int idx = pos0 * multipliers[0] + pos1 * multipliers[1] + pos2;
#ifdef DEBUG
         assert(idx >= 0 && idx < matSize);
#endif
         mat[idx] = val;
     }

     void addVal(double val, int pos0, int pos1) {
         int idx=0;
#ifdef DEBUG
         assert(ndim == 2);
#endif
         idx = pos0 * multipliers[0] + pos1;
#ifdef DEBUG
         assert(idx >= 0 && idx < matSize);
#endif
         mat[idx] += val;
     }
     void addVal(double val, int pos0, int pos1, int pos2) {
         int idx=0;
#ifdef DEBUG
         assert(ndim == 3);
#endif
         idx = pos0 * multipliers[0] + pos1 * multipliers[1] + pos2;
#ifdef DEBUG
         assert(idx >= 0 && idx < matSize);
#endif
         mat[idx] += val;
     }

     void logAddVal(double lnVal, int pos0, int pos1) {
         int idx=0;
#ifdef DEBUG
         assert(ndim == 2);
#endif
         idx = pos0 * multipliers[0] + pos1;
#ifdef DEBUG
         assert(idx >= 0 && idx < matSize);
#endif
         mat[idx] = logadd(mat[idx], lnVal);
     }
     void logAddVal(double lnVal, int pos0, int pos1, int pos2) {
         int idx=0;
#ifdef DEBUG
         assert(ndim == 3);
#endif
         idx = pos0 * multipliers[0] + pos1 * multipliers[1] + pos2;
#ifdef DEBUG
         assert(idx >= 0 && idx < matSize);
#endif
         mat[idx] = logadd(mat[idx], lnVal);
     }

     double get(int pos0, int pos1) {
#ifdef DEBUG
         assert(ndim ==2);
#endif
         if (pos0 < 0 || pos1 < 0) return defaultVal;
         int idx = pos0 * multipliers[0] + pos1;
         return mat[idx];
     }

     double get(int pos0, int pos1, int pos2) {
#ifdef DEBUG
         assert(ndim == 3);
#endif
         if (pos0 < 0 || pos1 < 0 || pos2 < 0) return defaultVal;
         //         printf("get %i %i %i %i %i %i\n", pos0, pos1, pos2, multipliers[0], multipliers[1], multipliers[2]);
         int idx = pos0*multipliers[0] + pos1 * multipliers[1] + pos2;
         return mat[idx];
     }

     double *getRow(int pos0) {
         assert(ndim == 2);
         return &mat[pos0*multipliers[0]];
     }

     double *getRow(int pos0, int pos1) {
         assert(ndim == 3);
         return &mat[pos0*multipliers[0] + pos1*multipliers[1]];
     }


     // this works for general case but makes debugging difficult; can't
     // call from gdb
     // should have the same number of arguments as the dimension of the matrix
     // returns mat[pos1][pos2][...]
     // if any of the positions are < 0, returns 0.0
     double get(int pos1, ...) {
         va_list ap;
         int idx = pos1 * multipliers[0];
         va_start(ap, pos1);
         for (int i=1; i < ndim; i++) {
             int pos = va_arg(ap, int);
             if (pos < 0) return 0.0;
             idx += pos * multipliers[i];
         }
         va_end(ap);
         return mat[idx];
     }

     int ndim;
     int *dimSize;
     int *multipliers;
     double *mat;
     int matSize;
     double defaultVal;
 };
}
#endif
