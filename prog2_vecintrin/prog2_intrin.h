// Define vector unit width here
#define VECTOR_WIDTH 5

#ifndef prog2_INTRIN_H_
#define prog2_INTRIN_H_

#include <cstdlib>
#include <cmath>
#include "logger.h"

//*******************
//* Type Definition *
//*******************

extern Logger prog2_Logger;

template <typename T>
struct __prog2_vec {
  T value[VECTOR_WIDTH];
};

// Declare a mask with __prog2_mask
struct __prog2_mask : __prog2_vec<bool> {};

// Declare a floating point vector register with __prog2_vec_float
#define __prog2_vec_float __prog2_vec<float>

// Declare an integer vector register with __prog2_vec_int
#define __prog2_vec_int   __prog2_vec<int>

//***********************
//* Function Definition *
//***********************

// Return a mask initialized to 1 in the first N lanes and 0 in the others
__prog2_mask _prog2_init_ones(int first = VECTOR_WIDTH);

// Return the inverse of maska
__prog2_mask _prog2_mask_not(__prog2_mask &maska);

// Return (maska | maskb)
__prog2_mask _prog2_mask_or(__prog2_mask &maska, __prog2_mask &maskb);

// Return (maska & maskb)
__prog2_mask _prog2_mask_and(__prog2_mask &maska, __prog2_mask &maskb);

// Count the number of 1s in maska
int _prog2_cntbits(__prog2_mask &maska);

// Set register to value if vector lane is active
//  otherwise keep the old value
void _prog2_vset_float(__prog2_vec_float &vecResult, float value, __prog2_mask &mask);
void _prog2_vset_int(__prog2_vec_int &vecResult, int value, __prog2_mask &mask);
// For user's convenience, returns a vector register with all lanes initialized to value
__prog2_vec_float _prog2_vset_float(float value);
__prog2_vec_int _prog2_vset_int(int value);

// Copy values from vector register src to vector register dest if vector lane active
// otherwise keep the old value
void _prog2_vmove_float(__prog2_vec_float &dest, __prog2_vec_float &src, __prog2_mask &mask);
void _prog2_vmove_int(__prog2_vec_int &dest, __prog2_vec_int &src, __prog2_mask &mask);

// Load values from array src to vector register dest if vector lane active
//  otherwise keep the old value
void _prog2_vload_float(__prog2_vec_float &dest, float* src, __prog2_mask &mask);
void _prog2_vload_int(__prog2_vec_int &dest, int* src, __prog2_mask &mask);

// Store values from vector register src to array dest if vector lane active
//  otherwise keep the old value
void _prog2_vstore_float(float* dest, __prog2_vec_float &src, __prog2_mask &mask);
void _prog2_vstore_int(int* dest, __prog2_vec_int &src, __prog2_mask &mask);

// Return calculation of (veca + vecb) if vector lane active
//  otherwise keep the old value
void _prog2_vadd_float(__prog2_vec_float &vecResult, __prog2_vec_float &veca, __prog2_vec_float &vecb, __prog2_mask &mask);
void _prog2_vadd_int(__prog2_vec_int &vecResult, __prog2_vec_int &veca, __prog2_vec_int &vecb, __prog2_mask &mask);

// Return calculation of (veca - vecb) if vector lane active
//  otherwise keep the old value
void _prog2_vsub_float(__prog2_vec_float &vecResult, __prog2_vec_float &veca, __prog2_vec_float &vecb, __prog2_mask &mask);
void _prog2_vsub_int(__prog2_vec_int &vecResult, __prog2_vec_int &veca, __prog2_vec_int &vecb, __prog2_mask &mask);

// Return calculation of (veca * vecb) if vector lane active
//  otherwise keep the old value
void _prog2_vmult_float(__prog2_vec_float &vecResult, __prog2_vec_float &veca, __prog2_vec_float &vecb, __prog2_mask &mask);
void _prog2_vmult_int(__prog2_vec_int &vecResult, __prog2_vec_int &veca, __prog2_vec_int &vecb, __prog2_mask &mask);

// Return calculation of (veca / vecb) if vector lane active
//  otherwise keep the old value
void _prog2_vdiv_float(__prog2_vec_float &vecResult, __prog2_vec_float &veca, __prog2_vec_float &vecb, __prog2_mask &mask);
void _prog2_vdiv_int(__prog2_vec_int &vecResult, __prog2_vec_int &veca, __prog2_vec_int &vecb, __prog2_mask &mask);


// Return calculation of absolute value abs(veca) if vector lane active
//  otherwise keep the old value
void _prog2_vabs_float(__prog2_vec_float &vecResult, __prog2_vec_float &veca, __prog2_mask &mask);
void _prog2_vabs_int(__prog2_vec_int &vecResult, __prog2_vec_int &veca, __prog2_mask &mask);

// Return a mask of (veca > vecb) if vector lane active
//  otherwise keep the old value
void _prog2_vgt_float(__prog2_mask &vecResult, __prog2_vec_float &veca, __prog2_vec_float &vecb, __prog2_mask &mask);
void _prog2_vgt_int(__prog2_mask &vecResult, __prog2_vec_int &veca, __prog2_vec_int &vecb, __prog2_mask &mask);

// Return a mask of (veca < vecb) if vector lane active
//  otherwise keep the old value
void _prog2_vlt_float(__prog2_mask &vecResult, __prog2_vec_float &veca, __prog2_vec_float &vecb, __prog2_mask &mask);
void _prog2_vlt_int(__prog2_mask &vecResult, __prog2_vec_int &veca, __prog2_vec_int &vecb, __prog2_mask &mask);

// Return a mask of (veca == vecb) if vector lane active
//  otherwise keep the old value
void _prog2_veq_float(__prog2_mask &vecResult, __prog2_vec_float &veca, __prog2_vec_float &vecb, __prog2_mask &mask);
void _prog2_veq_int(__prog2_mask &vecResult, __prog2_vec_int &veca, __prog2_vec_int &vecb, __prog2_mask &mask);

// Adds up adjacent pairs of elements, so
//  [0 1 2 3] -> [0+1 0+1 2+3 2+3]
void _prog2_hadd_float(__prog2_vec_float &vecResult, __prog2_vec_float &vec);

// Performs an even-odd interleaving where all even-indexed elements move to front half
//  of the array and odd-indexed to the back half, so
//  [0 1 2 3 4 5 6 7] -> [0 2 4 6 1 3 5 7]
void _prog2_interleave_float(__prog2_vec_float &vecResult, __prog2_vec_float &vec);

// Add a customized log to help debugging
void addUserLog(const char * logStr);

#endif
