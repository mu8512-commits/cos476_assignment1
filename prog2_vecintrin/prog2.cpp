#include "prog2_intrin.h"

// computes the absolute value of all elements in the input array
// values, stores result in output
void absSerial(float* values, float* output, int N) {
  for (int i = 0; i < N; i++) {
    float x = values[i];
    if (x < 0) {
      output[i] = -x;
    } else {
      output[i] = x;
    }
  }
}

// implementation of absSerial() above, but it is vectorized using the fake
// intrinsics
void absVector(float* values, float* output, int N) {
  __prog2_vec_float x;
  __prog2_vec_float result;
  __prog2_vec_float zero = _prog2_vset_float(0.f);
  __prog2_mask maskAll, maskIsNegative, maskIsNotNegative;

  //  Note: Take a careful look at this loop indexing.  This example
  //  code is not guaranteed to work when (N % VECTOR_WIDTH) != 0.
  //  Why is that the case?
  for (int i = 0; i < N; i += VECTOR_WIDTH) {
    // All ones
    maskAll = _prog2_init_ones();

    // All zeros
    maskIsNegative = _prog2_init_ones(0);

    // Load vector of values from contiguous memory addresses
    _prog2_vload_float(x, values + i, maskAll);  // x = values[i];

    // Set mask according to predicate
    _prog2_vlt_float(maskIsNegative, x, zero, maskAll);  // if (x < 0) {

    // Execute instruction using mask ("if" clause)
    _prog2_vsub_float(result, zero, x, maskIsNegative);  //   output[i] = -x;

    // Inverse maskIsNegative to generate "else" mask
    maskIsNotNegative = _prog2_mask_not(maskIsNegative);  // } else {

    // Execute instruction ("else" clause)
    _prog2_vload_float(result, values + i,
                       maskIsNotNegative);  //   output[i] = x; }

    // Write results back to memory
    _prog2_vstore_float(output + i, result, maskAll);
  }
}

// accepts an array of values and an array of exponents
//
// For each element, compute values[i]^exponents[i] and clamp value to
// 9.999.  Store result in output.
void clampedExpSerial(float* values, int* exponents, float* output, int N) {
  for (int i = 0; i < N; i++) {
    float x = values[i];
    int y = exponents[i];
    if (y == 0) {
      output[i] = 1.f;
    } else {
      float result = x;
      int count = y - 1;
      while (count > 0) {
        result *= x;
        count--;
      }
      if (result > 9.999999f) {
        result = 9.999999f;
      }
      output[i] = result;
    }
  }
}


void clampedExpVector(float* values, int* exponents, float* output, int N) {

    const float kClamp = 9.999999f;

    __prog2_vec_float clampVec = _prog2_vset_float(kClamp);
    __prog2_vec_int   zeroI    = _prog2_vset_int(0);
    __prog2_vec_int   oneI     = _prog2_vset_int(1);

    for (int i = 0; i < N; i += VECTOR_WIDTH) {

        // ----- Tail mask -----
        int remaining = N - i;
        __prog2_mask maskAll = (remaining >= VECTOR_WIDTH)
            ? _prog2_init_ones(VECTOR_WIDTH)
            : _prog2_init_ones(remaining);

        // ----- Load x and y -----
        __prog2_vec_float x;
        __prog2_vec_int   y;
        _prog2_vload_float(x, values + i, maskAll);
        _prog2_vload_int(y, exponents + i, maskAll);

        // ----- Identify lanes with y == 0 -----
        __prog2_mask y_eq_0 = _prog2_init_ones(0);
        _prog2_veq_int(y_eq_0, y, zeroI, maskAll);

        // Start output with 1.0 everywhere (then overwrite nonzero lanes)
        __prog2_vec_float outVec = _prog2_vset_float(1.0f);

        // nonzero mask = maskAll & !(y == 0)
        __prog2_mask not_y_eq_0 = _prog2_mask_not(y_eq_0);
        __prog2_mask nonzero = _prog2_mask_and(maskAll, not_y_eq_0);

        // For y != 0: result = x, count = y - 1
        __prog2_vec_float result = outVec;  // will be overwritten under nonzero
        _prog2_vmove_float(result, x, nonzero);

        __prog2_vec_int count = y;          // count = y - 1 under nonzero
        _prog2_vsub_int(count, y, oneI, nonzero);

        // active = (count > 0) & nonzero
        // store result in count_gt_0
        __prog2_mask count_gt_0 = _prog2_init_ones(0);
        _prog2_vgt_int(count_gt_0, count, zeroI, nonzero);
        __prog2_mask active = _prog2_mask_and(nonzero, count_gt_0);

        // while (count > 0) { result *= x; count--; }
        while (_prog2_cntbits(active) > 0) {
            // store multiplcation in result 
            _prog2_vmult_float(result, result, x, active);

            // subtract count in active lanes by one 
            _prog2_vsub_int(count, count, oneI, active);
            
            // get the count_gt_0 vector and combine it with active mask 
            _prog2_vgt_int(count_gt_0, count, zeroI, nonzero);
            active = _prog2_mask_and(nonzero, count_gt_0);
        }

        // clamp for nonzero lanes: if (result > clamp) result = clamp
        __prog2_mask over = _prog2_init_ones(0);
        _prog2_vgt_float(over, result, clampVec, nonzero);
        _prog2_vmove_float(result, clampVec, over);

        // write results into outVec for nonzero lanes
        _prog2_vmove_float(outVec, result, nonzero);

        // store
        _prog2_vstore_float(output + i, outVec, maskAll);
    }
}

// returns the sum of all elements in values
float arraySumSerial(float* values, int N) {
  float sum = 0;
  for (int i = 0; i < N; i++) {
    sum += values[i];
  }

  return sum;
}

// returns the sum of all elements in values
// You can assume N is a multiple of VECTOR_WIDTH
float arraySumVector(float* values, int N) {
    // N is a multiple of VECTOR_WIDTH, so all lanes always valid.
    __prog2_mask maskAll = _prog2_init_ones(VECTOR_WIDTH);

    __prog2_vec_float vecSum = _prog2_vset_float(0.0f);

    // 1) Vector accumulate
    // this is absed on the assumption that VECTOR_WIDTH is a factor of N 
    // all of these elements are VECTOR_WIDTH vectors
    for (int i = 0; i < N; i += VECTOR_WIDTH) {
        __prog2_vec_float x;
        //always loading VECTOR_WIDTH length 
        _prog2_vload_float(x, values + i, maskAll);
        _prog2_vadd_float(vecSum, vecSum, x, maskAll);
    }

    // 2) Store vector lanes to scalar temp array
    float tmp[VECTOR_WIDTH];
    
    // Initialize to avoid any chance of uninitialized read (paranoid-safe)
    for (int k = 0; k < VECTOR_WIDTH; k++) tmp[k] = 0.0f;
    _prog2_vstore_float(tmp, vecSum, maskAll);

    // 3) Scalar sum across lanes
    float sum = 0.0f;
    for (int k = 0; k < VECTOR_WIDTH; k++) {
        sum += tmp[k];
    }

    return sum;
}