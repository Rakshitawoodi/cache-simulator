#ifndef SIM_CACHE_H
#define SIM_CACHE_H

typedef 
struct {
   uint32_t BLOCKSIZE;
   uint32_t L1_SIZE;
   uint32_t L1_ASSOC;
   uint32_t L2_SIZE;
   uint32_t L2_ASSOC;
   uint32_t PREF_N;
   uint32_t PREF_M;
} cache_params_t;
typedef
struct {
	 uint32_t L1_BLOCK_OFFSET;
	 uint32_t L1_NO_OF_SET;
         uint32_t L1_INDEX;
         uint32_t L1_TAG;

	 uint32_t L2_BLOCK_OFFSET;
	 uint32_t L2_NO_OF_SET;
         uint32_t L2_INDEX;
         uint32_t L2_TAG;
}cache_values_t;
// Put additional data structures here as per your requirement.

#endif
