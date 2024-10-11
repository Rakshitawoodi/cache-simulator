#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include "sim.h"

// Define cache block structure
typedef struct {
    uint32_t tag;
    bool valid;
    int lru_counter;// For LRU policy
    int dirty_bit;
} cache_block_t;

// Define cache structure
typedef struct {
    cache_block_t *blocks;
    int associativity;
    int num_sets;
} cache_t;

// Function prototypes
int access_cache(cache_t *cache, uint32_t addr, uint32_t block_size, uint32_t num_sets, uint32_t assoc,char rw);
int update_cache(cache_t *cache, uint32_t addr, uint32_t block_size, uint32_t num_sets, uint32_t assoc,char rw,int update_db);
int get_lru_block(cache_t *cache, uint32_t index, int assoc);
int evict_cache(cache_t *cache, uint32_t addr, uint32_t block_size, uint32_t num_sets, uint32_t assoc, char rw,int dirty_bit);
int find_max_lru_cache(cache_t *cache, int num_sets, int assoc);
bool is_cache_full(cache_t *cache, uint32_t num_sets, uint32_t assoc);
void print_simulation_results();
// Global LRU counters (one per set)
int l1_read = 0, l1_readmiss = 0, l1_write = 0, l1_writemiss = 0;
int l2_read = 0, l2_readmiss = 0, l2_write = 0, l2_writemiss = 0, mem_traffic = 0,l1_to_l2_wb = 0,l2_to_mem_wb = 0,blocks_fetched_from_memory = 0;
int dirty_bits;
int update_dirty_bit;
//cache_t l1_cache;
//cache_t l2_cache;
cache_t l1_check;
cache_t l2_check;
uint32_t l1_assoc,l2_assoc;
int main (int argc, char *argv[]) {
    FILE *fp;            // File pointer.
    char *trace_file;     // Trace file name.
    cache_params_t params;
    cache_values_t cache_params;
    char rw;              // Request's type (read or write)
    uint32_t addr;        // Request's address
    uint32_t val[1000000];
    char rw_1[100000];
    uint32_t l1_block_offset_mask, l1_index_mask;
    uint32_t l2_block_offset_mask, l2_index_mask;
    int i = 0;
    update_dirty_bit = 0;
    if (argc != 9) {
        printf("Error: Expected 8 command-line arguments but was provided %d.\n", (argc - 1));
        exit(EXIT_FAILURE);
    }

    params.BLOCKSIZE = (uint32_t) atoi(argv[1]);
    params.L1_SIZE   = (uint32_t) atoi(argv[2]);
    params.L1_ASSOC  = (uint32_t) atoi(argv[3]);
    params.L2_SIZE   = (uint32_t) atoi(argv[4]);
    params.L2_ASSOC  = (uint32_t) atoi(argv[5]);
    params.PREF_N    = (uint32_t) atoi(argv[6]);
    params.PREF_M    = (uint32_t) atoi(argv[7]);
    trace_file       = argv[8];

    fp = fopen(trace_file, "r");
    if (fp == NULL) {
        printf("Error: Unable to open file %s\n", trace_file);
        exit(EXIT_FAILURE);
    }

    // Print simulator configuration.
    printf("===== Simulator configuration =====\n");
    printf("BLOCKSIZE:  %u\n", params.BLOCKSIZE);
    printf("L1_SIZE:    %u\n", params.L1_SIZE);
    printf("L1_ASSOC:   %u\n", params.L1_ASSOC);
    printf("L2_SIZE:    %u\n", params.L2_SIZE);
    printf("L2_ASSOC:   %u\n", params.L2_ASSOC);
    printf("PREF_N:     %u\n", params.PREF_N);
    printf("PREF_M:     %u\n", params.PREF_M);
    printf("trace_file: %s\n", trace_file);
    printf("\n");

    l1_assoc = params.L1_ASSOC;
    l2_assoc = params.L2_ASSOC;
    // L1 cache parameters
    cache_params.L1_NO_OF_SET = params.L1_SIZE / (params.BLOCKSIZE * params.L1_ASSOC);
    cache_params.L1_BLOCK_OFFSET = log2(params.BLOCKSIZE);
    cache_params.L1_INDEX = log2(cache_params.L1_NO_OF_SET);
    cache_params.L1_TAG = 32 - (cache_params.L1_BLOCK_OFFSET + cache_params.L1_INDEX);
    l1_block_offset_mask = (1 << cache_params.L1_BLOCK_OFFSET) - 1;
    l1_index_mask = ((1 << cache_params.L1_INDEX) - 1) << cache_params.L1_BLOCK_OFFSET;

    // L2 cache parameters
    cache_params.L2_NO_OF_SET = params.L2_SIZE / (params.BLOCKSIZE * params.L2_ASSOC);
    cache_params.L2_BLOCK_OFFSET = log2(params.BLOCKSIZE);
    cache_params.L2_INDEX = log2(cache_params.L2_NO_OF_SET);
    cache_params.L2_TAG = 32 - (cache_params.L2_BLOCK_OFFSET + cache_params.L2_INDEX);
    l2_block_offset_mask = (1 << cache_params.L2_BLOCK_OFFSET) - 1;
    l2_index_mask = ((1 << cache_params.L2_INDEX) - 1) << cache_params.L2_BLOCK_OFFSET;

    // Allocate caches
    cache_t l1_cache = { malloc(cache_params.L1_NO_OF_SET * params.L1_ASSOC * sizeof(cache_block_t)), params.L1_ASSOC, cache_params.L1_NO_OF_SET };
    cache_t l2_cache = { malloc(cache_params.L2_NO_OF_SET * params.L2_ASSOC * sizeof(cache_block_t)), params.L2_ASSOC, cache_params.L2_NO_OF_SET };
    l1_check = l1_cache;
    l2_check = l2_cache;

    // Initialize caches to invalid state
    for (int set = 0; set < l1_cache.num_sets; set++) {
        for (int assoc = 0; assoc < l1_cache.associativity; assoc++) {
            l1_cache.blocks[set * l1_cache.associativity + assoc].valid = false;
        }
    }
    for (int set = 0; set < l2_cache.num_sets; set++) {
        for (int assoc = 0; assoc < l2_cache.associativity; assoc++) {
            l2_cache.blocks[set * l2_cache.associativity + assoc].valid = false;
        }
    }

    // Read requests from the trace file
    while (fscanf(fp, "%c %x\n", &rw, &addr) == 2) {
        val[i] = addr;
        rw_1[i] = rw;
	printf("read tracefile %c %x\n",rw_1[i],val[i]); 
        i++;
    }
    printf("size of val = %d\n", i);

    for (int j = 0; j < i; j++) {
        if (rw_1[j] == 'r') {
	    evict_cache(&l1_check, val[j], params.BLOCKSIZE, cache_params.L1_NO_OF_SET, params.L1_ASSOC,rw_1[j],update_dirty_bit);
            if (access_cache(&l1_check, val[j], params.BLOCKSIZE, cache_params.L1_NO_OF_SET, params.L1_ASSOC,rw_1[j])) {
                l1_read++;
		printf("l1_read hit = %d\n",l1_read);
            } else {
                l1_read++;
                l1_readmiss++;
		printf("l1_read incrementing = %d,read miss %d\n",l1_read,l1_readmiss);
                if (access_cache(&l2_check, val[j], params.BLOCKSIZE, cache_params.L2_NO_OF_SET, params.L2_ASSOC,rw_1[j])) {
                    l2_read++;
		    printf("l2_read hit = %d\n",l2_read);
                    update_cache(&l1_check, val[j], params.BLOCKSIZE, cache_params.L1_NO_OF_SET, params.L1_ASSOC,rw_1[j],update_dirty_bit);
                } else {
                    l2_read++;
                    l2_readmiss++;
		    mem_traffic++;blocks_fetched_from_memory++;
		    printf("mem traffic = %d\n",mem_traffic); 
		    printf("l2_read incrementing = %d,read miss %d\n",l2_read,l2_readmiss); 
                    update_cache(&l2_check, val[j], params.BLOCKSIZE, cache_params.L2_NO_OF_SET, params.L2_ASSOC,rw_1[j],update_dirty_bit);
                    update_cache(&l1_check, val[j], params.BLOCKSIZE, cache_params.L1_NO_OF_SET, params.L1_ASSOC,rw_1[j],update_dirty_bit);
                }
            }
        } else {  // Write request
             evict_cache(&l1_check, val[j], params.BLOCKSIZE, cache_params.L1_NO_OF_SET, params.L1_ASSOC,rw_1[j],update_dirty_bit);
            if (access_cache(&l1_check, val[j], params.BLOCKSIZE, cache_params.L1_NO_OF_SET, params.L1_ASSOC,rw_1[j])) {
                l1_write++;
		printf("l1_write hit = %d\n",l1_write);
            } else {
                l1_write++;
                l1_writemiss++;
		printf("l1_write incrementing = %d,write miss %d\n",l1_write,l1_writemiss);
                if (access_cache(&l2_check, val[j], params.BLOCKSIZE, cache_params.L2_NO_OF_SET, params.L2_ASSOC,rw_1[j])) {
                    l2_read++;
		    update_dirty_bit = 1;
		    printf("l2_read hit = %d\n",l2_write);
                    update_cache(&l1_check, val[j], params.BLOCKSIZE, cache_params.L1_NO_OF_SET, params.L1_ASSOC,rw_1[j],update_dirty_bit);
		    update_dirty_bit = 0;
		    printf("update_dirty_bit made to 0 , value = %d\n",update_dirty_bit);
                } else {
                    l2_read++;
                    l2_readmiss++;
		    printf("update_dirty_bit = %d\n",update_dirty_bit);
		    printf("Because of l2 read miss after l1_write miss,entering this loop,values of l2_read = %d\n l2_readmiss = %d\n",l2_read,l2_readmiss);
		    mem_traffic++;blocks_fetched_from_memory++;
		    printf("mem traffic = %d\n",mem_traffic);
                    update_cache(&l2_check, val[j], params.BLOCKSIZE, cache_params.L2_NO_OF_SET, params.L2_ASSOC,rw_1[j],update_dirty_bit);
		    update_dirty_bit = 1;
                    update_cache(&l1_check, val[j], params.BLOCKSIZE, cache_params.L1_NO_OF_SET, params.L1_ASSOC,rw_1[j],update_dirty_bit);
		    //update_dirty_bit = 0;
                }
            }
        }
	//update_dirty_bit = 0;
    }

    // Cleanup
    free(l1_check.blocks);
    free(l2_check.blocks);
    fclose(fp);
    print_simulation_results();
    return 0;
}

// Cache access function
int access_cache(cache_t *cache, uint32_t addr, uint32_t block_size, uint32_t num_sets, uint32_t assoc,char rw) {
    uint32_t block_offset = addr & ((1 << (int)log2(block_size)) - 1);
    uint32_t index = (addr >> (int)log2(block_size)) % num_sets;
    uint32_t tag = addr >> (int)(log2(num_sets) + (int)log2(block_size));

    for (int i = 0; i < assoc; i++) {
        cache_block_t *block = &cache->blocks[index * assoc + i];
        if (block->valid && block->tag == tag) {
            block->lru_counter = 0;
	    if(rw == 'w')
		{
			block->dirty_bit = 1;
			printf("dirty_bit_incremented = %d\n",block->dirty_bit);
		}

            for (int j = 0; j < assoc; j++) {
                if (i != j) {
                    cache->blocks[index * assoc + j].lru_counter++;
                }
            }
            return true;  // Hit
        }
    }
    return false;  // Miss
}

// Update cache on miss
int update_cache(cache_t *cache, uint32_t addr, uint32_t block_size, uint32_t num_sets, uint32_t assoc,char rw,int update_db) {
    uint32_t index = (addr >> (int)log2(block_size)) % num_sets;
    uint32_t tag = addr >> (int)(log2(num_sets) + (int)log2(block_size));
    printf("update_db = %d\n",update_db);
    int lru_index = get_lru_block(cache, index, assoc);
    cache->blocks[index * assoc + lru_index].valid = true;
    cache->blocks[index * assoc + lru_index].tag = tag;
    cache->blocks[index * assoc + lru_index].lru_counter = 0;
    // Update LRU counters
    if(update_db ==1)
    {
	    cache->blocks[index * assoc + lru_index].dirty_bit = 1;
	    printf("dirty bit incremented = %d\n",cache->blocks[index * assoc + lru_index].dirty_bit);
    }

    for (int j = 0; j < assoc; j++) {
        if (lru_index != j) {
            cache->blocks[index * assoc + j].lru_counter++;
        }
    }
    return 0;
}

// Get LRU block index
int get_lru_block(cache_t *cache, uint32_t index, int assoc) {
    int lru_index = 0;
    int max_lru = -1;
    for (int i = 0; i < assoc; i++) {
        if (cache->blocks[index * assoc + i].lru_counter > max_lru) {
            max_lru = cache->blocks[index * assoc + i].lru_counter;
            lru_index = i;
        }
    }
    return lru_index;
}
int find_max_lru_cache(cache_t *cache, int num_sets, int assoc) {
    int max_lru = -1;  // Initialize the maximum LRU value to -1
    int max_lru_set = -1;  // Store the set index of the block with the maximum LRU
    int max_lru_block = -1;  // Store the block index (within the set) of the block with the maximum LRU

    // Loop over all sets in the cache
    for (int set = 0; set < num_sets; set++) {
        // Loop over all blocks in the current set
        for (int block = 0; block < assoc; block++) {
            int block_index = set * assoc + block;
            if (cache->blocks[block_index].valid && cache->blocks[block_index].lru_counter > max_lru) {
                max_lru = cache->blocks[block_index].lru_counter;  // Update the maximum LRU value
                max_lru_set = set;  // Track the set index
                max_lru_block = block;  // Track the block index
            }
        }
    }

    // If you want to return the block index (as a single value)
    if (max_lru_block != -1 && max_lru_set != -1) {
        printf("Block with max LRU is in set %d, block %d with LRU value %d\n", max_lru_set, max_lru_block, max_lru);
        return max_lru_set * assoc + max_lru_block;
    } else {
        printf("Cache is empty or invalid\n");
        return -1;  // If no valid block was found
    }
}
  int evict_cache(cache_t *cache, uint32_t addr, uint32_t block_size, uint32_t num_sets, uint32_t assoc, char rw,int dirty_bit) {
    uint32_t index = (addr >> (int)log2(block_size)) % num_sets;
    uint32_t tag = addr >> (int)(log2(num_sets) + (int)log2(block_size));
    
    // Get the LRU block within the set
    int lru_index = get_lru_block(cache, index, assoc);
    printf("update_db = %d\n", dirty_bit);
 if(is_cache_full(cache,num_sets,assoc == true))
 {
    // If the LRU block is valid and dirty, write back to l2_cache or l1_cache as needed
    if (cache->blocks[index * assoc + lru_index].valid == true && dirty_bit == 1) 
    {
        
        // Write back to L2 cache if the block is from L1
        if (cache == &l1_check) {
            l1_to_l2_wb++;  // Increment L1 to L2 writeback counter
            printf("l1_to_l2_writeback = %d\n", l1_to_l2_wb);
            l1_check.blocks[index * assoc + lru_index].valid = true;
            l1_check.blocks[index * assoc + lru_index].tag = cache->blocks[index * assoc + lru_index].tag;
            l1_check.blocks[index * assoc + lru_index].lru_counter = 0;  // Reset LRU counter
            l1_check.blocks[index * assoc + lru_index].dirty_bit = 0;

            // Update the LRU counters for the other blocks in the set
            for (int j = 0; j < assoc; j++) {
                if (lru_index != j) {
                    l1_check.blocks[index * assoc + j].lru_counter++;
                }
            }

            // Evict the LRU block from L2 cache
            int max_lru_block_index = find_max_lru_cache(&l2_check,num_sets,assoc);
            if (max_lru_block_index != -1) {
                printf("Max LRU block index in L2: %d\n", max_lru_block_index);

                // Write the value from L1 to the evicted LRU block in L2
                l2_check.blocks[max_lru_block_index].tag = cache->blocks[index * assoc + lru_index].tag;
                l2_check.blocks[max_lru_block_index].valid = true;
                l2_check.blocks[max_lru_block_index].lru_counter = 0;  // Reset LRU counter in L2
                l2_check.blocks[max_lru_block_index].dirty_bit = 1;

                for (int j,index = 0; j < l2_assoc; j++) {
                if (max_lru_block_index != j) {
                    l2_check.blocks[index * assoc + j].lru_counter++;
                }
            }
            l2_write++;
            printf("l2 write will get incremented due to l1_to_l2_writeback");
            }   
       }
	    else 
        {
            l2_to_mem_wb++;
            printf("l1_to_l2_writeback = %d\n", l1_to_l2_wb);
            l2_check.blocks[index * assoc + lru_index].valid = true;
            l2_check.blocks[index * assoc + lru_index].tag = cache->blocks[index * assoc + lru_index].tag;
            l2_check.blocks[index * assoc + lru_index].lru_counter = 0;  // Reset LRU counter
            l2_check.blocks[index * assoc + lru_index].dirty_bit = 0;

            // Update the LRU counters for the other blocks in the set
            for (int j = 0; j < assoc; j++) {
                if (lru_index != j) {
                    l2_check.blocks[index * assoc + j].lru_counter++;
                }
            }
            mem_traffic++;
            printf("mem_traffic will get incremented due to l2_to_l1_writeback");
        }
    }}
   else {
	   if(cache == &l1_check)
	   {
            printf("no dirty bit, no eviction required, updating l1_cache");
            l1_check.blocks[index * assoc + lru_index].valid = true;
            l1_check.blocks[index * assoc + lru_index].tag = cache->blocks[index * assoc + lru_index].tag;
	   }
	   else
           {printf("no dirty bit, no eviction required, updating l1_cache");
            l2_check.blocks[index * assoc + lru_index].valid = true;
            l2_check.blocks[index * assoc + lru_index].tag = cache->blocks[index * assoc + lru_index].tag;
	   }
   }
 
}

bool is_cache_full(cache_t *cache, uint32_t num_sets, uint32_t assoc) {
    // Loop through each set
    for (uint32_t set_index = 0; set_index < num_sets; set_index++) {
        bool set_full = true;  // Assume the set is full until proven otherwise

        // Loop through each block in the set
        for (uint32_t block_index = 0; block_index < assoc; block_index++) {
            // Check if any block in the set is invalid (not full)
            if (cache->blocks[set_index * assoc + block_index].valid == 0) {
                set_full = false;  // The set is not full
		                break;
	      	// No need to check further, go to the next set
            }
        }

        // If we find a set that is not full, the cache is not full
        if (!set_full) {
            return false;
	    printf("set is not full");

        }
    }

    // If all sets are full, return true
    return true;
}

void print_simulation_results() {
    // Print simulation results
    printf("===== Simulation Results =====\n");

    // L1 Statistics
    printf("L1 Reads:          %d\n",l1_read);
    printf("L1 Read Misses:    %d\n",l1_readmiss);
    printf("L1 Writes:         %d\n",l1_write);
    printf("L1 Write Misses:   %d\n",l1_writemiss);
    printf("L1 Writebacks to L2: %d\n",l1_to_l2_wb);
    printf("\n");

    // L2 Statistics
    printf("L2 Reads:          %d\n",l2_read);
    printf("L2 Read Misses:    %d\n",l2_readmiss);
    printf("L2 Writes:         %d\n",l2_write);
    printf("L2 Write Misses:   %d\n",l2_writemiss);
    printf("L2 Writebacks to Memory: %d\n",l2_to_mem_wb);
    printf("\n");

    // Overall Statistics
    printf("Blocks Fetched from Memory: %d\n", blocks_fetched_from_memory);
    printf("Memory Traffic (Writes to Memory): %d\n", mem_traffic);

    printf("================================\n");
}


    


