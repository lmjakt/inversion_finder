#include <R.h>
#include <Rinternals.h>
#include <map>

// defines one function that does a running kmer-count
// with both forward and reverse complement sequences
// reverse complement sequences decrease k-mer counts
// forward reads increases. Hence counts go up, and
// and then down.. 

// Sequences will be represented in two bits as unsigned integers
// i.e. we will only count 16-mers.
extern "C" {
  // seqs_r, is simply a character vector of sequence types
  SEXP find_inversions( SEXP seqs_r){
    if(TYPEOF( seqs_r ) != STRSXP)
      error("seqs_r should be a character vector");
    int seq_n = length(seqs_r);
    if( seq_n < 1 )
      error("seqs_r must have positive length");
    // Return a list containing the running sum for each one of these..
    SEXP ret_data = PROTECT(allocVector( VECSXP, seq_n ));
    for(int i=0; i < seq_n; ++i){
      SEXP seq_r = STRING_ELT(seqs_r, i);
      int length = length(seq_r);
      const char *seq = CHAR(seq_r);
      SET_VECTOR_ELT( ret_data, i, allocVector(INTSXP, length) );
      int *sum_vector = INTEGER( VECTOR_ELT(ret_data, i));
      memset( (void*)sum_vector, 0, sizeof(int) * length );
      // note that counts can be negative!
      std::map<unsigned int, int> k_counts;
      std::map<unsigned int, int>::iterator map_it;
      int running_sum = 0;  // a running sum updated as we increase values
      unsigned int forward = 0;
      unsigned int reverse = 0;
      unsigned int word_length = sizeof(unsigned int) * 8;
      // first encode the initial k-mer;
      // This uses bitwise operations to avoid conditionals and memory lookups
      // Note that this will only do the correct thing for A,a,C,c,G,g,T,t
      // Ambiguity symbols will do undefined but not dangrous things.
      // The ASCII for A, C, G, T are
      // A: 0x41     0100 0001  -> 00
      // C: 0x43     0100 0011  -> 01
      // T: 0x54     0101 0100  -> 10
      // G: 0x47     0100 0111  -> 11
      // The forward operations do the following (where c is seq[j])
      // c & 0x07        sets all higher bits to 0
      //   >> 1           shifts remaining bits giving 00, 01, 10, 11.           
      // 
      // To convert to complement we see that we have the complementary
      // pairs: A/T : 0/2 and C/G : 1/3
      // And to note that +2 % 4 will do the correct thing:
      // (0 + 2) % 4  -> 2 (A -> T)
      // (1 + 2) % 4  -> 3 (C -> G)
      // (2 + 2) % 4  -> 0 (T -> A)
      // (3 + 2) % 4  -> 1 (G -> C)
      // Should really make an inline function for this:
      size_t j = 0;
      while(j < word_length / 2 && seq[j]){
	forward <<= 2;
	forward |= ( ((unsigned int)seq[j] & 0x07) >> 1);
	reverse >>= 2;
	reverse |= (((( forward & 3 ) + 2) % 4) << (word_length-2) );
	++j;
      }
      if(forward != reverse)
	k_counts[ forward ] = 1;
      while( seq[j] ){
	// c++ is almost as easy as perl!
	map_it = k_counts.find(reverse);
	k_counts[forward]++;
	if(map_it == k_counts.end() || map_it->second == 0){
	  // k_counts[ forward ]++;
	  running_sum++;
	}else{
	  k_counts[ forward ]--;
	  map_it->second--;
	  running_sum--;
	}
	sum_vector[j] = running_sum;
	forward <<= 2;
	forward |= ( ((unsigned int)seq[j] & 0x07) >> 1);
	reverse >>= 2;
	reverse |= (((( forward & 3 ) + 2) % 4) << (word_length-2) );
	++j;
      }
    }
    UNPROTECT(1);
    return( ret_data );
  }
}
