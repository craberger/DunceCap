#ifndef UINTEGER_H
#define UINTEGER_H
/*

THIS CLASS IMPLEMENTS THE FUNCTIONS ASSOCIATED WITH AN UNCOMPRESSED SET LAYOUT.
QUITE SIMPLE THE LAYOUT JUST CONTAINS UNSIGNED INTEGERS IN THE SET.

*/

#include "utils/utils.hpp"

class uinteger{
  public:
    static type::layout get_type();
    static std::tuple<size_t,type::layout> build(uint8_t *r_in, const uint32_t *data, const size_t length);

    static long find(uint32_t key, const uint8_t *data_in, const size_t number_of_bytes, const type::layout t);

    template<typename F>
    static void foreach(
        F f,
        const uint8_t *data_in,
        const size_t cardinality,
        const size_t number_of_bytes,
        const type::layout t);

    template<typename F>
    static void foreach_index(
        F f,
        const uint8_t *data_in,
        const size_t cardinality,
        const size_t number_of_bytes,
        const type::layout t);

    template<typename F>
    static void foreach_until(
        F f,
        const uint8_t *data_in,
        const size_t cardinality,
        const size_t number_of_bytes,
        const type::layout t);

    template<typename F>
    static size_t par_foreach(
      F f,
      const uint8_t *data_in,
      const size_t cardinality,
      const size_t number_of_bytes,
      const type::layout t);

    template<typename F>
    static size_t par_foreach_index(
      F f,
      const uint8_t *data_in,
      const size_t cardinality,
      const size_t number_of_bytes,
      const type::layout t);

    static std::tuple<size_t,size_t,type::layout> intersect(uint8_t *C_in, const uint8_t *A_in, const uint8_t *B_in, const size_t A_cardinality, const size_t B_cardinality, const size_t A_num_bytes, const size_t B_num_bytes, const type::layout a_t, const type::layout b_t);
};

inline type::layout uinteger::get_type(){
  return type::UINTEGER;
}
//Copies data from input array of ints to our set data r_in
inline std::tuple<size_t,type::layout> uinteger::build(uint8_t *r_in, const uint32_t *data, const size_t length){
  uint32_t *r = (uint32_t*) r_in;
  std::copy(data,data+length,r);
  return std::make_pair(length*sizeof(uint32_t),type::UINTEGER);
}

//Iterates over set applying a lambda.
template<typename F>
inline void uinteger::foreach(
    F f,
    const uint8_t *data_in,
    const size_t cardinality,
    const size_t number_of_bytes,
    const type::layout t) {
 (void) number_of_bytes; (void) t;

 uint32_t *data = (uint32_t*) data_in;
 for(size_t i=0; i<cardinality;i++){
  f(data[i]);
 }
}

//Iterates over set applying a lambda.
template<typename F>
inline void uinteger::foreach_index(
    F f,
    const uint8_t *data_in,
    const size_t cardinality,
    const size_t number_of_bytes,
    const type::layout t) {
 (void) number_of_bytes; (void) t; (void) data_in;
  const uint32_t * const data = (uint32_t*) data_in;
  for(size_t i=0; i<cardinality;i++){
    f(i,data[i]);
  }
}

//Iterates over set applying a lambda.
template<typename F>
inline void uinteger::foreach_until(
    F f,
    const uint8_t *data_in,
    const size_t cardinality,
    const size_t number_of_bytes,
    const type::layout t) {
 (void) number_of_bytes; (void) t;

  const uint32_t * const data = (uint32_t*) data_in;
  for(size_t i=0; i<cardinality;i++){
    if(f(data[i]))
      break;
  }
}

// Iterates over set applying a lambda in parallel.
template<typename F>
inline size_t uinteger::par_foreach(
      F f,
      const uint8_t *data_in,
      const size_t cardinality,
      const size_t number_of_bytes,
      const type::layout t) {
   (void) number_of_bytes; (void) t;

   const uint32_t * const data = (uint32_t*) data_in;
   return par::for_range(0, cardinality, 128,
     [&](size_t tid, size_t i) {
        f(tid, data[i]);
     });
}

// Iterates over set applying a lambda in parallel.
template<typename F>
inline size_t uinteger::par_foreach_index(
      F f,
      const uint8_t *data_in,
      const size_t cardinality,
      const size_t number_of_bytes,
      const type::layout t) {
   (void) number_of_bytes; (void) t; (void) data_in;

   const uint32_t * const data = (uint32_t*) data_in;
   return par::for_range(0, cardinality, 128,
     [&](size_t tid, size_t i) {
        f(tid, (uint32_t)i, data[i]);
     });
}

inline long uinteger::find(uint32_t key, 
  const uint8_t *data_in, 
  const size_t number_of_bytes,
  const type::layout t){
  (void) t;

  if(number_of_bytes == 0) return -1;
  else return binary_search((uint32_t*)data_in,0,(number_of_bytes/sizeof(uint32_t))-1,key,[&](uint32_t d){return d;});
}

#endif
