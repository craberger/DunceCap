#ifndef _TRIEBLOCK_H_
#define _TRIEBLOCK_H_

#include "set/ops.hpp"

template<class T, class R>
struct TrieBlock{
  Set<T> set;
  TrieBlock<T,R>** next_level;
  R* values;
  bool is_sparse;

  TrieBlock(TrieBlock *init){
    set = init->set;
    next_level = init->next_level;
    is_sparse = init->is_sparse;
    values = init->values;
  }
  TrieBlock(Set<T> setIn){
    set = setIn;
    values = NULL;
  }
  TrieBlock(){
    values = NULL;
  }
  TrieBlock(bool sparse){
    is_sparse = sparse;
    values = NULL;
  }

  /*
  * Write to a binary file. The prev_index and prev_data allow us to set the pointers
  * from the previous level to the current level. We write out these first, followed
  * by the set data then next the pointers to the next level.
  */
  void to_binary(std::ofstream* outfile, const uint32_t prev_index, const uint32_t prev_data){
    outfile->write((char *)&prev_index, sizeof(prev_index));
    outfile->write((char *)&prev_data, sizeof(prev_data));
    outfile->write((char *)&is_sparse, sizeof(is_sparse));
    set.to_binary(outfile);
  }

  static std::tuple<TrieBlock<T,R>*,uint32_t,uint32_t> from_binary(std::ifstream* infile, 
      allocator::memory<uint8_t> *allocator_in, 
      const size_t tid){
    
    uint32_t prev_index; uint32_t prev_data;

    infile->read((char *)&prev_index, sizeof(prev_index));
    infile->read((char *)&prev_data, sizeof(prev_data));

    //pull from an allocator and create a trieblock
    TrieBlock<T,R>* output_block = new (
              allocator_in->get_next(tid, sizeof(TrieBlock<T, R>)))
              TrieBlock<T, R>();

    infile->read((char *)&output_block->is_sparse, sizeof(output_block->is_sparse));
    
    output_block->set = *Set<T>::from_binary(infile,allocator_in,tid);
    /*
    std::cout << "reading binary: " << tid << std::endl;
    output_block->set.foreach([&](uint32_t data){
      std::cout << data << std::endl;
    });
    */

    return std::tuple<TrieBlock<T,R>*,uint32_t,uint32_t>(output_block,prev_index,prev_data);
  }

  //refactor this code

  void init_pointers(const size_t tid, allocator::memory<uint8_t> *allocator_in){
    is_sparse = (set.range == 0) ? ((double)set.cardinality/(double)set.range) > (1.0/256.0) : true;
    if(!is_sparse){
      next_level = (TrieBlock<T,R>**)allocator_in->get_next(tid, sizeof(TrieBlock<T,R>*)*(set.range+1) );
    } else{
      next_level = (TrieBlock<T,R>**)allocator_in->get_next(tid, sizeof(TrieBlock<T,R>*)*set.cardinality);
    }
  }

  void alloc_data(size_t tid, allocator::memory<uint8_t> *allocator_in, const size_t cardinality, const size_t range){
    if(!is_sparse){
      values = (R*)allocator_in->get_next(tid, sizeof(R)*(range+1));
    } else{
      values = (R*)allocator_in->get_next(tid, sizeof(R)*cardinality);
    }
  }

  void init_data(const size_t tid, allocator::memory<uint8_t> *allocator_in, const size_t cardinality, const size_t range, const R value){
    if(!is_sparse){
      values = (R*)allocator_in->get_next(tid, sizeof(R)*(range+1));
      std::fill(values,values+range+1,value);
    } else{
      values = (R*)allocator_in->get_next(tid, sizeof(R)*(cardinality));
      std::fill(values,values+cardinality,value);
    }
  }

  void set_data(uint32_t index, uint32_t data, R value){
    if(!is_sparse){
      (void) index;
      values[data] = value;
    } else{
      (void) data;
      values[index] = value;
    }
  }

  R get_data(uint32_t data){
    if(!is_sparse){
      return values[data];
    } else{
      //something like get the index from the set then move forward.
      const long index = set.find(data);
      if(index != -1)
        return values[data];
    }
    //FIXME
    return (R)0;
  }

  R get_data(uint32_t index, uint32_t data){
    if(!is_sparse){
      (void) index;
      return values[data];
    } else{
      (void) data;
      return values[index];
    }
  }

  void set_block(uint32_t index, uint32_t data, TrieBlock<T,R> *block){
    if(!is_sparse){
      (void) index;
      next_level[data] = block;
    } else{
      (void) data;
      next_level[index] = block;
    }
  }

  inline std::tuple<size_t,TrieBlock<T,R>*> get_block_forward(size_t index, uint32_t data) const{
    if(!is_sparse){
      return std::make_tuple(index,next_level[data]);
    } else{
      //something like get the index from the set then move forward.
      auto tup = set.find(index,data);
      size_t find_index = std::get<0>(tup);
      if(std::get<1>(tup))
        return std::make_tuple(find_index,next_level[find_index]);
      else
        return std::make_tuple(find_index,(TrieBlock<T,R>*)NULL);
    }
  }
  inline TrieBlock<T,R>* get_block(uint32_t data) const{
    TrieBlock<T,R>* result = NULL;
    if(!is_sparse){
      result = next_level[data];
    } else{
      //something like get the index from the set then move forward.
      const long index = set.find(data);
      if(index != -1)
        result = next_level[index];
    }
    return result;
  }
  inline TrieBlock<T,R>* get_block(uint32_t index, uint32_t data) const{
    if(!is_sparse){
      return next_level[data];
    } else{
      return next_level[index];
    }
    return NULL;
  }

};

template<class T,class R>
struct TrieBlockIterator{
  size_t pointer_index;
  TrieBlock<T,R>* trie_block;

  TrieBlockIterator(TrieBlock<T,R> *init){
    pointer_index = 0;
    trie_block = init;
  }

  TrieBlockIterator<T,R> get_block(uint32_t data) {
    auto tup = trie_block->get_block_forward(0,data);
    pointer_index = std::get<0>(tup);
    return TrieBlockIterator(std::get<1>(tup));
  }
};

#endif
