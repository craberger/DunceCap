#include "main.hpp"

template<class T>
struct barbell_listing: public application<T> {
  uint64_t result = 0;
  void run(std::string path){
    //create the relation (currently a column wise table)
    Relation<uint64_t,uint64_t> *R_ab = new Relation<uint64_t,uint64_t>();

//////////////////////////////////////////////////////////////////////
    //File IO (for a tsv, csv should be roughly the same)
    auto rt = debug::start_clock();
    //tsv_reader f_reader("simple.txt");
    tsv_reader f_reader(path);
    char *next = f_reader.tsv_get_first();
    R_ab->num_rows = 0;
    while(next != NULL){
      //have to code generate number of attributes here
      //maybe can accomplish with variadic templates? Might be hard.
      R_ab->get<0>().append_from_string(next);
      next = f_reader.tsv_get_next();
      R_ab->get<1>().append_from_string(next);
      next = f_reader.tsv_get_next();
      R_ab->num_rows++;
    }
    debug::stop_clock("Reading File",rt);
//////////////////////////////////////////////////////////////////////
    //Encoding
    auto et = debug::start_clock();
    //encode A
    std::vector<Column<uint64_t>> *a_attributes = new std::vector<Column<uint64_t>>();
    a_attributes->push_back(R_ab->get<0>());
    a_attributes->push_back(R_ab->get<1>());
    Encoding<uint64_t> *a_encoding = new Encoding<uint64_t>(a_attributes);
    debug::stop_clock("Encoding",et);
//////////////////////////////////////////////////////////////////////
    //Trie building

    //after all encodings are done, set up encoded relation (what is passed into the Trie)
    //You can switch the ordering here to be what you want it to be in the Trie.
    //A mapping will need to be kept in the query compiler so it is known what
    //encoding each level in the Trie should map to.
    
    auto bt = debug::start_clock();
    std::vector<size_t> *ranges_ab = new std::vector<size_t>();
    std::vector<Column<uint32_t>> *ER_ab = new std::vector<Column<uint32_t>>();
    ER_ab->push_back(a_encoding->encoded.at(0)); //perform filter, selection
    ranges_ab->push_back(a_encoding->num_distinct);
    ER_ab->push_back(a_encoding->encoded.at(1));
    ranges_ab->push_back(a_encoding->num_distinct);

    //add some sort of lambda to do selections 
    Trie<T> *TR_ab = Trie<T>::build(ER_ab,ranges_ab,[&](size_t index){
      (void) index;
      return true;
    });
    
    debug::stop_clock("Build",bt);


//////////////////////////////////////////////////////////////////////
    //Prints the relation
    //R(a,b) join T(b,c) join S(a,c)

  TrieBlock<T>* a1_block; // bag BCD
  {
      //allocate memory
      allocator::memory<uint8_t> output_buffer(R_ab->num_rows * sizeof(uint64_t) * sizeof(TrieBlock<T>));

      auto qt = debug::start_clock();

      const TrieBlock<T> H = *TR_ab->head;
      const Set<T> A = H.set;
      a1_block = new(output_buffer.get_next(0, sizeof(TrieBlock<T>))) TrieBlock<T>(H);
      a1_block->init_pointers(0,&output_buffer,TR_ab->ranges->at(0)); 

      par::reducer<size_t> num_triangles(0,[](size_t a, size_t b){
        return a + b;
      });

      A.par_foreach([&](size_t tid, uint32_t a_d){
        const Set<T> matching_b = H.get_block(a_d)->set;
        //build output B block
        TrieBlock<T>* b_block = new(output_buffer.get_next(tid, sizeof(TrieBlock<T>))) TrieBlock<T>(true);
        const size_t alloc_size = sizeof(uint64_t)*TR_ab->ranges->at(0)*2;

        Set<T> B(output_buffer.get_next(tid, alloc_size));
        b_block->set = ops::set_intersect(&B, &matching_b, &A); //intersect the B
        output_buffer.roll_back(tid, alloc_size - b_block->set.number_of_bytes);
        b_block->init_pointers(tid, &output_buffer, TR_ab->ranges->at(1)); //find out the range of level 1

        //Set a block pointer to new b block
        a1_block->set_block(a_d,a_d,b_block); 

        //Next attribute to peel off
        b_block->set.foreach_index([&](uint32_t b_i, uint32_t b_d){ // Peel off B attributes
          const TrieBlock<T>* matching_c = H.get_block(b_d);
          // Placement new!!
          TrieBlock<T>* c_block = new(output_buffer.get_next(tid, sizeof(TrieBlock<T>))) TrieBlock<T>(true);
          c_block->set = Set<T>(output_buffer.get_next(tid, alloc_size));

          const size_t count = ops::set_intersect(&c_block->set, &matching_c->set, &matching_b)->cardinality;
          num_triangles.update(tid,count);

          output_buffer.roll_back(tid, alloc_size - c_block->set.number_of_bytes);
          b_block->set_block(b_i,b_d,c_block);
          //assert(count == b_block->get_block(b_d)->set.cardinality);
        });
      });

      std::cout << num_triangles.evaluate(0) << std::endl;
      debug::stop_clock("Query",qt);
  }

  TrieBlock<T>* a2_block; // bag AEF
  {
      //allocate memory
      allocator::memory<uint8_t> output_buffer(R_ab->num_rows * sizeof(uint64_t) * sizeof(TrieBlock<T>));

      auto qt = debug::start_clock();

      const TrieBlock<T> H = *TR_ab->head;
      const Set<T> A = H.set;
      a2_block = new(output_buffer.get_next(0, sizeof(TrieBlock<T>))) TrieBlock<T>(H);
      a2_block->init_pointers(0,&output_buffer,TR_ab->ranges->at(0)); 

      par::reducer<size_t> num_triangles(0,[](size_t a, size_t b){
        return a + b;
      });

      A.par_foreach([&](size_t tid, uint32_t a_d){
        const Set<T> matching_b = H.get_block(a_d)->set;
        //build output B block
        TrieBlock<T>* b_block = new(output_buffer.get_next(tid, sizeof(TrieBlock<T>))) TrieBlock<T>(true);
        const size_t alloc_size = sizeof(uint64_t)*TR_ab->ranges->at(0)*2;

        Set<T> B(output_buffer.get_next(tid, alloc_size));
        b_block->set = ops::set_intersect(&B, &matching_b, &A); //intersect the B
        output_buffer.roll_back(tid, alloc_size - b_block->set.number_of_bytes);
        b_block->init_pointers(tid, &output_buffer, TR_ab->ranges->at(1)); //find out the range of level 1

        //Set a block pointer to new b block
        a2_block->set_block(a_d,a_d,b_block); 

        //Next attribute to peel off
        b_block->set.foreach_index([&](uint32_t b_i, uint32_t b_d){ // Peel off B attributes
          const TrieBlock<T>* matching_c = H.get_block(b_d);
          // Placement new!!
          TrieBlock<T>* c_block = new(output_buffer.get_next(tid, sizeof(TrieBlock<T>))) TrieBlock<T>(true);
          c_block->set = Set<T>(output_buffer.get_next(tid, alloc_size));

          const size_t count = ops::set_intersect(&c_block->set, &matching_c->set, &matching_b)->cardinality;
          num_triangles.update(tid,count);

          output_buffer.roll_back(tid, alloc_size - c_block->set.number_of_bytes);
          b_block->set_block(b_i,b_d,c_block);
          //assert(count == b_block->get_block(b_d)->set.cardinality);
        });
      });

      std::cout << num_triangles.evaluate(0) << std::endl;
      debug::stop_clock("Query",qt);
  }
  
  TrieBlock<T>* a3_block; // this is the edgelist R(A,B) joined with the projection of the previous 2 bags on A and B
  {      
    auto qt = debug::start_clock();
    allocator::memory<uint8_t> output_buffer(R_ab->num_rows * sizeof(uint64_t) * sizeof(TrieBlock<T>));
    const size_t alloc_size = sizeof(uint64_t)*TR_ab->ranges->at(0)*2; // why this magic number ?
    const TrieBlock<T> H = *TR_ab->head;
    const Set<T> A = H.set;

    a3_block = new(output_buffer.get_next(0, sizeof(TrieBlock<T>))) TrieBlock<T>(true);
    Set<T> new_A(output_buffer.get_next(0, alloc_size));
    a3_block->set = ops::set_intersect(&new_A, &A, &a1_block->set);
    output_buffer.roll_back(0, alloc_size - a3_block->set.number_of_bytes);
    a3_block->init_pointers(0, &output_buffer, TR_ab->ranges->at(1));
    a3_block->set.par_foreach_index([&](size_t tid, uint32_t a3_i, uint32_t a3_d) {
      const TrieBlock<T>* matching_b = H.get_block(a3_d);
      TrieBlock<T>* b_block = new(output_buffer.get_next(tid, sizeof(TrieBlock<T>))) TrieBlock<T>(true);
      Set<T> B(output_buffer.get_next(tid, alloc_size));
      b_block->set = ops::set_intersect(&B, &a2_block->set, &matching_b->set);
      output_buffer.roll_back(tid, alloc_size - b_block->set.number_of_bytes);
      a3_block->set_block(a3_i, a3_d, b_block); 
    });
    debug::stop_clock("Query",qt);
  }

  auto qt = debug::start_clock();
  // now do the topdown pass of yannakakis
  allocator::memory<uint8_t> output_buffer(100 * R_ab->num_rows * sizeof(uint64_t) * sizeof(TrieBlock<T>));
  TrieBlock<T>* a_block = a3_block;
  a_block->set.par_foreach([&](size_t tid, uint32_t a_d) {
    TrieBlock<T>* b_block = a_block->get_block(a_d);
    TrieBlock<T>* e_block = a1_block->get_block(a_d);
    if (b_block && e_block) {
      b_block->init_pointers(tid, &output_buffer, TR_ab->ranges->at(1));
      b_block->set.foreach_index([&](uint32_t b_i, uint32_t b_d) {
        TrieBlock<T>* old_c_block = a2_block->get_block(b_d);
        if (old_c_block) {
          TrieBlock<T>* c_block = new(output_buffer.get_next(tid, sizeof(TrieBlock<T>))) TrieBlock<T>(old_c_block);
          c_block->init_pointers(tid, &output_buffer, TR_ab->ranges->at(1));
          b_block->set_block(b_i, b_d, c_block);
          c_block->set.foreach_index([&](uint32_t c_i, uint32_t c_d) {
            TrieBlock<T>* old_d_block = old_c_block->get_block(c_d);
            if (old_d_block) {
              TrieBlock<T>* d_block = new(output_buffer.get_next(tid, sizeof(TrieBlock<T>))) TrieBlock<T>(old_d_block);
              c_block->set_block(c_i, c_d, d_block);
              d_block->init_pointers(tid, &output_buffer, TR_ab->ranges->at(1));
              d_block->set.foreach_index([&](uint32_t d_i, uint32_t d_d) {
                d_block->set_block(d_i, d_d, e_block);
                // e_block->init_pointers(tid, &output_buffer, TR_ab->ranges->at(1));
                // e_block->set.foreach_index([&](uint32_t e_i, uint32_t e_d) {
                //  TrieBlock<T>* f_block = e_block->get_block(e_d);
                //  f_block->set.foreach_index([&](uint32_t f_i, uint32_t f_d) {
                //    std::cout << a_d << " " << b_d << " " << c_d << " " << d_d << " " << e_d << " " << f_d << std::endl;
                //  });
                // });
              });
            }
          });
        }
      });
    }
  });
  debug::stop_clock("topdown pass of yannakakis",qt);

  std::cout << "Checking answer now" << std::endl;
  unsigned long size = 0;
  a_block->set.foreach([&](uint32_t a_d) {
      // std::cout << a_d << std::endl;
      TrieBlock<T>* b_block = a_block->get_block(a_d);
      if (b_block) {
        b_block->set.foreach([&](uint32_t b_d) {
            TrieBlock<T>* c_block = b_block->get_block(b_d);
            if (c_block) {
                c_block->set.foreach([&](uint32_t c_d) {
                TrieBlock<T>* d_block = c_block->get_block(c_d);
                  if (d_block) {
                    d_block->set.foreach([&](uint32_t d_d){
                      TrieBlock<T>* e_block = d_block->get_block(d_d);
                        if (e_block) {
                          e_block->set.foreach([&](uint32_t e_d){
                            TrieBlock<T>* f_block = e_block->get_block(e_d);
                            if (f_block) {
                              size += f_block->set.cardinality;
                              //f_block->set.foreach([&](uint32_t f_d) {
                              //  std::cout << a_d << " " << b_d << " " << c_d << " " << d_d << " " << e_d << " " << f_d << std::endl;
                              //});
                            }
                          });
                        }
                    });
                  }
              });
            }
        });
      }
    });

    std::cout << size << std::endl;
   //////////////////////////////////////////////////////////////////////
  }
};

template<class T>
application<T>* init_app(){
  return new barbell_listing<T>(); 
}
