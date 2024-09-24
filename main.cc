/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

// @sect3{Include files}

#include <boost/lambda/lambda.hpp>
#include <examples.h>



// The final step in importing deal.II is this: All deal.II functions and
// classes are in a namespace <code>dealii</code>, to make sure they don't
// clash with symbols from other libraries you may want to use in conjunction
// with deal.II. One could use these functions and classes by prefixing every
// use of these names by <code>dealii::</code>, but that would quickly become
// cumbersome and annoying. Rather, we simply import the entire deal.II
// namespace for general use:
using namespace dealii;
using namespace dealt;


int main(int argc, char* argv[])
{
  try {
    MultithreadInfo::set_thread_limit();
    MultithreadInfo::initialize_multithreading();
  } catch(std::exception &exc) {
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl;
  }

  int ref = 10, order = 0; 
  std::string num = "disc";
  if (argc > 1){ 
    for (int i = 1; i < argc; i+=2){
      if (std::string(argv[i]) == "-n" || std::string(argv[i]) == "--example"){
        num = argv[i+1]; 
//        char* p; 
//        long conv = strtol(argv[i+1], &p, 10);
//        if(*p != '\0' || conv > INT_MAX || conv < INT_MIN){
//          std::cout << "Invalid input argument " << argv[i+1] << ", default value used for option -n, resp --example" << std::endl; 
//          num = "disc";
//        } else 
//          num = conv; 
      } else if (std::string(argv[i]) == "-r" || std::string(argv[i]) == "--refine_steps"){
        char* p;
        long conv = strtol(argv[i+1], &p, 10);
        if (*p != '\0' || conv > INT_MAX || conv < INT_MIN) {
          std::cout << "Invalid input argument " << argv[i+1] << ", default value used for option -r, resp --refine_steps" << std::endl; 
          ref = 10; 
        } else 
          ref = conv;
      } else if (std::string(argv[i]) == "-o" || std::string(argv[i]) == "--order_ele") {
        char* p; 
        long conv = strtol(argv[i+1], &p, 10);
        if(*p != '\0' || conv > INT_MAX || conv < INT_MIN){
          std::cout << "Invalid input argument " << argv[i+1] << ", default value used for option -o, resp --order_ele" << std::endl; 
          order = 2;
        } else 
          order = conv; 
      }
      else {
        std::cout << "Invalid argument specifier " << argv[i] << std::endl;
        exit(1);
      }
    }
  }


  print_example(num, ref, order);
}







