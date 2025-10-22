// Source code adapted from “SystemC Training: The Definitive Guide to SystemC” by Doulos Ltd., 2015
#include "mon.h"
#include <iostream>
#include <cassert>

void Mon::monitor(){
  while (true)
  {
    wait();

    int a_val = A.read();
    int b_val = B.read();
    int f_val = F.read();
    int expected = a_val * b_val;

    std::cout << sc_time_stamp() << " | "
              << "A=" << a_val << " "
              << "B=" << b_val << " "
              << "F=" << f_val
              << " (expected " << expected << ")" << std::endl;

    sc_assert(f_val == expected);
  }
}
