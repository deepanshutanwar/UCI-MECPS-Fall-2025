// Source code adapted from “SystemC Training: The Definitive Guide to SystemC” by Doulos Ltd., 2015
#include "stim.h"

void Stim::stimulus(){
  const int a_vals[] = {1, 2, 3, 6, 7, 14, 21, 42};
  const int b_vals[] = {42, 21, 14, 7, 6, 3, 2, 1};

  wait();

  for (int i = 0; i < 8; ++i){
    A.write(a_vals[i]);
    B.write(b_vals[i]);
    wait();
  }

  sc_stop();
}
