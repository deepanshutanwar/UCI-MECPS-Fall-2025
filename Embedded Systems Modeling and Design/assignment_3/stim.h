// Source code adapted from “SystemC Training: The Definitive Guide to SystemC” by Doulos Ltd., 2015
#include "systemc.h"

SC_MODULE(Stim){
  sc_in<bool> Clk;
  sc_out<int> A;
  sc_out<int> B;

  void stimulus();

  SC_CTOR(Stim)
  {
    SC_THREAD(stimulus);
    sensitive << Clk.pos();
  }
};