// Source code adapted from “SystemC Training: The Definitive Guide to SystemC” by Doulos Ltd., 2015
#include "systemc.h"
#include "stim.h"
#include "mult.h"
#include "mon.h"

SC_MODULE(Top){
  sc_signal<int> asig, bsig, fsig;
  sc_clock       testclk;

  Stim stim1;
  Mult mult1;
  Mon  mon1;

  SC_CTOR(Top)
  : testclk("testclk", 10, SC_NS),
    stim1("stim1"),
    mult1("mult1"),
    mon1("mon1")
  {
    stim1.Clk(testclk);
    stim1.A(asig);
    stim1.B(bsig);

    mult1.a(asig);
    mult1.b(bsig);
    mult1.f(fsig);

    mon1.Clk(testclk);
    mon1.A(asig);
    mon1.B(bsig);
    mon1.F(fsig);
  }
};