// Source code adapted from “SystemC Training: The Definitive Guide to SystemC” by Doulos Ltd., 2015
#include "mult.h"

void Mult::action(){
  f = a.read() * b.read();
}
