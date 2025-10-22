Topic: Introduction to SystemC language and simulation

Overview->
In this project, we used three modules: Stimulus, Multiplier and Monitor
> Stimulus: It sends a set of test vectors on each positive clock edge.
we have used these pairs
(1,42), (2,21), (3,14), (6,7), (7,6), (14,3), (21,2), and (42,1).

> Multiplier: It performs multiplication of input A and B

> Monitor: It checks the results on the falling edge of the clock, it prints A, B and F values.

Issue faced->
Faced errors because of the port name (A vs a, Clk vs clk), after fixing this, we got our expected result

Result->
        Copyright (c) 1996-2018 by all Contributors,
        ALL RIGHTS RESERVED
5 ns | A=1 B=42 F=42 (expected 42)
15 ns | A=2 B=21 F=42 (expected 42)
25 ns | A=3 B=14 F=42 (expected 42)
35 ns | A=6 B=7 F=42 (expected 42)
45 ns | A=7 B=6 F=42 (expected 42)
55 ns | A=14 B=3 F=42 (expected 42)
65 ns | A=21 B=2 F=42 (expected 42)
75 ns | A=42 B=1 F=42 (expected 42)

