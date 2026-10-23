* 3rd order, 2dB ripple, 5GHz cutoff, Z0 = 50
* consult Ngspice user manual for component formats

Vin 1 0 pulse(0 2 0 50p 50p 1u 10u)       $ complete this line

Rin 1 2 50          $ a resistor connected between nodes 1 and 2, with value 50 Ohm
T1  2 0 3 0 Z0=50 TD=250ps                $ complete the following five lines
L1 3 4 4.314n
C2 4 0 0.5301p
L3 4 5 4.314n
Rout 5 0 50

.control
tran 1p 3n          $ perform transient analysis with step 1p and run to 2n
*.probe v(5) v(2)     $ only the probed point will be stored
*.print par'v(2)-v(1)/2'
*.op
*.option post = 3
*.option delmax = 1p
*.option interp       $ interpolate the waveform to uniform time step
linearize v(2)
set wr_singlescale
wrdata prbs.txt v(2)
.endc
.end                 $ every circuit needs a .end statement at the end

