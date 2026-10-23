* 3rd order, 2dB ripple, 5GHz cutoff, Z0 = 50
* consult Ngspice user manual for component formats

Vin 1 0 pulse       $ complete this line

Rin 1 2 50          $ a resistor connected between nodes 1 and 2, with value 50 Ohm
T1                  $ complete the following five lines
L1 
C2 
L3 
Rout 

.TRAN 1p 3n          $ perform transient analysis with step 1p and run to 2n
.probe v(1) v(2)     $ only the probed point will be stored
.option interp       $ interpolate the waveform to uniform time step
.end                 $ every circuit needs a .end statement at the end

