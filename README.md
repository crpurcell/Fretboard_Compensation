# Fret, Nut and Saddle Compensation for Correct Intonation

The notebook here seeks to calculate measurements for a compensated
nut and saddle for a stringed instrument. A compensated nut is
especially important for a short instrument like a ukulele, which can
suffer from extreme intonation problems on the first few frets.

The two jupiter notebooks contain detailed explanations and
fundamental code exploring the physics of intonation offsets due to 1)
depressing a string to a fret and 2) the effect of string
stiffness. The results of these investigations have been used to
implement a python class that produces a scale diagram of a
compensated fretboard for use by luthiers.

## Using fretCalc.py

The ```fretCalc.py``` code defines classes that can be used to
determine optimal nut and bridge offsets that minimise the intonation
error of each string of an instrument. The code is meant to be used
interactively, e.g., from an iPython terminal or Jupyter notebook, and
is capable of producing a scale fretboard diagram.

### Note on string parameters:

Values of elasticity (Young's modulus)
for nylon strings range from 4.41 GPa - 5.40 GPa, with the former
being the value used in the literature for Aquila branded strings
(Wong 2014).

d_A4 = 0.56
d_E4 = 0.73
d_C4 = 0.89
d_G4 = 0.62

```
# Import the instrument ans string classes
from fretCalc import *

# Define an instance of the instrument class
myUke = instrument(nFrets=18, scaleLength_m=0.367, actionNut_m=0.0015, actionFret12_m=0.0035, fingerWidth_m=0.005, stringSep_m=0.007)

# Add A4 string
myUke.add_string(freqOpen_Hz=440.0, massPerLength_kgm=1e-4, stringDiameter_m=0.00050, elasticity_Pa=4.41e9)

# Add E4 string
myUke.add_string(freqOpen_Hz=329.63, massPerLength_kgm=2e-4, stringDiameter_m=0.00073, elasticity_Pa=4.41e9)

# Add C4 string
myUke.add_string(freqOpen_Hz=261.63, massPerLength_kgm=4e-4, stringDiameter_m=0.00089, elasticity_Pa=4.41e9)

# Add G4 string
myUke.add_string(freqOpen_Hz=392.0, massPerLength_kgm=1e-4, stringDiameter_m=0.00062, elasticity_Pa=4.41e9)

# Fit the fretboard
myUke.fit_fretboard()


#
a  = calc_delta()
b  = calc_delta([-2.5/1000, 2.5/1000])
plt.plot(a[1:])
plt.plot(b[1:])
plt.show()

# Get the chisq function
calc_chisq = myUke._get_chisq_func(myUke.stringLst[0])
print(calc_chisq())
print(calc_chisq([-2.51/1000, 2.5/1000]))






# Fit the fretboard
myUke.fit_fretboard()