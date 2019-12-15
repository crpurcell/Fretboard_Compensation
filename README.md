# Fret, Nut and Saddle Compensation for Correct Intonation

The Jupyter notebooks here describe how to calculate nut and saddle
offsets to optimise the intonation for a stringed instrument. The
error in frequency tends to increase (sharper notes) towards higher
frets because of increasing distance between the string and
fretboard. Almost perfect intonation can be achieved by 1) offsetting
the saddle further from the nut, which tends to flatten the error
curve, 2) offsetting the nut, which tends to shift the error curve down
in frequency and 3) minimising the action at the nut, which prevents
intonation issues at the first few frets.

The two Jupyter notebooks contain detailed explanations and python
code exploring the physics of intonation offsets due to 1) depressing
a string to a fret and 2) the effect of string stiffness. The results
of these investigations have been used to implement a python code
```fretCalc.py``` that produces a scale diagram of a compensated
fretboard for use by luthiers.

## Using fretCalc.py

The ```fretCalc.py``` code defines classes that can be used to
fit for optimal nut and bridge offsets that minimise the intonation
error of each string of an instrument. The code is meant to be used
interactively, e.g., from an iPython terminal or Jupyter notebook, and
is capable of producing a scale fretboard diagram.

Values of elasticity (Young's modulus) for nylon strings range from
4.41 GPa - 5.40 GPa, with the former being the value used in the
literature for Aquila branded strings (Wong 2014). However, string
parameters should ideally be measured directly (notes coming soon!).

The following shows how to run ```fretCalc.py``` from an interactive
python shell or notebook:

```
# Import the instrument and string classes
from fretCalc import *

# Define an instance of the instrument class
myUke = instrument(nFrets=18, scaleLength_m=420.0/1000, actionNut_m=1.5/1000, \
                   actionFret12_m=3.5/1000, fretHeight_m=1.0/1000, \
                   depressHeight_m=0.5/1000, depressMod=0.2, stringSep_m=7/1000)

# Add A4 string
myUke.add_string(freqOpen_Hz=440.00, massPerLength_kgm=1e-4, \
                 stringDiameter_m=0.56/1000, elasticity_Pa=4.41e9)

# Add E4 string
myUke.add_string(freqOpen_Hz=329.63, massPerLength_kgm=2e-4, \
                 stringDiameter_m=0.73/1000, elasticity_Pa=4.41e9)

# Add C4 string
myUke.add_string(freqOpen_Hz=261.63, massPerLength_kgm=4e-4, \
                 stringDiameter_m=0.89/1000, elasticity_Pa=4.41e9)

# Add G4 string
myUke.add_string(freqOpen_Hz=392.00, massPerLength_kgm=1e-4, \
                 stringDiameter_m=0.62/1000, elasticity_Pa=4.41e9)

# Fit the fretboard
myUke.fit_fretboard()
```