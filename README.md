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

# Necessary imports
from fretCalc import *

# Define an instance of the instrument class
myUke = instrument(nFrets=18, scaleLength_m=0.367, actionNut_m=0.0015, actionFret12_m=0.0035, fingerWidth_m=0.005, stringSep_m=0.007)

# Add a string
myUke.add_string(freqOpen_Hz=440.0, massPerLength_kgm=1e-4, stringDiameter_m=0.0005, elasticity_Pa=4e9)

# Fit the fretboard
myUke.fit_fretboard()

# Get the model function
calc_delta = myUke._get_model_func(myUke.stringLst[0])

#
a  = calc_delta()
b  = calc_delta([-2.5/1000, 2.5/1000])
plt.plot(a[1:])
plt.plot(b[1:])
plt.show()
