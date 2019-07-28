#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     fretCalc.py                                                       #
#                                                                             #
# PURPOSE:  Code to calculate the fretboard layout of a stringed instrument.  #
#                                                                             #
# MODIFIED: 28-Jul-2010 by C. Purcell                                         #
#                                                                             #
#=============================================================================#
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import scipy.optimize as op


#-----------------------------------------------------------------------------#
class instrument:
    """
    Class that defines the fundamental properties of a stringed instrument.
    """

    def __init__(self, nFrets=18, scaleLength_m=0.367, actionNut_m=0.0015,
                 actionFret12_m=0.0035, fingerWidth_m=0.005,
                 stringSep_m=0.007):
        """
        Constructor for the instrument class.
        """

        # Geometric parameters the instrument
        self.nFrets = nFrets
        self.scaleLength_m = scaleLength_m
        self.actionNut_m  = actionNut_m
        self.actionFret12_m = actionFret12_m
        self.fingerWidth_m = fingerWidth_m
        self.stringSep_m = stringSep_m

        # Calculate the fret positions for an ideal instrument (no offsets)
        self.vibeLenIdealArr_m = self._calc_vibelen(verbose=True)

        # List to store the string objects added via the add_string method
        self.stringLst = []

    def fit_fretboard(self):
        """
        Fit for nut and saddle offset for each string to minimise the
        intonation errror.
        """

        for strObj in self.stringLst:
            self._calc_fretted_length(strObj)

    def print_fretboard():
        """
        Print a scale diagram of the fretboard.
        """

        pass

    def add_string(self, freqOpen_Hz=440.0, massPerLength_kgm=1e-4,
                   stringDiameter_m=0.0005, elasticity_Pa=4e9):
        """
        Add a string to the instrument.
        """

        # Create a new string object
        strObj = string(freqOpen_Hz, massPerLength_kgm, stringDiameter_m,
                        elasticity_Pa)

        # Initialise the string object with ideal fret spacing
        strObj.vibeLenArr_m = self.vibeLenIdealArr_m.copy()

        # Calculate the ideal frequencies for this string
        strObj.freqIdealArr_Hz = self._calc_freqs(strObj)

        # Append to the list of string objects
        self.stringLst.append(strObj)

    def _calc_vibelen(self, strObj=None, verbose=False):
        """
        Return array of vibrating string length for an ideal instrument
        given a scale length. In the arrays returned, the 1st entry [0]
        is the nut and the last entry [nFrets + 1] is the saddle.
        """

        # If called with no stringObj, default to no offsets
        if strObj is None:
            nutOffset_m = 0
            saddleOffset_m = 0
        else:
            nutOffset_m = strObj.nutOffset_m
            saddleOffset_m = strObj.saddleOffset_m

        # Calculate the ideal vibrating length from each fret to saddle
        fretIndxArr = np.arange(self.nFrets + 2)
        vibeLenArr_m = self.scaleLength_m / (2**(fretIndxArr / 12))
        vibeLenArr_m[0] += nutOffset_m
        vibeLenArr_m += saddleOffset_m
        vibeLenArr_m[-1] = 0.0

        # Calculate the distance from the nut
        distNutArr_m = vibeLenArr_m[0] - vibeLenArr_m

        # Print the values if verbose flag set
        if verbose:
            print("Fret  Nut_Distance  Vibrate_Length")
            print("              (mm)            (mm)")
            print("-" * 36)
            for i in range(len(distNutArr_m)):
                print(" {:2d}: {:13.1f} {:15.1f}".format(
                    i, distNutArr_m[i]*1000, vibeLenArr_m[i]*1000))

        return vibeLenArr_m

    def _calc_freqs(self, strObj):
        """
        Return an array of frequencies for an ideal string, i.e., one
        without any tension increase when depressed.
        """

        # Calculate the fundamental frequency at each fret
        with np.errstate(divide='ignore', invalid='ignore'):
            freqArr_Hz = (strObj.freqOpen_Hz * strObj.vibeLenArr_m[0]
                          / strObj.vibeLenArr_m)

        return freqArr_Hz

    def _calc_fretted_length(self, strObj):
        """
        Return an array of open string lengths corresponding to the string
        being depressed behind each fret. This is the fundamental geometric
        model based on the 'clothes-line' effect.
        """

        # Calculate the slope of the open string and action at saddle
        x1 = strObj.vibeLenArr_m[12]
        x2 = strObj.vibeLenArr_m[0]
        y1 = self.actionFret12_m
        y2 = self.actionNut_m
        m = (y2 - y1) / (x2 - x1)
        actionSaddle_m = self.actionFret12_m - m * strObj.vibeLenArr_m[12]

        # Calculate Lprime: string length  projected to the fingerboard
        Lprime = ( np.sqrt(strObj.vibeLenArr_m[0]**2.0
                           - (actionSaddle_m - self.actionNut_m)**2.0) )

        # Calculate the clipped finger width (if overlaps previous fret)
        fretWidthArr_m = np.zeros_like(strObj.vibeLenArr_m)
        fretWidthArr_m[1:] = -1 * np.diff(strObj.vibeLenArr_m)
        fingerWidthArr_m = np.where(fretWidthArr_m <= self.fingerWidth_m,
                                    fretWidthArr_m, self.fingerWidth_m)
        fingerWidthArr_m[-1] = 0

        # Calculate L1, the string length above the finger press
        L1 = np.sqrt(self.actionNut_m**2.0
                     + (Lprime - fingerWidthArr_m
                        - strObj.vibeLenArr_m)**2.0)

        # Calculate L2, the string length below the finger press
        L2 = np.sqrt(actionSaddle_m**2.0 + strObj.vibeLenArr_m**2.0)

        # Deflected length: nut-to-finger + fingerWidth + fret-to-saddle
        # Set the nut to the open length and saddle to np.nan
        newLArr = L1 + fingerWidthArr_m + L2
        newLArr[0] = strObj.vibeLenArr_m[0]
        newLArr[-1] = np.nan

        return newLArr


#-----------------------------------------------------------------------------#
class string:
    """
    Class to store the properties of a string.
    """

    def __init__(self, freqOpen_Hz, massPerLength_kgm, stringDiameter_m,
                 elasticity_Pa):
        """
        Constructor method for the string class.
        """

        # Fundamental physical properties of the string
        self.freqOpen_Hz = freqOpen_Hz
        self.massPerLength_kgm = massPerLength_kgm
        self.stringDiameter_m = stringDiameter_m
        self.elasticity_Pa = elasticity_Pa

        # Initialise nut and saddle offsets
        self.nutOffset_m = 0.0
        self.saddleOffset_m = 0.0

        # Array of ideal frequencies (no offsets)
        self.freqIdealArr_Hz = None

        # Array of vibrating lengths (including offsets)
        self.vibeLenArr_m = None

