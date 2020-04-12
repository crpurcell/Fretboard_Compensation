#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     fretCalc.py                                                       #
#                                                                             #
# PURPOSE:  Code to calculate the fretboard layout of a stringed instrument.  #
#                                                                             #
# MODIFIED: 12-Apr-2020 by C. Purcell                                         #
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

    def __init__(self, nFrets=18, scaleLength_m=420.0/1000,
                 actionNut_m=1.5/1000, actionFret12_m=3.5/1000,
                 fretHeight_m=1.0/1000, depressHeight_m=0.5/1000,
                 depressMod=0.2, stringSep_m=7/1000):
        """
        Constructor for the instrument class.
        """

        # Fixed geometric parameters of the instrument
        self.nFrets = nFrets
        self.scaleLength_m = scaleLength_m
        self.actionNut_m  = actionNut_m
        self.actionFret12_m = actionFret12_m
        self.fretHeight_m = fretHeight_m
        self.depressHeight_m = depressHeight_m
        self.depressMod = depressMod
        self.stringSep_m = stringSep_m

        # Calculate the fret positions for an ideal instrument (no offsets)
        self.vibeLenIdealArr_m = self._calc_vibelen()

        # List to store the string objects added via the add_string method
        self.stringLst = []

    def add_string(self, freqOpen_Hz=440.0, massPerLength_kgm=1e-4,
                   stringDiameter_m=0.5/1000, elasticity_Pa=4e9):
        """
        Add a string to the instrument and calculate basic parameters.
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

    def fit_fretboard(self):
        """
        Fit for nut and saddle offset for each string to minimise the
        intonation errror.
        """

        # Loop through each of the strings in turn and fit
        for strObj in self.stringLst:

            # Get a function to evaluate chi-squared for current string
            # Free parameters are p = [nutOffset_m, saddleOffset_m]
            calc_chisq = self._get_chisq_func(strObj)

            # Optimise p to minimise chi-squared
            p0 = [0.0, 0.0]
            retMatrix = op.fmin(calc_chisq, p0, full_output=1, disp=False,
                                retall=False)

            # Parse the return values and store in the string object
            p = retMatrix[0]
            chiSq = retMatrix[1]
            chiSqRed = chiSq/(self.nFrets-len(p)-1)
            nIter = retMatrix[2]
            warnFlag = retMatrix[4]
            strObj.nutOffset_m = p[0]
            strObj.saddleOffset_m = p[1]

            # Feedback to user
            print("-" * 80)
            print("Frequency           = {:.2f} Hz".format(strObj.freqOpen_Hz))
            print("Num. Iterations     = {:d}".format(nIter))
            print("Warning Flag        = {:d}".format(warnFlag))
            print("chi-Squared         = {:.1f}".format(chiSq))
            print("chi-Squared Reduced = {:.1f}".format(chiSqRed))
            print("Offsets:            = {:.2f} mm, {:.2f} mm".format(
                p[0]*1000, p[1]*1000))

    def print_fretboard(self):
        """
        Print the fretboard parmeters.
        """

        print("Action Nut (mm)     = {:.1f} Hz".format(self.actionNut_m * 1000))

    def print_fretboard_diagram(self):
        """
        Print a scale diagram of the fretboard.
        """

        pass
        # Print common parameters

    def _calc_vibelen(self, strObj=None, verbose=False):
        """
        Return array of vibrating string length for an ideal instrument
        given a scale length. In the arrays returned, the 1st entry [0]
        is the nut and the last entry [nFrets + 1] is the saddle.
        """

        # If called with no stringObj, default to no offsets
        if strObj is None:
            nutOffset_m = 0.0
            saddleOffset_m = 0.0
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
        strObj.actionSaddle_m = (self.actionFret12_m
                                 - m * strObj.vibeLenArr_m[12])

        # Intermediate variables
        a = self.fretHeight_m
        b = self.actionNut_m - self.fretHeight_m
        c = strObj.actionSaddle_m - self.actionNut_m

        # We assume that vibeLenArr contains only ideal, projectected lengths
        # i.e., already projected onto the fretboard coord system.
        # Need to calculate Lopen as the projection onto the string
        Lproj = strObj.vibeLenArr_m[0]
        Lopen = np.sqrt(strObj.vibeLenArr_m[0]**2.0 + c**2.0)

        # Calculate the fret width array
        fretWidthArr = np.zeros_like(strObj.vibeLenArr_m)
        fretWidthArr[1:] = -1 * np.diff(strObj.vibeLenArr_m)
        fretWidthArr[-1] = np.nan

        # Calculate l1, the string length above the finger press.
        # At 1st fret l1 = 0 because there is no fret closer to nut.
        l1 = np.sqrt(b**2.0 + (Lproj - fretWidthArr - strObj.vibeLenArr_m)**2.0)
        l1[0:2] = 0

        # Calculate l3, the string length below the finger press
        # At nut (index=0) l3 = Lopen
        l3 = np.sqrt((b + c)**2.0 + strObj.vibeLenArr_m**2.0)
        l3[0] = Lopen

        # Calculate the depressHeight array, which becomes less as the
        # frets get closer together. The depressMod term decreases the
        # strength of the correction with fret-number, with 0 = light
        # finger (less depression at higher frets) and 1 = heavy finger
        # (full depression at higher frets).
        depressHeightArr = np.ones_like(fretWidthArr) * self.depressHeight_m
        depressHeightArr[-1] = np.nan
        depressFrac = fretWidthArr/fretWidthArr[1]
        depressFrac += (1-depressFrac) * self.depressMod
        depressHeightArr *= depressFrac

        # Calculate g, the offset from the fret
        gArr = fretWidthArr / 2

        # Calculate l2, the sharp finger deflection between frets
        l2 = (np.sqrt(depressHeightArr**2.0 + (fretWidthArr - gArr)**2.0) +
              np.sqrt(depressHeightArr**2.0 + gArr**2.0))
        l2[0] = 0

        # Calculate l2 at the 1st fret, which runs directly to the nut
        l2[1] = (np.sqrt((depressHeightArr[1] + b)**2.0 +
                         (fretWidthArr[1] - gArr[1])**2.0) +
                 np.sqrt(depressHeightArr[1]**2.0 + gArr[1]**2.0))

        # Deflected length: nut-to-finger + fingerWidth + fret-to-saddle
        # Set the nut to the open length and saddle to np.nan
        frettedLenArr_m = l1 + l2 + l3
        frettedLenArr_m[-1] = np.nan

        return frettedLenArr_m

    def _calc_tension_mult(self, strObj):
        """
        Return an array of (tau + delta_tau) / given the open string
        frequency, scale-length, string properties and the increase in
        string length due to bending.
        """

        # Calculate the tension in the open string
        lenOpen_m = strObj.frettedLenArr_m[0]
        tauOpen_N = (strObj.massPerLength_kgm
                     * np.power(2 * lenOpen_m * strObj.freqOpen_Hz, 2))

        # Calculate the slack length
        slackL_m = lenOpen_m / (tauOpen_N
                                / (strObj.elasticity_Pa * strObj.area_m2) + 1)

        # Calculate the tension in the deflected string
        strObj.frettedTauArr_N = (strObj.elasticity_Pa * strObj.area_m2
                                  * (strObj.frettedLenArr_m - slackL_m)
                                  / slackL_m)

        # The frequency multiplier is proportional to sqrt(tension)
        return np.sqrt(strObj.frettedTauArr_N / tauOpen_N)

    def _calc_stiff_mult(self, strObj, n=1):
        """
        Calculate the stiffness multiplier: the increase (f + delta_f) / f
        of the frequency of the nth harmonic due to the resistance of the
        string to bending motions. Assume n=1 dominates perception of tuning.
        """

        # Calculate the tension in the open string
        lenOpen_m = strObj.vibeLenArr_m[0]
        tau_N = (strObj.massPerLength_kgm
                     * np.power(2 * lenOpen_m * strObj.freqOpen_Hz, 2))

        # Note: More correct to use the deflected tension, but small effect
        #tau_N =  strObj.frettedTauArr_N

        # Calculate the alpha and beta stretch terms
        alpha = 4 + (n**2 * np.pi**2) / 2
        kappa = strObj.stringRadius_m / 2
        beta = np.sqrt(strObj.elasticity_Pa * strObj.area_m2
                       * kappa**2 / tau_N)

        # Calculate stiffness multiplier
        with np.errstate(divide='ignore', invalid='ignore'):
            stiffFreqMultArr = (1 + 2 * beta / strObj.vibeLenArr_m
                                + alpha * beta**2 / strObj.vibeLenArr_m**2)

        # Assuming tuned to open string, so normalise to 1 at nut
        stiffFreqMultArr /= stiffFreqMultArr[0]
        stiffFreqMultArr[stiffFreqMultArr == np.inf] = np.nan

        return stiffFreqMultArr

    def _get_model_func(self, strObj):
        """
        Return a function to calculate the frequecny offset for the
        instrument given a string object. The returned function takes a
        vector of p = [nutOffset_m, saddleOffset_m] as an argument.
        """

        # Return function takes only a vector of free parameters
        def calc_deltafreq_cent(p=[0, 0]):
            """
            Model function to calculate intonation given nut & bridge offset.

            p = [nutOffset_m, saddleOffset_m]
            """

            # Update the offsets in the string object
            strObj.nutOffset_m, strObj.saddleOffset_m = p

            # Calculate the fret positions including offsets
            strObj.vibeLenArr_m = self._calc_vibelen(strObj)

            # Calculate the ideal frequencies including offsets
            strObj.freqArr_Hz = self._calc_freqs(strObj)

            # Calculate the deflected length
            strObj.frettedLenArr_m = self._calc_fretted_length(strObj)

            # Calculate the tension multiplier
            strObj.tauFreqMultArr = self._calc_tension_mult(strObj)

            # Calculate the stiffness multiplier
            strObj.stiffFreqMultArr = self._calc_stiff_mult(strObj)

            # Calculate the deviation from the ideal frequencies in cents
            freqNewArr_Hz = (strObj.freqArr_Hz * strObj.stiffFreqMultArr
                             * strObj.tauFreqMultArr)
            deltaFreq_cent = 1200 * np.log2(freqNewArr_Hz
                                            / strObj.freqIdealArr_Hz)

            return deltaFreq_cent

        return calc_deltafreq_cent

    def _get_chisq_func(self, strObj, weightArr=None):
        """
        Return a function to calculate the chi-squared of the current
        instrument given a string object. The returned function takes a
        vector of p = [nutOffset_m, saddleOffset_m] as an argument.
        """

        # Return function takes only a vector of free parameters
        def calc_chisq(p=[0, 0]):
            """
            Function to caclulate chi-squared given nut & bridge offset.

            p = [nutOffset_m, saddleOffset_m]
            """

            # Get the model function
            model = self._get_model_func(strObj)

            # Calculate the intonation error array
            deltaFreq_cent = model(p)

            # Set the weight to unity by default
            if weightArr is None:
                weights = np.ones_like(deltaFreq_cent)
            else:
                weights = weightArr

            # Calculate chi-squared
            chiSq = np.nansum( (deltaFreq_cent[1:]**2/weights[1:])**2 )

            return chiSq

        return calc_chisq

#-----------------------------------------------------------------------------#
class string:
    """
    Class to store the fixed and dynamic properties of a string.
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
        self.stringRadius_m = stringDiameter_m / 2
        self.area_m2 = np.pi * self.stringRadius_m**2.0

        # Array of ideal frequencies (no offsets)
        self.freqIdealArr_Hz = None

        #------- Variables below here are dynamic; depend on offsets ---------#

        # Initialise nut and saddle offsets
        self.nutOffset_m = 0.0
        self.saddleOffset_m = 0.0
        self.actionSaddle_m = None

        # Array of vibrating lengths that define fretboard
        self.vibeLenArr_m = None

        # Array of ideal frequencies for current fretboard
        self.freqArr_Hz = None

        # Array of open fretted lengths from the clothesline effect
        self.frettedLenArr_m = None

        # Array of fretted tension
        self.frettedTauArr_N = None

        # Tension multiplier array
        self.tauFreqMultArr = None

        # Stiffness multiplier array
        self.stiffFreqMultArr = None
