# pl_testfunctions.py
#
# Defines several tests functions to
# integrate over with known answers.

import numpy as np
from scipy.special import airy
from scipy.misc import derivative

# CURVE FUNCTIONS:
#
# Define curves to integrate over (currently
# only the real line)

def real_line(x):
    return x

# INTEGRATION FUNCTIONS:
#
# Define functions to be integrated
# for which there are known answers

def AiryWittenexp(x, lamb):
    return lamb*1j*(x**3/3 - x)

def AiryWittenint(x, lamb):
    return np.exp(AiryWittenexp(x, lamb))

def Airyexp(x, lamb):
    return 1j*((1/3)*x**3 + lamb*x)

def Airyint(x, lamb):
    return np.exp(Airyexp(x, lamb))

def Gaussexp(x, lamb):
    return 1j*lamb*x**2

def Gaussint(x, lamb):
    return np.exp(Gaussexp(x, lamb))

# EXACT FUNCTIONS:
#
# Define answer functions for the given integrals.
# Note we still don't have an answer function for
# Witten's ""Airy"" integral.

def Gaussexact(lamb):
    return np.sqrt(np.pi*1j/lamb)

def Airyexact(lamb):
    return airy(lamb)[0]*2*np.pi
