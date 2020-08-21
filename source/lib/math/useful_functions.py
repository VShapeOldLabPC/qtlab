#  file with useful functions to send waveform to the AWG Tabor
import numpy as np

# class TaborUsefulFunctions:
def volt2bit(volt):
    """
        Return the bit code corresponding to the entered voltage value in uint16
    """
    full = 4. # in volt
    resolution = 2**14. - 1.

    return  np.array(np.round((volt + full/2.)*resolution/full, 0),
                     dtype ='uint16')

def volt2bit_2(volt):
    """
        Return the bit code corresponding to the entered voltage value in uint16
    """
    full = 2. # in volt
    resolution = 2**14. - 1.

    return  np.array(np.round((volt + full/2.)*resolution/full, 0),
                     dtype ='uint16')

def pulse(p, x, type='DC'):
    start, width, amplitude = p

    # from second to sample
    time_step = x[1] - x[0]
    start = int(round(start/time_step))
    width = int(round(width/time_step))

    after  = len(x) - width - start

    pulse = np.concatenate((np.zeros(start),
                            np.ones(width),
                            np.zeros(after)))

    return amplitude*pulse

def cos(p, x):
    start, duration, amplitude, frequency = p

    # from second to sample
    time_step = x[1] - x[0]
    start = int(round(start/time_step))
    width = int(round(duration/time_step))

    after  = len(x) - width - start

    pulse = np.concatenate((np.zeros(start),
                            amplitude*np.cos(2.*np.pi*frequency*x[start:-after]),
                            np.zeros(after)))

    return pulse
