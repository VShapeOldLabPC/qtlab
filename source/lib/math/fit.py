 # This Python file uses the following encoding: utf-8
# fit.py, functions for fitting including uncertainty estimation
# and option to keep parameters fixed.
# Reinier Heeres <reinier@heeres.eu>, 2011
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

import numpy as np
from scipy.optimize import leastsq
from numpy.random import rand
from scipy.misc import factorial as factoriel
import code
import copy

WEIGHT_EQUAL    = 0
WEIGHT_10PCT    = 1
WEIGHT_20PCT    = 2
WEIGHT_SQRTN    = 3
WEIGHT_SQRTN2   = 4
WEIGHT_N        = 5
WEIGHT_LOGN     = 6
WEIGHT_LEN      = 7

# Symbolic description of functions for the future.
def eval_func(codestr, **kwargs):
    ret = eval(codestr, kwargs)
    return ret

class Function:
    def __init__(self, xdata=None, ydata=None, xerr=None, yerr=None,
                    weight=WEIGHT_EQUAL, minerr=None, nparams=None):
        '''
        Fitting function class.

        Parameters:
        - xdata, ydata: data vectors
        - xerr, yerr: error vectors; xerr is not used
        - weight: weighting mechanism
        - minerr: minimum absolute error value in case of automatic weight
        generation using e.g. WEIGHT_SQRTN
        - nparams: number of parameters, not checked if not specified
        '''

        self._fixed = {}
        self._nparams = nparams
        self._weight = weight
        self._minerr = minerr
        self.set_data(xdata, ydata, xerr, yerr)

    def set_weight(self, w):
        self._weight = w

    def set_data(self, x, y, xerr=None, yerr=None):
        self._xdata = x
        self._xerr = xerr

        self._ydata = y
        if yerr is not None:
            self._yerr = yerr
        elif self._xdata is None or self._ydata is None:
            self._yerr = None
        elif self._weight == WEIGHT_EQUAL:
            self._yerr = np.ones_like(self._xdata)
        elif self._weight == WEIGHT_10PCT:
            self._yerr = 0.10 * np.abs(self._ydata)
        elif self._weight == WEIGHT_20PCT:
            self._yerr = 0.20 * np.abs(self._ydata)
        elif self._weight == WEIGHT_SQRTN:
            self._yerr = np.sqrt(np.abs(self._ydata))
        elif self._weight == WEIGHT_SQRTN2:
            self._yerr = 0.5 * np.sqrt(np.abs(self._ydata))
        elif self._weight == WEIGHT_N:
            self._yerr = np.abs(self._ydata)
        elif self._weight == WEIGHT_LOGN:
            self._yerr = np.log(np.abs(self._ydata))
        elif self._weight == WEIGHT_LEN:
            self._yerr = (1.+np.arange(len(self._xdata)))

        # Set minimum errors
        if self._yerr is not None and self._minerr is not None:
            self._yerr[self._yerr < self._minerr] = self._minerr

    def set_nparams(self, n):
        if self._nparams not in (n, None):
            raise ValueError('Different number of parameters expected')
        self._nparams = n

    def get_parameters(self, p):
        '''
        Return set of parameters including fixed parameters when given
        either a complete set or a reduced set of only free parameters.
        '''
        if len(p) == self._nparams:
            return p

        p = copy.copy(p)
        for i, v in self._fixed.iteritems():
            p = np.insert(p, i, v)
        return p

    def get_px(self, p, x=None):
        '''
        Return tuple of parameter and x value vector
        '''
        p = self.get_parameters(p)
        if x is None:
            x = self._xdata
        return p, x

    def func(self, p, x=None):
        '''
        Should be implemented in derived classes.
        '''
        pass

    def err_func(self, p):
        residuals = np.abs(self._ydata - self.func(p)) / self._yerr
        return residuals

    def fit(self, p0, fixed=[]):
        '''
        Fit the function using p0 as starting parameters.

        Fixed is a list of numbers specifying which parameter to keep fixed.
        '''

        self.set_nparams(len(p0))

        # Get free parameters
        p1 = []
        for i in range(len(p0)):
            if i not in fixed:
                p1.append(p0[i])

        # Store fixed parameters
        for i in fixed:
            self._fixed[i] = p0[i]

        out = leastsq(self.err_func, p1, full_output=1)
        params = out[0]
        covar = out[1]
        self._fit_params = self.get_parameters(params)
        self._fit_err = np.zeros_like(params)

        if covar is not None:
            dof = len(self._xdata) - len(p1)
            chisq = np.sum(self.err_func(params)**2)
            for i in range(len(params)):
                self._fit_err[i] = np.sqrt(covar[i][i]) * np.sqrt(chisq / dof)

        # Set error of fixed parameters to 0
        for i in fixed:
            self._fit_err = np.insert(self._fit_err, i, 0)

        return self._fit_params

    def fit_odr(self, p0):
        from scipy import odr
        model = odr.Model(self.func)
        data = odr.Data(self._xdata, self._ydata)
        fit = odr.ODR(data, model, p0, maxit=100)
        fit.set_job(fit_type=0) #0 = LSQ, 1 = EXPlicit
        out = fit.run()

        self._fit_params = out.beta
        self._fit_err = out.sd_beta
        return self._fit_params

    def get_fit_params(self):
        return self._fit_params

    def get_fit_errors(self):
        return self._fit_err

    def test_random(self, x0, x1, nx, params, noise, p0=None, logy=False, yerr=None, weight=None):
        '''
        Perform a test with a random data set.
        x0, x1, nx: range x0 to x1 in n steps
        params: parameters to generate y data
        noise: noise level to add
        p0: starting fit parameters, if None use params * 0.8 or params * 1.2
        '''

        x = np.linspace(x0, x1, nx)
        y0 = self.func(params, x)
        y = y0 + (np.random.random([nx]) - 0.5) * noise
        if p0 is None:
            if rand() > 0.5:
                p0 = 0.8 * np.array(params)
            else:
                p0 = 1.2 * np.array(params)

        if weight is not None:
            self.set_weight(weight)
        self.set_data(x, y, yerr=yerr)
        p = self.fit(p0)

        print '\tRandom par: %s' % (pr, )
        s = ''
        for val, err in zip(p, self.get_fit_errors()):
            s += ' %f (+-%f)' % (val, err)
        print '\tResult:%s' % (s, )

        yfit = self.func(p, x)
        plt.errorbar(x, y, yerr=self._yerr, fmt='ks')
        plt.plot(x, yfit, 'r')
        if logy:
            plt.yscale('log')

class Fit3D(Function):

    def __init__(self, x, y, z, *args, **kwargs):
        self._weight = weight
        self._minerr = minerr
        self.set_data(xdata, ydata, xerr, yerr)

    def set_data(self, x, y, z, xerr=None, yerr=None, zerr=None):
        # Pretend z = y for 2D case
        Function.set_data(x, z, yerr=zerr)
        self._zdata = self._ydata
        self._zerr = self._yerr
        self._ydata = y
        self._yerr = yerr

    def func(self, p, x=None, y=None):
        pass

    def err_func(self, p):
        residuals = np.abs(self._zdata - self.func(p)) / self._zerr
        return residuals

class Polynomial(Function):
    '''
    Polynomial fit function of order n

    Polynomial(order=2) creates fit:
                     2
        a + b⋅x + c⋅x
    '''

    def __init__(self, *args, **kwargs):
        self._order = kwargs.pop('order', 2)
        kwargs.setdefault('nparams', self._order + 1)
        Function.__init__(self, *args, **kwargs)

    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret = np.ones_like(x) * p[0]
        for n in range(1, self._order + 1):
            ret += p[n] * x**n

        return ret

class Linear(Polynomial):
    '''
    Linear fit function:
        a + b⋅x
    '''

    def __init__(self, *args, **kwargs):
        kwargs['order'] = 1
        Polynomial.__init__(self, *args, **kwargs)

class Gaussian(Function):
    '''
    Gaussian fit function:
                       ⎛           2   ⎞
              ___      ⎜ -2⋅(-c + x)   ⎜
            ╲╱ 2 ⋅b⋅exp⎜ ───────────── ⎜
                       ⎜        2      ⎜
                       ⎝      d        ⎠
        a + ─────────────────────────────
                            ___
                          ╲╱ π ⋅d

     parameters:
        background
        area
        position
        full width at (exp^(-0.5) = 0.607)
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('nparams', 4)
        Function.__init__(self, *args, **kwargs)

    def get_fwhm(self, p=None):
        if p is None:
            p = self._fit_params
        return np.sqrt(2 * np.log(2)) * p[3]

    def get_height(self, p=None):
        '''Return height (from background).'''
        if p is None:
            p = self._fit_params
        return p[1] / p[3] / np.sqrt(np.pi/2)

    def get_area(self, p=None):
        if p is None:
            p = self._fit_params
        return p[1]

    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret = p[0] + p[1] / p[3] / np.sqrt(np.pi / 2) * np.exp(-2*(x - p[2])**2 / p[3]**2)
        return ret

class DoubleGaussian(Function):
    '''
    Gaussian fit function:
                       ⎛           2   ⎞
              ___      ⎜ -2⋅(-c + x)   ⎜
            ╲╱ 2 ⋅b⋅exp⎜ ───────────── ⎜
                       ⎜        2      ⎜
                       ⎝      d        ⎠
        a + ───────────────────────────── +
                            ___
                          ╲╱ π ⋅d

     parameters:
        background
        area
        position
        full width at (exp^(-0.5) = 0.607)
        area2
        position2
        full width at (exp^(-0.5) = 0.607)2
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('nparams', 7)
        Function.__init__(self, *args, **kwargs)

    def get_fwhm(self, p=None):
        if p is None:
            p = self._fit_params
        return np.sqrt(2 * np.log(2)) * p[3]

    def get_height(self, p=None):
        '''Return height (from background).'''
        if p is None:
            p = self._fit_params
        return p[1] / p[3] / np.sqrt(np.pi/2)

    def get_area(self, p=None):
        if p is None:
            p = self._fit_params
        return p[1]

    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret = p[0] + p[1] / p[3] /\
            np.sqrt(np.pi / 2) * np.exp(-2*(x - p[2])**2 / p[3]**2) + \
             p[4] / p[6] / np.sqrt(np.pi / 2) * np.exp(-2*(x - p[5])**2 / p[6]**2)

        return ret

class GaussianPlain(Function):
    '''
    Gaussian fit function:
                 ⎛           2          ⎞
                 ⎜ -4⋅(-c + x) ⋅ln(2)   ⎜
        a + b⋅exp⎜ ──────────────────── ⎜
                 ⎜           2          ⎜
                 ⎝         d            ⎠

     parameters:
        background
        height
        position
        full width at half max
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('nparams', 4)
        Function.__init__(self, *args, **kwargs)

    def get_fwhm(self, p=None):
        if p is None:
            p = self._fit_params
        return p[3]

    def get_height(self, p=None):
        '''Return height (from background).'''
        if p is None:
            p = self._fit_params
        return p[1]

    def get_area(self, p=None):
        if p is None:
            p = self._fit_params
        return p[1] * np.sqrt(np.pi) * p[3] / np.sqrt(4*np.log(2))

    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret = p[0] + p[1] * np.exp(-4 * np.log(2) * (x - p[2])**2 / p[3]**2)
        return ret

class Lorentzian(Function):
    '''
    Lorentzian fit function:
                   2⋅b⋅d
        a + ───────────────────
              ⎛ 2            2⎞
            π⋅⎝d  + 4⋅(x - c) ⎠

     parameters:
        background
        area
        position
        width
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('nparams', 4)
        Function.__init__(self, *args, **kwargs)

    def get_fwhm(self, p=None):
        if p is None:
            p = self._fit_params
        return p[3]

    def get_height(self, p=None):
        '''Return height (from background).'''
        if p is None:
            p = self._fit_params
        return 2 / np.pi / p[3] * p[1]

    def get_area(self, p=None):
        if p is None:
            p = self._fit_params
        return p[1]

    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret = np.ones_like(x) * p[0] + 2 * p[1] / np.pi * p[3] / (4*(x - p[2])**2 + p[3]**2)
        return ret

class FanoLorentzian(Function):
    '''
    FanoLorentzian fit function:
                               2
                  ⎛c⋅q        ⎞
                b⋅⎜─── + 2⋅π⋅x⎟
                  ⎝ 2         ⎠
        a + ───────────────────────
            ⎛ 2          ⎞
            ⎜c       2  2⎟ ⎛ 2    ⎞
            ⎜── + 4⋅π ⋅x ⎟⋅⎝q  + 1⎠
            ⎝4           ⎠

     parameters:
        background
        amplitude
        width
        Fano factor
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_EQUAL)
        kwargs.setdefault('nparams', 4)
        Function.__init__(self, *args, **kwargs)

    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret =  p[0] + (p[1] * ( p[3]*np.pi*p[2] + 2.*np.pi*(x - p[4]) )**2)/((1 + p[3]**2)*( (2 *np.pi *(x - p[4]))**2 + (np.pi*p[2])**2))
        return ret

class LinearSin(Function):
    '''
    Sine fit function with a linear slop:
        a + b⋅x + c⋅sin(d⋅x + e)

     parameters:
        background
        slop
        amplitude
        frequency
        phi0
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_EQUAL)
        kwargs.setdefault('nparams', 5)
        Function.__init__(self, *args, **kwargs)

    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret =  p[0] + p[1] * x + p[2] * np.sin(x * p[3] + p[4])
        return ret

class Exponential(Function):
    '''
    Exponential fit function:
        a + b⋅exp(-d⋅(x -c))

     parameters:
        background
        amplitude
        displacement
        exponent
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_SQRTN)
        kwargs.setdefault('nparams', 4)
        Function.__init__(self, *args, **kwargs)

    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret = np.ones_like(x) * p[0] + p[1] * np.exp(-(x - p[2]) * p[3])
        return ret

class Sine(Function):
    '''
    Sine fit function:
        a + b⋅sin(x⋅c + d)

     parameters:
        background
        amplitude
        pulsation
        phi0
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_EQUAL)
        kwargs.setdefault('nparams', 4)
        Function.__init__(self, *args, **kwargs)

    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret = np.ones_like(x) * p[0] + p[1] * np.sin(x * p[2] + p[3])
        return ret

class S21dB_dip_asymetric(Function):
    '''
    Asymetrical S21 function.
    See Étienne Dumur thesis manuscript appendix B.
                        ⎛2⋅ⅈ⋅Qi⋅(-f₀ + x)    ⎞
                     Z₀⋅⎜       ───────── + 1⎟
                        ⎝            f₀      ⎠
        ───────────────────────────────────────────────────
                    ⎛2⋅ⅈ⋅Qi⋅(-f₀ + x)       Qi⋅(ⅈ⋅Xe + Z₀)⎞
        (ⅈ⋅Xe + Z₀)⋅⎜       ───────── + 1 + ──────────────⎟
                    ⎝            f₀             Qc⋅Z₀     ⎠

     parameters:
        Internal quality factor
        External quality factor
        Resonance frequency [GHz]
        Asymetry (could be positive or negative) [ohm]
        background [dB]
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_EQUAL)
        kwargs.setdefault('nparams', 6)
        Function.__init__(self, *args, **kwargs)

    def func(self, p, x=None, Z0=50., output='dB'):
        p, x = self.get_px(p, x)

        Qi = p[0]
        Qc = p[1]
        f0 = p[2]
        Xe = p[3]
        backgrounddB = p[4]
        backgroundPhase = p[5]

        a = Z0/(Z0 + 1j*Xe)
        b = (1. + 2.*1j*Qi*(x-f0)/f0)/\
            (1. + 2.*1j*Qi*(x-f0)/f0 + (Z0 + 1j*Xe)*Qi/Qc/Z0)
        y = a*b

        if output == 'dB':
            return 20.*np.log10(abs(y)) + backgrounddB
        elif output == 'phase':

            return np.angle(y.real + 1j*y.imag) + backgroundPhase
        elif output == 'amplitude':

            return abs(y)
        elif output == 'real':

            return y.real
        elif output == 'imag':

            return y.imag

class S21dB_pic_amplitude(Function):
    '''
    resonance S21 pic function. Only in amplitude.
    parameters:
        Internal quality factor
        External quality factor
        Resonance frequency [GHz]
        background
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_EQUAL)
        kwargs.setdefault('nparams', 4)
        Function.__init__(self, *args, **kwargs)

    def func(self, p, x=None):
        p, x = self.get_px(p, x)

        Qi = p[0]
        Qc = p[1]
        f0 = p[2]
        background = p[3]

        y = 1. /(1. + 2.*1j*Qc*(x-f0)/f0 + Qc/Qi)
        # we have to check about signs
        return abs(y) +background

class Transmission_PL(Function):
    '''
    resonance S21 pic function. Only in amplitude.
    parameters:
        kappa_ext [GHz]
        kappa_loss [GHz]
        Resonance frequency [GHz]
        attenuation
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_EQUAL)
        kwargs.setdefault('nparams', 5)
        Function.__init__(self, *args, **kwargs)

    def func(self, p, x=None):
        p, x = self.get_px(p, x)

        kappa_1 = p[0]
        kappa_2 = p[0]
        kappa_loss = p[2]
        f0 = p[3]
        bg = p[4]

        y = 2.*np.sqrt(kappa_1*kappa_2) /(kappa_1+kappa_2+kappa_loss -2j*(x-f0))
        return abs(y)+bg

class dip_lorentzian_amplitude(Function):
    '''
    resonance S21 pic function. Only in amplitude.
    parameters:
        Internal quality factor
        External quality factor
        Resonance frequency [GHz]
        background
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_EQUAL)
        kwargs.setdefault('nparams', 4)
        Function.__init__(self, *args, **kwargs)

    def func(self, p, x=None):
        p, x = self.get_px(p, x)

        Qi = p[0]
        Qc = p[1]
        f0 = p[2]
        background = p[3]

        y = 1. /(1. + 2.*1j*Qc*(x-f0)/f0 + Qc/Qi)
        # we have to check about signs
        return -abs(y) +background

class dip_lorentzian_cauchy(Function):
    '''
    resonance dip function. Only in amplitude.
    parameters:
        FWHM
        Resonance frequency [GHz]
        background
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_EQUAL)
        kwargs.setdefault('nparams', 4)
        Function.__init__(self, *args, **kwargs)

    def func(self, p, x=None):
        p, x = self.get_px(p, x)

        fwhm = p[0]
        f0 = p[1]
        background = p[2]
        amplitude = p[3]

        y = amplitude*2. /np.pi/fwhm/(1. + 4.*(x-f0)**2/fwhm**2 )
        # we have to check about signs
        return -y + background

class ExponentialDecaySine(Function):
    '''
    Sine fit function:
        a + b⋅sin(x⋅c + d)⋅exp(-x/e)

     parameters:
        background
        amplitude
        pulsation
        phi0
        decay time
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_EQUAL)
        kwargs.setdefault('nparams', 5)
        Function.__init__(self, *args, **kwargs)

    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret = np.ones_like(x) * p[0] + p[1] * np.sin(x * p[2] + p[3])*np.exp(-x/p[4])
        return ret

class ExponentialDecayDoubleSine(Function):
    '''
    Sine fit function:
        a + (b⋅sin(x⋅c + d)+b1⋅sin(x⋅c1 + d1))⋅exp(-x/e)

     parameters:
        background
        amplitude0
        pulsation0
        phi0
        amplitude1
        pulsation1
        phi1
        decay time
    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_EQUAL)
        kwargs.setdefault('nparams', 8)
        Function.__init__(self, *args, **kwargs)

    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret = np.ones_like(x) * p[0] + (p[1] * np.sin(x * p[2] + p[3])+ p[4] * np.sin(x * p[5] + p[6]))*np.exp(-x/p[7])
        return ret

class Fit3D_lineartest(Function):

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_EQUAL)
        kwargs.setdefault('nparams', 3)
        Function.__init__(self, *args, **kwargs)
        # self._weight = weight
        # self._minerr = minerr
        # self.set_data(x, y, z)

    def set_data3D(self, x, y, z, xerr=None, yerr=None, zerr=None):
        # Pretend z = y for 2D case
        self.set_data(x, z, yerr=zerr)
        self._zdata = self._ydata
        self._zerr = self._yerr
        self._ydata = y
        self._yerr = yerr

    def func(self, p, x=None, y=None):
        p, x = self.get_px(p, x)
        y = self._ydata

        ret = p[0] + p[1]*x + p[2]*y
        return ret

    def err_func(self, p):
        residuals = np.abs(self._zdata - self.func(p)) / self._zerr
        return residuals

class Fit3D_lineshape(Function):
    '''
    parameters:
        background
        qubit frequency [Hz]
        T_phi [s]
        kappa [Hz]
        Xi [Hz]
        alpha (conversion between power and photon number)
        amp
    '''
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_EQUAL)
        kwargs.setdefault('nparams', 5)
        Function.__init__(self, *args, **kwargs)
        self.N_sum = 1

    def set_data3D(self, x, y, z, xerr=None, yerr=None, zerr=None):
        # Pretend z = y for 2D case
        self.set_data(x, z, yerr=zerr)
        self._zdata = self._ydata
        self._zerr = self._yerr
        self._ydata = y
        self._yerr = yerr

    def set_N_sum(self,n):
        self.N_sum = n

    def func(self, p, x=None, y=None):
        p, x = self.get_px(p, x)
        y = self._ydata
        ret = 0.

        f0 = p[0]
        T2 = p[1]
        kappa = p[2]
        Xi = p[3]
        Lambda = p[4]

        conv = kappa/2 /((kappa/2)**2 + Xi**2) *Lambda
        n = conv*y
        theta2 = (2*Xi/kappa)**2
        for j in np.arange(self.N_sum):
            # ret += (-1)**j /np.math.factorial(j) \
            #         *(4*p[3]/p[2]*np.sqrt(conv*y + p[5]))**(2*j) \
            #         *(1/p[1] + np.pi*p[2]*( 4*p[3]/p[2]*np.sqrt(conv*y +p[5]) )**2 + j )\
            #         /( 4*np.pi**2*(x - p[0] + (2*conv*y +2*p[5]+ 1)*p[3] )**2  \
            #         + (1/p[1] + np.pi*p[2]*(4*p[3]/p[2]*np.sqrt(conv*y+p[5]))**2 + j)**2  )

            ret = ret + (-4.*theta2*n)**j/np.math.factorial(j)*(1./T2 + kappa*(2*theta2*n + j/2.) )/( (2*np.pi*(x - f0) -(2*n+1)*Xi )**2 + (1./T2 + kappa*(2*theta2*n + j/2.) )**2  )

        return - ret/np.pi

    def err_func(self, p):
        residuals = np.abs(self._zdata - self.func(p)) / self._zerr
        return residuals

# class acstark_lineshape(Function):
#     def __init__(self, *args, **kwargs):
#         kwargs.setdefault('weight', WEIGHT_EQUAL)
#         kwargs.setdefault('nparams', 6)
#         Function.__init__(self, *args, **kwargs)
#         self.N_sum = 1
#
#     def set_N_sum(self,n):
#         self.N_sum = n
#
#     def func(self, p, x=None):
#         p, x = self.get_px(p, x)
#         ret = 0.
#         f0 = p[0]
#         T2 = 8.5e-6
#         # p[1]
#         kappa = 2*np.pi*5e6
#         # p[2]
#         Xi = p[3]
#         n = p[4]
#         theta2 = (2*Xi/kappa)**2
#
#         # version 1: computing factorial can be pretty slow
#         for j in np.arange(self.N_sum):
#
#             ret = ret + (-4.*theta2*n)**j/factoriel(j)*(1./T2 + kappa*(2*theta2*n + j/2.) )/( (2*np.pi*(x - f0) +(2*n+1)*Xi )**2 + (1./T2 + kappa*(2*theta2*n + j/2.) )**2  )
#
#         # version 2: testing a mathematical trick
#         # for j in np.arange(self.N_sum):
#         #     if j > 1:
#         #         Log_fact = 0.
#         #         for jj in np.arange(j):
#         #             Log_fact += np.log(jj+1)
#         #         facto = np.exp(Log_fact)
#         #         ret = ret + (-4.*theta2*n)**j/facto*(1./T2 + kappa*(2*theta2*n + j/2.) )/( (2*np.pi*(x - f0) +(2*n+1)*Xi )**2 + (1./T2 + kappa*(2*theta2*n + j/2.) )**2  )
#         #
#         #     else:
#         #         ret = ret + (-4.*theta2*n)**j*(1./T2 + kappa*(2*theta2*n + j/2.) )/( (2*np.pi*(x - f0) +(2*n+1)*Xi )**2 + (1./T2 + kappa*(2*theta2*n + j/2.) )**2  )
#
#         return - ret/np.pi/p[5]/p[1]

class acstark_lineshape(Function):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_EQUAL)
        kwargs.setdefault('nparams', 3)
        Function.__init__(self, *args, **kwargs)
        self.N_sum = 1
        self.f0 = 2.05015e9
        self.T2 = 8.5e-6
        self.kappa = 2*np.pi*5e6
    def set_f0(self,f0):
        self.f0 = f0

    def set_T2(self,T2):
        self.T2 = T2
    def set_kappa(self,kappa):
        self.kappa = kappa
    def set_N_sum(self,n):
        self.N_sum = n

    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret = 0.
        f0 = self.f0
        T2 = self.T2
        # p[1]
        kappa = self.kappa
        # p[2]
        Xi = p[0]
        n = p[1]
        theta2 = (2*Xi/kappa)**2

        # version 1: computing factorial can be pretty slow
        for j in np.arange(self.N_sum):

            ret = ret + (-4.*theta2*n)**j/factoriel(j)*(1./T2 + kappa*(2*theta2*n + j/2.) )/( (2*np.pi*(x - f0) +(2*n+1)*Xi )**2 + (1./T2 + kappa*(2*theta2*n + j/2.) )**2  )

        # version 2: testing a mathematical trick
        # for j in np.arange(self.N_sum):
        #     if j > 1:
        #         Log_fact = 0.
        #         for jj in np.arange(j):
        #             Log_fact += np.log(jj+1)
        #         facto = np.exp(Log_fact)
        #         ret = ret + (-4.*theta2*n)**j/facto*(1./T2 + kappa*(2*theta2*n + j/2.) )/( (2*np.pi*(x - f0) +(2*n+1)*Xi )**2 + (1./T2 + kappa*(2*theta2*n + j/2.) )**2  )
        #
        #     else:
        #         ret = ret + (-4.*theta2*n)**j*(1./T2 + kappa*(2*theta2*n + j/2.) )/( (2*np.pi*(x - f0) +(2*n+1)*Xi )**2 + (1./T2 + kappa*(2*theta2*n + j/2.) )**2  )

        return - ret/np.pi/p[2]/T2
        # return - ret/np.pi/0.1/T2

class acstark_lineshape_v2(Function):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('weight', WEIGHT_EQUAL)
        kwargs.setdefault('nparams', 2)
        Function.__init__(self, *args, **kwargs)
        self.N_sum = 1
        self.f0 = 2.05015e9
        self.T2 = 8.5e-6
        self.kappa = 2*np.pi*5e6
        self.n = 0
    def set_f0(self,f0):
        self.f0 = f0

    def set_T2(self,T2):
        self.T2 = T2
    def set_kappa(self,kappa):
        self.kappa = kappa
    def set_N_sum(self,n):
        self.N_sum = n

    def set_n(self, n):
        self.n = n
    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret = 0.
        f0 = self.f0
        T2 = self.T2

        # p[1]
        kappa = self.kappa
        # p[2]
        Xi = p[0]
        n = self.n
        theta2 = (2*Xi/kappa)**2

        # version 1: computing factorial can be pretty slow
        for j in np.arange(self.N_sum):

            ret = ret + (-4.*theta2*n)**j/factoriel(j)*(1./T2 + kappa*(2*theta2*n + j/2.) )/( (2*np.pi*(x - f0) +(2*n+1)*Xi )**2 + (1./T2 + kappa*(2*theta2*n + j/2.) )**2  )

        # version 2: testing a mathematical trick
        # for j in np.arange(self.N_sum):
        #     if j > 1:
        #         Log_fact = 0.
        #         for jj in np.arange(j):
        #             Log_fact += np.log(jj+1)
        #         facto = np.exp(Log_fact)
        #         ret = ret + (-4.*theta2*n)**j/facto*(1./T2 + kappa*(2*theta2*n + j/2.) )/( (2*np.pi*(x - f0) +(2*n+1)*Xi )**2 + (1./T2 + kappa*(2*theta2*n + j/2.) )**2  )
        #
        #     else:
        #         ret = ret + (-4.*theta2*n)**j*(1./T2 + kappa*(2*theta2*n + j/2.) )/( (2*np.pi*(x - f0) +(2*n+1)*Xi )**2 + (1./T2 + kappa*(2*theta2*n + j/2.) )**2  )

        return - ret/np.pi/p[1]/T2
        # return - ret/np.pi/0.1/T2

class Lorentzian_sum(Function):
    '''
    Remy
    Lorentzian fit function:
                   2⋅b⋅d
        a + ───────────────────
              ⎛ 2            2⎞
            π⋅⎝d  + 4⋅(x - c) ⎠

     parameters:
        background
        area
        position
        width

    '''

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('nparams', 4)
        Function.__init__(self, *args, **kwargs)

    def get_fwhm(self, p=None):
        if p is None:
            p = self._fit_params
        return p[3]

    def get_height(self, p=None):
        '''Return height (from background).'''
        if p is None:
            p = self._fit_params
        return 2 / np.pi / p[3] * p[1]

    def get_area(self, p=None):
        if p is None:
            p = self._fit_params
        return p[1]

    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret = np.ones_like(x) * p[0] + 2 * p[1] / np.pi * p[3] / (4*(x - p[2])**2 + p[3]**2)
        return ret

class NISTRationalHahn(Function):
    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        ret = (p[0] + p[1] * x + p[2] * x**2 + p[3] * x**3) / (1 + p[4] * x + p[5] * x**2 + p[6] * x**3)
        return ret

class NISTGauss(Function):
    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        p = self.get_parameters(p)
        ret = p[0] * np.exp(-p[1] * x) + p[2] * np.exp(-(x-p[3])**2/p[4]**2) +p[5] * np.exp(-(x - p[6])**2/p[7]**2)
        return ret

class FunctionFit(Function):

    def __init__(self, func, *args, **kwargs):
        self._func = func
        Function.__init__(self, *args, **kwargs)

    def func(self, p, x=None):
        p, x = self.get_px(p, x)
        return self._func(p, x)

def fit(f, xdata, ydata, p0, fixed=[], yerr=None, weight=WEIGHT_EQUAL):
    '''
    Fit function 'f' using p0 as starting parameters. The function should
    take a parameter vector and an x data vector as input, e.g.:

        lambda p, x: p[0] + p[1] * x

    Fixed is a list of numbers specifying which parameter to keep fixed.
    weight specifies the weithing method if no y error vector is specified.

    Returns the fitting class.
    '''

    ff = FunctionFit(f, xdata=xdata, ydata=ydata, yerr=yerr, weight=weight)
    result = ff.fit(p0, fixed)
    return ff

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    plt.figure()
    lin = Linear()
    pr = [rand(), (rand() - 0.5) * 5]
    print 'Linear fit:'
    lin.test_random(-10, 10, 20, pr, 2, [0, 1])

    plt.figure()
    gauss = Gaussian()
    # BG, height, pos, width
    pr = [rand(), 5 * (rand() + 0.1), 10 * (rand() - 0.5), 3 * (rand() + 0.1)]
    print 'Gaussian fit:'
    gauss.test_random(-10, 10, 50, pr, 1)

    plt.figure()
    exp = Exponential()
    # BG, height, pos, exponent
    pr = [rand(), 5 * (rand() + 0.1), 10 * (rand() - 0.5), 2 * (rand() + 0.1)]
    print 'Exponential fit:'
    exp.test_random(-10, 10, 50, pr, 1, logy=True)

    plt.figure()
    sine = Sine()
    # BG, amplitude, frequency, phi0
    pr = [rand(), 3 * (rand() + 0.5), 0.5 * np.pi * (rand() + 0.1), 2 * np.pi * rand()]
    print 'Sine fit:'
    sine.test_random(-10, 10, 50, pr, 2)

    plt.figure()
    data = np.loadtxt('/qtlab/source/lib/math/data/gauss_ref.dat')
    print 'Gauss ref:'
    gauss = Gaussian(data[:,0], data[:,1], weight=WEIGHT_EQUAL)
    p0 = [-1, 10, 2, 0.7]
    p = gauss.fit(p0, fixed=(0,))
    print '\tStart par: %s' % (p0, )
    s = ''
    for val, err in zip(p, gauss.get_fit_errors()):
        print '\t\t%e (+-%e)' % (val, err)

    f = lambda p, x: p[0] + p[1] / p[3] / np.sqrt(np.pi / 2) * np.exp(-2*(x - p[2])**2 / p[3]**2)
    fc = fit(f, data[:,0], data[:,1], p0)
    p = fc.get_fit_params()
    print '\tStart par: %s' % (p0, )
    s = ''
    for val, err in zip(p, gauss.get_fit_errors()):
        print '\t\t%e (+-%e)' % (val, err)

    plt.errorbar(data[:,0], data[:,1], yerr=gauss._yerr, fmt='ks')
    plt.plot(data[:,0], gauss.func(p))

    # Reference data sets from NIST:
    # http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
    plt.figure()
    data = np.loadtxt('/qtlab/source/lib/math/data/Hahn1.dat')
    print 'NIST Hahn:'
    hahn = NISTRationalHahn(data[:,1], data[:,0])
    p0 = [1e1, -1e-1, 5e-3, -1e-6, -5e-3, 1e-4, -1e-7]
    pa = [1.08e0,-1.23e-1,4.09e-3,-1.43e-6,-5.76e-3,2.41e-4,-1.23e-7]
    p = hahn.fit(p0)
    print '\tStart par: %s' % (p0, )
    s = ''
    for val, err in zip(p, hahn.get_fit_errors()):
        print '\t\t%e (+-%e)' % (val, err)

    plt.plot(data[:,1], data[:,0], 'ks')
    plt.plot(data[:,1], hahn.func(p), 'r+')

    plt.figure()
    data = np.loadtxt('/qtlab/source/lib/math/data/Gauss1.dat')
    print 'NIST Gauss:'
    gauss = NISTGauss(data[:,1], data[:,0])
    p0 = [9.7e1,9e-3,1e2,6.5e1,2e1,7e1,1.78e2,1.65e1]
    p = gauss.fit(p0)
    print '\tStart par: %s' % (p0, )
    s = ''
    for val, err in zip(p, gauss.get_fit_errors()):
        print '\t\t%e (+-%e)' % (val, err)

    plt.plot(data[:,1], data[:,0], 'ks')
    plt.plot(data[:,1], gauss.func(p), 'r+')
