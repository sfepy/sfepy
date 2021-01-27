from __future__ import print_function
import numpy as nm

from sfepy.base.base import pause, Struct
from sfepy.homogenization.utils import integrate_in_time

def compute_mean_decay(coef):
    r"""
    Compute mean decay approximation of a non-scalar fading memory
    coefficient.
    """
    n_step = coef.shape[0]
    weights = nm.abs(coef[int(nm.fix(n_step / 2))])
    weights /= weights.sum()

    ## coef_flat = nm.reshape(coef, (n_step, nm.prod(coef.shape[1:])))
    ## maxs = nm.abs(coef_flat).max(axis=-1)
    ## decay = avgs + maxs

    aux = weights[None,...] * coef
    aux_flat = nm.reshape(aux, (n_step, nm.prod(coef.shape[1:])))
    avgs = aux_flat.sum(axis=-1)

    decay = avgs # Gives better results than the above.
    decay /= decay[0]

    return decay

def eval_exponential(coefs, x):
    c1, c2 = coefs
    return c1 * nm.exp(-c2 * x)

def approximate_exponential(x, y):
    r"""
    Approximate :math:`y = f(x)` by :math:`y_a = c_1 exp(- c_2 x)`.

    Initial guess is given by assuming y has already the required exponential
    form.
    """
    from scipy.optimize import leastsq

##     weights = nm.abs(y)
##     weights = weights / weights.sum()

    weights = nm.ones_like(y)
    
    def fun(c, x, y):
        val = weights * (y - eval_exponential(c, x))
        return val

    c1 = y[0]
    c2 = - nm.log(y[1] / c1) / (x[1] - x[0])

    coefs, ier = leastsq(fun, nm.array([c1, c2]), args=(x, y))

    if ier != 1:
        print(c1, c2)
        print(coefs, ier)
        pause('exponential fit failed!')

    return coefs

def fit_exponential(x, y, return_coefs=False):
    """
    Evaluate :math:`y = f(x)` after approximating :math:`f` by an exponential.
    """
    coefs = approximate_exponential(x, y)

    if return_coefs:
        return coefs, eval_exponential(coefs, x)
    else:
        return eval_exponential(coefs, x)

class ConvolutionKernel(Struct):
    r"""
    The convolution kernel with exponential synchronous decay approximation
    approximating the original kernel represented by the array :math:`c[i]`,
    :math:`i = 0, 1, \dots`.

    .. math::
        \begin{split}
        & c_0 \equiv c[0] \;, c_{e0} \equiv c_0 c^e_0 \;, \\
        & c(t) \approx c_0 d(t) \approx c_0 e(t) = c_{e0} e_n(t) \;,
        \end{split}
        
    where :math:`d(0) = e_n(0) = 1`, :math:`d` is the synchronous decay and
    :math:`e` its
    exponential approximation, :math:`e = c^e_0 exp(-c^e_1 t)`.
    """
    def __init__(self, name, times, kernel, decay=None,
                 exp_coefs=None, exp_decay=None):
        Struct.__init__(self, name=name, times=times, c=kernel,
                        d=decay, ec=exp_coefs, e=exp_decay)

        if decay is None:
            self.d = compute_mean_decay(kernel)

        if exp_decay is None:
            self.ec, self.e = fit_exponential(times, self.d,
                                              return_coefs=True)

        self.en = self.e / self.e[0]
        self.c0 = self.c[0]
        self.e_c0 = self.c0 * self.e[0]
        self.e_d1 = self.e[1] / self.e[0]

        self.c_slice = (slice(None),) + ((None,) * (self.c.ndim - 1))

    def diff_dt(self, use_exp=False):
        """
        The derivative of the kernel w.r.t. time.
        """
        if use_exp:
            val = - self.ec[1] * self.c0 * self.e[self.c_slice]

        else:
            zz = nm.zeros(self.c[0:1].shape, dtype=self.c.dtype)
            val = nm.r_[nm.diff(self.c, axis=0) /
                        nm.diff(self.times, axis=0)[self.c_slice], zz]

        return val

    def int_dt(self, use_exp=False):
        """
        The integral of the kernel in time.
        """
        if use_exp:
            val = (self.e_c0/self.ec[1]) * (1.0 - self.en[-1])
        else:
            val = integrate_in_time(self.c, self)

        return val

    def get_exp(self):
        """
        Get the exponential synchronous decay kernel approximation.
        """
        return self.c0 * self.e[self.c_slice]

    def get_full(self):
        """
        Get the original (full) kernel.
        """
        return self.c

    def __call__(self, use_exp=False):
        """
        Get the kernel or its approximation.
        """
        if use_exp:
            return self.get_exp()

        else:
            return self.get_full()
