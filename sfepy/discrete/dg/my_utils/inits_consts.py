# -*- coding: utf-8 -*-
"""
Different initial conditions and environment constants as
functions to enable different samplig in solvers
"""
import numpy as np
from scipy import signal

#---------------------------------------#
#   Some useful piecewise functions     #
#---------------------------------------#
a = 0
b = 1


def left_par_q(x):
    """
    q(x) =  -(x - a) * (x - (a + 1)) for  a < x <a+1
             0 elswhere
    i.e. concave parabola with roots a and a + 1
    :param x:
    :return:
    """
    return np.piecewise(x, [x < a, 0.1 <= x , 0.3 < x],
                           [0, lambda t: -100*(t - 0.1) * (t - 0.3), 0])


def right_par_q(x):
    """
    q(x) =  -(x - b) * (x - (b - 2)) for  b - 2  < x < b
             0 elswhere
    i.e. concave parabola with roots b - 2 and b
    :param x:
    :return:
    """
    return np.piecewise(x, [x < b - .2, b - .2 <= x,  b < x],
                           [0, lambda t: -100*(t - b) * (t - (b - .2)), 0])


def middle_par_q(x):
    return np.piecewise(x, [x <= a, x <= a + .3, a + .3 < x, a + .4 <= x],
                           [0, 0, lambda t: -100*(t - (a + .4)) * (t - (a + .3)), 0])


def left_cos(x):
    """

    :param x:
    :return:
    """
    return np.piecewise(x, [x <= a, x <= a + .3, a + .3 < x, a + .4 <= x],
                        [0, 0, lambda t: (np.cos(np.pi * 20 * (t - .35)) + 1) / 2, 0])


#----------------------------------#
#   Piecewise constant functions   #
#----------------------------------#
def three_step_q(x):
    """
    piecewise constant (-inf, a],(a, a + 5](a+5, inf)
    :param x:
    :return:
    """
    return np.piecewise(x, [x <= a, x <= a + .5, a + .5 < x], [0, 1, 0])


def three_step_u(x):
    """
    piecewise constant (-inf, a],(a, a + 5](a+5, inf)

    :param x:
    :return:
    """
    return np.piecewise(x, [x <= 0, x <= .7, .7 < x], [0, .5, 2])


def four_step_u(x):
    """
    piecewise constant (-inf, 1.8],(1.8, a + 4](a+4, a + 5](a + 5, inf)

    :param x:
    :return:
    """
    return np.piecewise(x, [x <= a, x <= a + .4, a + .4 < x, a + .5 <= x], [0, 1, 1.5, 1])


def four_step_q(x):
    """
    piecewise constant (-inf, 1.8],(1.8, a + 4](a+4, a + 5](a + 5, inf)

    :param x:
    :return:
    """
    return np.piecewise(x, [x <= a, x <= a + .4, a + .4 < x, a + .5 <= x], [0, 0, 1, 0])

#------------------------#
#   Constant functions   #
#------------------------#
def const_u(x):
    """
    piecewise constant (-inf, a],(a, a + 5](a+5, inf)
    :param x:
    :return:
    """
    return np.ones((np.size(x), 1)) * 0


def const_q(x):
    """
    piecewise constant (-inf, a],(a, a + 5](a+5, inf)
    :param x:
    :return:
    """
    return np.ones(np.size(x)) * .5


#------------------------------#
#   Bell shaped initial cond   #
#------------------------------#
def ghump(x):
    """
    :param x:
    :return:
    """
    return np.exp(-200 * x**2)


def gauss_init(x):
    """

    :param x:
    :return:
    """
    return np.array(ghump(x-.2))

def gsmooth(x):
    """
    :param x:
    :return:
    """
    return np.piecewise(x, [x[:, 0] <= - 1, x[:, 0] >= -1, 1 < x[:, 0]], [0, lambda x: np.exp(1/(x**2 - 1)), 0])


#----------------------------------#
#    Piecewise linear function     #
#----------------------------------#
def sawtooth_q(x):
    """
    piecewise linear
    :param x:
    :return:
    """
    return signal.sawtooth(2 * np.pi * 5 * x)

#----------------------------------#
#        Sinus and constant        #
#----------------------------------#
def cos_const_q(x):
    """
    :param x:
    :return:
    """
    return np.piecewise(x, [x <= a + .3, a + .3 < x, a + .35 < x, a + .45 < x, a + .5 <= x],
                        [0, lambda t: (np.cos(np.pi * 20 * (t - .35)) + 1) / 2, 1,
                            lambda t: (np.cos(np.pi * 20 * (t - .45)) + 1) / 2, 0])

#----------------------------------#
#   Quadratic and cubic function   #
#----------------------------------#
def quadr_cub(x):
    """
    :param x:
    :return:
    """
    return np.piecewise(x, [x <= a + .1, a + .1 < x, a + .3 < x, a + .4 <= x],
                        [0, lambda t: -25*(t - 0.1) * (t - 0.5),
                            lambda t: -250*(t - 0.1) * (t - 0.1) * (t-0.4), 0])


#-------------------------#
#   System initial cond   #
#-------------------------#
def qsysinit(x):
    """
    initial values for nonconservative system
    :param x:
    :return:
    """
    return np.array([ghump(x),
                     np.zeros(len(x))]).swapaxes(0, 1)


def wsysinit(x):
    """
    transformed innitial value for noncnons system solver
    :param x:
    :return:
    """
    return np.stack((np.sum(invRfunc(x)[:, 0, :] * gauss_init(x), 1),
                     np.sum(invRfunc(x)[:, 1, :] * gauss_init(x), 1)), 1)


#---------------------------------------#
#    System parameters from article     #
#---------------------------------------#
arta = -5.
artb = 3.


def Kfunc(x):
    """
    bulk modulus function
    :param x:
    :return:
    """
    return np.piecewise(x, [x <= 0, x <= 5, 5 < x], [1, 1, 1])


def rhofunc(x):
    """
    density function
    :param x:
    :return:
    """
    return np.piecewise(x, [x <= 0, x > 0], [1, 4])


def cfunc(x):
    """
    sound velocity function
    :param x:
    :return:
    """
    return np.sqrt(Kfunc(x) / rhofunc(x))


def Zfunc(x):
    """
    impendance function
    :param x:
    :return:
    """
    return np.sqrt(Kfunc(x) * rhofunc(x))


def Afunc(x):
    """
    enviroment matrix
    A(x) =  [   0    , K(x) ]
            [1/rho(x),  0   ]
    :param x:
    :return:
    """
    return np.array([[np.zeros(len(x)), 1 / rhofunc(x)],
                     [Kfunc(x), np.zeros(len(x))]]).swapaxes(0, 2)


def ATfunc(x):
    """
    transposed enviroment matrix,
    for adjoint conservative system
    :param x:
    :return:
    """
    return np.transpose(Afunc(x), (0, 2, 1))


def eigAfunc(x):
    """
    Jordan shape of A
    :param x:
    :return:
    """
    return np.array([[-cfunc(x), np.zeros(len(x))],
                     [np.zeros(len(x)), cfunc(x)]]).swapaxes(0, 2)


def Rfunc(x):
    """
    eigen vectors of A
    :param x:
    :return:
    """
    return np.array([[-Zfunc(x), np.ones(len(x))],
                     [Zfunc(x), np.ones(len(x))]]).swapaxes(0, 2)


def invRfunc(x):
    """
    inverse egigen vectors of A
    :param x:
    :return:
    """
    return np.array([[-1 / (2 * Zfunc(x)), np.ones(len(x)) / 2],
                     [-1 / (2 * Zfunc(x)), np.ones(len(x)) / 2]]).swapaxes(0, 2)
