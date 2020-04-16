# -*- coding: utf-8 -*-
"""
Script for various function doing analysis of solutions, maybe also calculation
of exact solution, hm?
"""

import numpy as np
from numpy import newaxis as nax


def get_variation(q, space_sampling, bc, time_sampling=None):
    """
    Computes total variation of q with regard to nodes of space_sampling,
    if time_sampling is provided, assumes q is a matrix and computes:
    variation of each row with regard to space_sampling and X x T variation
    :param q:
    :param space_sampling:
    :param time_sampling:
    :return:
    """
    if time_sampling is None:
        if bc == 'zeroDirichlet':
            Qp = np.append(q[1:], 0)  # update zprava
            Ql = np.append(0, q[:-1])  # update zleva
        elif bc == 'periodic':
            Qp = np.append(q[1:], q[0])  # update zprava
            Ql = np.append(q[-1], q[:-1])  # update zleva
        else:
            Qp = np.append(q[1:], q[-1])  # update zprava
            Ql = np.append(q[0], q[:-1])  # update zleva
        varq = np.sum(np.abs(Qp - Ql))
    else:
        if bc == 'zeroDirichlet':
            Qp = np.hstack((q[:, 1:], np.zeros((np.shape(q)[0], 1))))  # update zprava
            Ql = np.hstack((np.zeros((np.shape(q)[0], 1)), q[:, :-1]))  # update zleva
        elif bc == 'periodic':
            Qp = np.hstack((q[:, 1:], q[:, 0, None]))  # update zprava
            Ql = np.hstack((q[:, -1, None], q[:, :-1]))  # update zleva
        else:
            # NOT TESTED!
            Qp = np.hstack((q[:, 1:], q[:, -1]))  # update zprava
            Ql = np.hstack((q[:, 0], q[:, -1]))  # update zleva
        varq = np.sum(np.abs(Qp - Ql), axis=1)
    return varq


def get_exact_solution(X, T, bc, q0, ux):
    exactY = np.zeros((len(T), len(X)))
    for n in range(len(T)):
        if bc == 'zeroDirichlet':
            exactY[n] = q0(X - ux(X) * (T[n]))
        elif bc == 'periodic':
            positions = X - ux(X) * (T[n])
            x_min = X[0]
            x_max = X[-1]
            x_range = x_max - x_min
            for i in range(len(X)):
                while positions[i] < x_min:
                    positions[i] = positions[i] + x_range
                while positions[i] > x_max:
                    positions[i] = positions[i] - x_range
            exactY[n] = q0(positions)
    # exactY[n] = q0(X - ux(X)*(T[n]%len(X))
    return exactY


def calc_FFT2D(q, X, T):
    # placeholder so far
    return np.fft.rfft2(q)


def calc_FFT1D(y):
    """
    Computes fourier transform ov y, returns only positive part of
    spectrum
    :param y: signal
    :return: frequency bins, intensities
    """
    n = len(y)
    Y = np.fft.fft(y) / n  # fft computing and normalization
    Y = Y[range(n // 2)]
    return np.abs(Y)


def calc_FFT1dfbins(n, Fs):
    """
    :param n: samples
    :param Fs: signal frequency
    :return: frequency bin intervals
    """
    k = np.arange(n)
    T = n / Fs
    frq = k / T  # two sides frequency range
    frq = frq[range(n // 2)]  # one side frequency range
    return frq


def broadcast_sols_1dir(Ys, Ts):
    """
    Broadcasts solutions to
    match the finest time discretization,
    use along with animate1D_dgsol
    :param Ys: list or tuple of solutions
    :param Ts: list or tuple of different time discretizations
    :return: bYs - broadcasted solutions, bTs - broadcasted times
    """
    lens = np.zeros(len(Ts), dtype=np.int64)
    for i, T in enumerate(Ts):
        lens[i] = np.size(T)

    max_len = np.max(lens)
    max_leni = np.where(lens == max_len)[0][0]
    max_T = Ts[max_leni]

    bTs = np.zeros((max_len, len(Ts)))
    bYs = np.zeros((max_len, np.shape(Ys[0])[1], len(Ts)))
    for i, T in enumerate(Ts):
        txt = np.abs(max_T[nax] - T[nax].T)
        tidx = np.argmin(txt, axis=0)
        bTs[:, i] = T[tidx]
        bYs[:, :, i] = Ys[i][tidx, :]

    return bYs, max_T


def cast_sols(Ys, Xs, Ts, fit="min"):
    """
    Broadcasts solutions to
    match the finest time and space discretization,
    use along with animate1D_dgsol
    :param Ys: list or tuple of solutions
    :param Xs: list or tuple of different space discretizations
    :param Ts: list or tuple of different time discretizations
    :return: bYs - broadcasted solutions, bTs - broadcasted times
    """
    Tlens = np.zeros(len(Ts), dtype=np.int64)
    for i, T in enumerate(Ts):
        Tlens[i] = np.size(T)

    if fit == "max":
        target_Tlen = np.max(Tlens)
    else:
        target_Tlen = np.min(Tlens)
    target_leni = np.where(Tlens == target_Tlen)[0][0]
    target_T = Ts[target_leni]

    bTs = np.zeros((target_Tlen, len(Ts)))

    Xlens = np.zeros(len(Xs), dtype=np.int64)
    for i, X in enumerate(Xs):
        Xlens[i] = np.size(X)

    if fit == "max":
        target_Xlen = np.max(Xlens)
    else:
        target_Xlen = np.min(Xlens)
    target_leni = np.where(Xlens == target_Xlen)[0][0]
    target_X = Xs[target_leni]

    bXs = np.zeros((target_Xlen, len(Xs)))

    bYs = np.zeros((target_Tlen, target_Xlen, len(Ts)))
    for i, (X, T) in enumerate(zip(Xs, Ts)):
        txt = np.abs(target_T[nax] - T[nax].T)
        tidx = np.argmin(txt, axis=0)
        bTs[:, i] = T[tidx]

        xxx = np.abs(target_X[nax] - X[nax].T)
        xidx = np.argmin(xxx, axis=0)
        bXs[:, i] = X[xidx]

        Xidx, Tidx = np.meshgrid(xidx, tidx)

        bYs[:, :, i] = Ys[i][Tidx, Xidx]

    return bYs, target_T, target_X
