# -*- coding: utf-8 -*-
"""
This module contains solver and various updaters for nonlinear
advection equation

q_t + [f(q)]_x = 0

Updaters corespond to different schemes which solver uses
to calculate solution between time levels.

Calling scheme is:

Y, X, T = some_solver(some_updater, [data: initial condition, enviroment variables], a, b, t_0, t_f
n_space_steps, n_time_steps)

All data are expected to be functions from which solver creates samples.
"""
import warnings
import numpy as np
apd = np.append
from matplotlib.ticker import MaxNLocator  # do not remove this import!
import time


#------------#
#   Solver   #
#------------#
def nonlin_solver(nlin_update, q_0x, fq, df, a, b, t0, tf, bc, space_steps_J=100, time_steps_N=None, c=1.0):
	"""
		Solver for 1D nonlinear advection equation
			q_t + [f(q)]_x = 0,
		best used with one of our updaters.
		Warning using with third party updaters may cause unstability.

		:param nlin_update: method formula
		:param q_0x: initial condition
		:param fq: function f, should be at least C^1([a,b])
		:param df: derivation f' of f
		:param a: spatial interval start
		:param b: spatial interval end
		:param t0: time interval start
		:param tf: time interval end
		:param space_steps_J: number of uniform steps in space variable, time stepping is determined to fulfill stability
				condition |max{f(q0(x))} * dt/dx| = c (default c=1)
		:param time_steps_N: number of time steps, if you don't want them to be calculated for you
		:param c: desired generalized "courant" number, see above, default is 1, if time_steps_N provided this is ignored
		:param bc: type of boundary condition
		:return: Y (solution), X(space interval discretization array), T (time interval discretization array)
		"""

	# TODO how do we get stability conditions for nonlinear case?
	space_sampling, dx = np.linspace(a, b, space_steps_J + 1, retstep=True)
	u = df(q_0x(space_sampling))  # to get courant number
	maxu = max(abs(u))
	if time_steps_N is None:
		dt = dx / maxu * c
		# time_steps_N = int((tf - t0) / dt) * 2
		time_steps_N = int(np.ceil((tf - t0) / dt))
	time_sampling, dt = np.linspace(t0, tf, time_steps_N + 1, retstep=True)
	dtdx = dt / dx

	print("--- Update: "+str(nlin_update.__name__) + " ---")
	print("Space divided into {0} nodes, {1} steps, step size is {2}".format(space_steps_J + 1, space_steps_J, dx))
	print("Time divided into {0} nodes, {1} steps, step size is {2}".format(time_steps_N + 1, time_steps_N, dt))
	print("Time Level 0 Courant number c = max(abs(f'(q_0))) * dt/dx = {0}".format(maxu * dtdx))

	# print("\nDiffusion parameter alph = diff*dt/dx^2 = {0}".format(diff * dtdx / dx))
	# with warnings.catch_warnings():
	# 	warnings.simplefilter("ignore")
	# 	print("Péclet number P_d = {0} \n".format(maxu / diff * dx))
	#
	# print("alph + c/4 = {0} <= 1/2".format(diff * dtdx / dx + maxu*dtdx/4))
	# print("c^2 = {0} <= 2 alph = {1}".format((maxu * dtdx) ** 2, diff * dtdx/dx))

	q = np.zeros((len(time_sampling), len(space_sampling)))
	q[0] = q_0x(space_sampling)  # sampling initial values
	q_R = np.append(q[0][1:], q[0][-1])
	vertices = (q[0] + q_R) / 2  # initial vertices - middle values of q[0]

	for n in range(0, len(time_sampling) - 1):
		q[n + 1], vertices = nlin_update(q[n], fq, dtdx, bc, dx, vertices)

	return q, space_sampling, time_sampling


#---------------#
#   Updaters    #
#---------------#
def laxwen_nonlin_1d(q, f, dtdx, bc, dx, vertices):
	"""
	Updater for equation
		q_t + [f(q)]_x = 0
	Uses Lax-Wendroff scheme

	:param q: solution in previous time step
	:param f: sampled velocity
	:param dtdx: time step to space step ratio
	:return: solution in current time step
	"""
	fq = f(q)
	if bc == 'zeroDirichlet':
		q_L = np.append(0, q[:-1])
		q_R = np.append(q[1:], 0)
	elif bc == 'periodic':
		q_L = np.append(q[-1], q[:-1])
		q_R = np.append(q[1:], q[0])
	else:
		q_L = np.append(q[0], q[:-1])
		q_R = np.append(q[1:], q[-1])

	ap = (q_R - q)
	ap[ap != 0] = (f(q_R) - fq)[ap != 0] / ap[ap != 0]
	Fp = 1 / 2 * (fq + f(q_R)) - 1 / 2 * ap**2 * dtdx*(q_R - q)
	al = (q - q_L)
	al[al != 0] = (fq - f(q_L))[al != 0] / al[al != 0]
	Fl = 1 / 2 * (f(q_L) + fq) - 1 / 2 * al**2 * dtdx*(q - q_L)

	res = q - dtdx * (Fp - Fl)
	return res, None


def laxfrie_nonlin_1d(q, f, dtdx, bc, dx, vertices):
	"""
	Updater for equation
		q_t + [u(x)q]_x = 0
	Uses Lax-Friedrichs scheme

	:param q: solution in previous time step
	:param f:
	:param dtdx: time step to space step ratio
	:return: solution in current time step
	"""

	if bc == 'zeroDirichlet':
		avgq = np.append(0, np.append((q[:-2] + q[2:]) / 2, 0)) #FIXME why q[:-2] + q[2:] not q[:-1] + q[1:] ??
		q_L = np.append(0, q[:-1])
		q_R = np.append(q[1:], 0)
	elif bc == 'periodic':
		avgq = np.append((q[-1] + q[1]) / 2, np.append((q[:-2] + q[2:]) / 2, (q[-2] + q[0]) / 2))
		q_L = np.append(q[-1], q[:-1])
		q_R = np.append(q[1:], q[0])
	else:
		avgq = np.append(0, np.append((q[:-2] + q[2:]) / 2, 0))
		q_L = np.append(q[0], q[:-1])
		q_R = np.append(q[1:], q[-1])

	res = avgq - dtdx / 2 * (f(q_R) - f(q_L))
	return res, None


def quickest_nonlin_1d(q, f, dtdx, bc, dx, vertices):
	"""

	QUICK with Estimated Streaming Terms for nonlinear advection problem:
		q_t + [f(q)]_x = 0

	After some wild generalizations made ad hoc, based on
	Leonard, B. P. (1979). A Stable And Accurate Convective Modelling Procedure
		Based on Quadraticu pstream Interpolation.
	Computer Methods in Applied Mechanics and Engineering, 19, 59–98.

	Stability conditions can be found in Fig.18 p.82 what is stability condition for nonlinear case, hm?

	:param q:
	:param f:
	:param dtdx:
	:param dx:
	:return:
	"""
	diff = 0
	if bc == 'zeroDirichlet':
		q_L = np.append(0, q[:-1])
		q_R = np.append(q[1:], 0)
		q_FL = np.append(0, np.append(0, q[:-2]))
	elif bc == 'periodic':
		q_L = np.append(q[-1], q[:-1])
		q_R = np.append(q[1:], q[0])
		q_FL = np.append(q[-2:], q[:-2])
	else:
		q_L = np.append(q[0], q[:-1])
		q_R = np.append(q[1:], q[-1])
		q_FL = np.append(q[:2], q[:-2])

	grad_l = (f(q) - f(q_L))/dx
	grad_l2 = (f(q)**2 - f(q_L)**2) / dx
	curv_l = (f(q_FL) + f(q) - 2*f(q_L))/dx**2  # [(52), p.79]
	curv_l3 = (f(q_FL)**3 + f(q)**3 - 2 * f(q_L)**3) / dx ** 2  # [(52), p.79]

	grad_r = (f(q_R) - f(q))/dx
	grad_r2 = (f(q_R)**2 - f(q)**2) / dx
	curv_r = (f(q_L) + f(q_R) - 2*f(q))/dx**2  # [(53), p.79]
	curv_r3 = (f(q_L)**3 + f(q_R)**3 - 2 * f(q)**3) / dx ** 2  # [(53), p.79]

	# diff is diffusion coefficient currently set to zero
	# TODO this seems to work - why?
	res = q - dtdx * ((1/2*(f(q) + f(q_R)) - dx/2*dtdx*grad_r2 - dx**2/6*(curv_r-dtdx**2*curv_r3-3*diff*curv_r))
				    - (1/2*(f(q_L) + f(q)) - dx/2*dtdx*grad_l2 - dx**2/6*(curv_l-dtdx**2*curv_l3-3*diff*curv_l))) \
		    + diff*((dx*grad_r-dx**2*dtdx*curv_r) - (dx*grad_l - dx**2/2*dtdx*curv_l))  # [(57), p.80]

	return res, None


def naiveupwind_nonlin_1d(q, f, dtdx, bc, dx, vertices):
	"""
	First order iterative upwind (prednasky SNM2)
	method for nonlinear equation
		q_t + [f(q)]_x = 0

	This is kinda naive version (hence the name) but works pretty good, at least for burgess :-)

	:param q: solution in prevorius time
	:param f: sampled velocity
	:param dtdx: dt/dx
	:return: solution in current time step
	"""

	fq = f(q)
	if min(fq) >= 0:
		res = q[1:] - dtdx * (fq[1:] - fq[:-1])
		if bc == 'zeroDirichlet':
			res = np.append(0, res) # zero boundary condition on left side
		elif bc == 'periodic':
			res = np.append(q[-1], res) # TODO are these right boundaries? I think so.
	else:
		res = q[:-1] - dtdx * (fq[1:] - fq[:-1])
		if bc == 'zeroDirichlet':
			res = np.append(res, 0) # zero boundary condition on right side
		elif bc == 'periodic':
			res = np.append(res, q[0])  # TODO are these right boundaries? I think so.
	return res, None


def upwind_nonlin_1d(q, f, dtdx, bc, dx, vertices):
	"""
	First order iterative "upwind" (prednasky SNM2)
	method for nonlinear equation
		q_t + [f(q)]_x = 0

	:param q: solution in prevorius time
	:param f: sampled velocity
	:param dtdx: dt/dx
	:return: solution in current time step
	"""

	fq = f(q)
	if bc == 'zeroDirichlet':
		q_L = np.append(0, q[:-1])
		q_R = np.append(q[1:], 0)
		q_FL = np.append(0, np.append(0, q[:-2]))
	elif bc == 'periodic':
		q_L = np.append(q[-1], q[:-1])
		q_R = np.append(q[1:], q[0])
		q_FL = np.append(q[-2:], q[:-2])
	else:
		q_L = np.append(q[0], q[:-1])
		q_R = np.append(q[1:], q[-1])
		q_FL = np.append(q[:2], q[:-2])

	ap = (q_R - q)
	ap[ap != 0] = (f(q_R) - fq)[ap != 0] / ap[ap != 0]
	Fp = 1/2 * (fq + f(q_R)) - 1/2*np.abs(ap)*(q_R - q)

	al = (q - q_L)
	al[al != 0] = (fq - f(q_L))[al != 0]/al[al != 0]
	Fl = 1/2 * (f(q_L) + fq) - 1/2*np.abs(al)*(q - q_L)

	res = q - dtdx * (Fp - Fl)
	return res, None


def harten_nonlin_1d(q, f, dtdx, bc, dx, vertices):
	"""
	A. Harten hybrid method (prednasky SNM2)
	for nonlinear equation
		q_t + (f(q))_x = 0

	It does not work yet

	:param q: solution in prevorius time
	:param f: sampled velocity
	:param dtdx: dt/dx
	:return: solution in current time step
	"""

	fq = f(q)
	if bc == 'zeroDirichlet':
		q_L = np.append(0, q[:-1])
		q_R = np.append(q[1:], 0)
		q_FL = np.append(0, np.append(0, q[:-2]))
	elif bc == 'periodic':
		q_L = np.append(q[-1], q[:-1])
		q_R = np.append(q[1:], q[0])
		q_FL = np.append(q[-2:], q[:-2])
	else:
		q_L = np.append(q[0], q[:-1])
		q_R = np.append(q[1:], q[-1])
		q_FL = np.append(q[:2], q[:-2])

	ap = (q_R - q)
	ap[ap != 0] = (f(q_R) - fq)[ap != 0] / ap[ap != 0]
	FpH = 1 / 2 * (fq + f(q_R)) - 1 / 2 * ap**2 * dtdx*(q_R - q)
	FpL = 1 / 2 * (fq + f(q_R)) - 1/2*np.abs(ap)*(q_R - q)
	thp = (q - q_L) / (q_R - q)  # FIXME get theta calculation for nonlinear case

	Fp = FpL + fluxlim1(thp) * (FpH - FpL)

	al = (q - q_L)
	al[al != 0] = (fq - f(q_L))[al != 0] / al[al != 0]
	FlH = 1 / 2 * (f(q_L) + fq) - 1 / 2 * al**2 * dtdx*(q - q_L)
	FlL = 1 / 2 * (f(q_L) + fq) - 1 / 2 * np.abs(al) * (q - q_L)
	thl = (q_L - q_FL)/(q - q_L)

	Fl = FpL + fluxlim1(thl) * (FlH - FlL)

	res = q - dtdx * (Fp - Fl)
	return res, None


def naive_afs_nonlin_1d(q, f, dtdx, bc, dx, vertices):
	"""
	AFS

	Eymann, Timothy Andrew (2013). Active Flux Schemes.
	Dissertation in the University of Michigan (Aerospace Engineering and Scientific Computing)

	Scheme p. 41, 46, 47

	:param q:
	:param u:
	:param dtdx:
	:param vertices:
	:return:
	"""

	if bc == 'zeroDirichlet':
		q_L_mid = np.append(0, vertices[:-1])
		q_R_mid = vertices
		q_L = np.append(0, q[:-1])
		q_LL_mid = np.append(0, np.append(0, vertices[:-2]))
		q_R = np.append(q[1:], 0)
		q_RR_mid = np.append(vertices[1:], 0)
	elif bc == 'periodic':
		q_L_mid = np.append(vertices[-1], vertices[:-1])
		q_R_mid = vertices
		q_L = np.append(q[-1], q[:-1])
		q_LL_mid = np.append(vertices[-2:], vertices[:-2])
		q_R = np.append(q[1:], q[0])
		q_RR_mid = np.append(vertices[1:], vertices[0])
	else:
		q_L_mid = np.append(vertices[0], vertices[:-1])
		q_R_mid = vertices
		q_L = np.append(q[0], q[:-1])
		q_LL_mid = np.append(vertices[:2], vertices[:-2])
		q_R = np.append(q[1:], q[-1])
		q_RR_mid = np.append(vertices[1:], vertices[-1])

	# c1 = q_L_mid
	# c2 = 1/4*(6*q - q_L_mid - q_R_mid)
	# c3 = q_R_mid
	# alpha = 2*(c1-2*c2+c3)
	# beta = -(3*c1-4*c2+c3)
	# def phi1(ksi):
	# 	return (2*ksi-1)*(ksi-1)
	# def phi2(ksi):
	# 	return 4*ksi*(1-ksi)
	# def phi3(ksi):
	# 	return ksi*(2*ksi-1)

	res = np.zeros(len(q))
	u = 1 / 2 * (q_L_mid + q_R_mid)
	for i in range(len(q)):
		if u[i] >= 0:
			ni = u[i] * dtdx
			res[i] = ni ** 2 * (ni - 1) * q_LL_mid[i] + ni ** 2 * (3 - 2 * ni) * q_L[i] + ni * (1 - ni) * q_L_mid[i] + (
						1 - ni) ** 2 * (1 + 2 * ni) * q[i] - ni * (1 - ni) ** 2 * q_R_mid[i]
			vertices[i] = ni * (3 * ni - 2) * q_L_mid[i] + 6 * ni * (1 - ni) * q[i] + (1 - ni) * (1 - 3 * ni) * q_R_mid[i]

			# ksi_0_0 = 2 * (1 - c1 * dtdx) / ((beta * dtdx + 1) - np.sqrt((beta * dtdx + 1) ** 2 + 4 * alpha * dtdx * (1 - c1 * dtdx)))
			# ksi_0_1 = 2 * (1 - c1 * dtdx) / ((beta * dtdx + 1) + np.sqrt((beta * dtdx + 1) ** 2 + 4 * alpha * dtdx * (1 - c1 * dtdx)))
			# vertex_0 = c1*phi1(ksi_0_0) + c2*phi2(ksi_0_0) + c3*phi3(ksi_0_0)
			# vertex_1 = c1*phi1(ksi_0_1) + c2*phi2(ksi_0_1) + c3*phi3(ksi_0_1)
		else:
			ni = u[i] * dtdx
			res[i] = - ni ** 2 * (ni + 1) * q_RR_mid[i] + ni ** 2 * (3 + 2 * ni) * q_R[i] - ni * (1 + ni) * q_R_mid[i] + (
						1 + ni) ** 2 * (1 - 2 * ni) * q[i] + ni * (1 + ni) ** 2 * q_L_mid[i]
			vertices[i] = ni * (3 * ni + 2) * q_R_mid[i] - 6 * ni * (1 + ni) * q[i] + (1 + ni) * (1 + 3 * ni) * q_L_mid[i]

	return res, vertices

def afs_nonlin_1d(q, f, dtdx, bc, dx, vertices):
	"""
	AFS

	Eymann, Timothy A.; Roey Phillip L.. Active Flux Schemes for Systems. 20th AIAA Computational Fluid Dynamics Conference, Honolulu, Hawaii, June 2011.
	Scheme p. 3, 4

	Eymann, Timothy Andrew (2013). Active Flux Schemes. Dissertation in the University of Michigan (Aerospace Engineering and Scientific Computing)
	Scheme p. 41, 46, 47)

	:param q:
	:param u:
	:param dtdx:
	:param vertices:
	:return:
	"""

	if bc == 'zeroDirichlet':
		q_L_mid = np.append(0, vertices[:-1])
		q_R_mid = vertices
	elif bc == 'periodic':
		q_L_mid = np.append(vertices[-1], vertices[:-1])
		q_R_mid = vertices
	else:
		q_L_mid = np.append(vertices[0], vertices[:-1])
		q_R_mid = vertices

	res = np.zeros(len(q))
	mid_vertices = np.zeros(len(q))
	u = 1 / 2 * (q_L_mid + q_R_mid)
	for i in range(len(q)):
		if u[i] >= 0:
			ni = u[i] * dtdx
			mid_vertices[i] = q_R_mid[i] - ni*((ni-1)*(q[i]-q_L_mid[i])+(2-ni)*(q_R_mid[i]-q[i]))
			res[i] = q[i] - dtdx*(f(mid_vertices[i]) - f(mid_vertices[i-1]))
			vertices[i] = q_R_mid[i] - ni*((3*ni-2)*(q[i]-q_L_mid[i])+(4-3*ni)*(q_R_mid[i]-q[i]))
		else:
			ni = u[i] * dtdx
			mid_vertices[i] = q_R_mid[i] + ni*((2 + ni) * (q[i] - q_L_mid[i]) + (-ni - 1) * (q_R_mid[i] - q[i]))
			res[i] = q[i] - dtdx * (f(mid_vertices[i]) - f(mid_vertices[i - 1]))
			vertices[i] = q_R_mid[i] + ni * ((4 + 3 * ni) * (q[i] - q_L_mid[i]) - (3 * ni + 2) * (q_R_mid[i] - q[i]))
	return res, vertices


#---------------#
#   Limiters    #
#---------------#
# so called flux limiter
def fluxlim1(th):
	return np.maximum(np.minimum(np.ones(np.size(th)), 2*th),
                      np.minimum(np.ones(np.size(th)), 2*th),
                      np.zeros(np.size(th)))






