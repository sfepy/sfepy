# -*- coding: utf-8 -*-
"""
This module contains solver and various updaters for conservatve and nonconservative
advection equation. Updaters corespond to different schemes which solver uses
to calculate solution between time levels.

Calling scheme is:

Y, X, T = some_solver(some_updater, [data: initial condition, enviroment variables], a, b, t_0, t_f
n_space_steps, n_time_steps)

All data are expected to be functions from which solver creates samples.
"""
import warnings

import numpy as np
from matplotlib.ticker import MaxNLocator  # do not remove this import!
import time


# diffusion parameter for QUICK and QUISKEST
# FIXME make into parameter of solver and updaters?
diff = 0


#------------#
#   Solver   #
#------------#
def lin_solver(update, q_0x, ux, a, b, t0, tf, bc, space_steps_J=100, c=1.):
	"""
	Solver for 1D hyperbolic equation, which uses update
	method

	:param update: method formula
	:param q_0x: initial condition
	:param ux: velocity function
	:param a: spatial interval start
	:param b: spatial interval end
	:param t0: time interval start
	:param tf: time interval end
	:param space_steps_J: number of uniform steps in space variable, time stepping is determined to fulfill stability
			condition |max{ux} * dt/dx| = c (default c=1)
	:param c: desired courant number, see above, default is 1
	:param bc: type of boundary condition
	:return: Y (solution), X(space interval discretization array), T (time interval discretization array)
	"""

	space_sampling, dx = np.linspace(a, b, space_steps_J + 1, retstep=True)
	u = ux(space_sampling)  # sampled velocity data
	maxu = max(abs(u))
	dt = dx / maxu * c
	# time_steps_N = int((tf - t0) / dt) * 2
	time_steps_N = int(np.ceil((tf - t0) / dt))
	time_sampling, dt = np.linspace(t0, tf, time_steps_N + 1, retstep=True)
	dtdx = dt/dx
	print("--- Update: "+str(update.__name__) + " ---")
	print("Space divided into {0} nodes, {1} steps, step size is {2}".format(space_steps_J + 1, space_steps_J, dx))
	print("Time divided into {0} nodes, {1} steps, step size is {2}".format(time_steps_N + 1, time_steps_N, dt))
	print("Courant number c = max(abs(u)) * dt/dx = {0}".format(maxu*dtdx))
	print("Diffusion parameter alph = diff*dt/dx^2 = {0}\n".format(diff * dtdx/dx))
	# with warnings.catch_warnings():
	# 	warnings.simplefilter("ignore")
	# 	print("Péclet number P_d = {0} \n".format(maxu / diff * dx))
	#
	# print("alph + c/4 = {0} <= 1/2".format(diff * dtdx / dx + maxu*dtdx/4))
	# print("c^2 = {0} <= 2 alph = {1}".format((maxu * dtdx) ** 2, diff * dtdx/dx))

	q = np.zeros((len(time_sampling), len(space_sampling)))
	q[0] = q_0x(space_sampling)  # sampling initial values
	q_R = np.append(q[0][1:], q[0][-1])
	vertices = (q[0]+q_R)/2  # initial vertices - middle values of q[0]

	for n in range(0, len(time_sampling)-1):
		q[n+1], vertices = update(q[n], u, dtdx, bc, dx, vertices)

	return q, space_sampling, time_sampling


#---------------#
#   Updaters    #
#---------------#
# Stability condition for central is different!!
def central_cons_1d(q, u, dtdx, bc, dx, vertices):
	"""
	Updater for equation
		q_t + (u(x)q)_x = 0 i.e. q_t + u(x)q_x + u(x)_xq = 0
	Uses central difference scheme

	:param q: solution in prevorius time
	:param u: sampled velocity
	:param dtdx: dt/dx
	:return:solution in current time step
	"""
	if bc == 'zeroDirichlet':
		Fp = np.append(u[1:] * q[1:], 0)  # update zprava
		Fl = np.append(0, u[:-1] * q[:-1])  # update zleva
	elif bc == 'periodic':
		Fp = np.append(u[1:] * q[1:], u[0]*q[0])  # update zprava
		Fl = np.append(u[-1]*q[-1], u[:-1] * q[:-1])  # update zleva
	else:
		Fp = np.append(u[1:] * q[1:], q[-1])  # update zprava
		Fl = np.append(q[0], u[:-1] * q[:-1])  # update zleva
	res = q - dtdx * (Fp - Fl)/2
	return res, None


def laxfrie_cons_1d(q, u, dtdx, bc, dx, vertices):
	"""
	Updater for equation
		q_t + (u(x)q)_x = 0 i.e. q_t + u(x)q_x + u(x)_xq = 0
	Uses Lax-Frie scheme

	:param q: solution in previous time step
	:param u: sampled velocity
	:param dtdx: time step to space step ratio
	:return: solution in current time step
	"""
	
	if bc == 'zeroDirichlet':
		avgq = np.append(0, np.append((q[:-2] + q[2:]) / 2, 0))
		Fp = np.append(u[1:] * q[1:], 0)  # update zprava
		Fl = np.append(0, u[:-1] * q[:-1])  # update zleva
	elif bc == 'periodic':
		avgq = np.append((q[-1] + q[1]) / 2, np.append((q[:-2] + q[2:]) / 2,(q[-2] + q[0]) / 2))
		Fp = np.append(u[1:] * q[1:], u[0]*q[0])  # update zprava
		Fl = np.append(u[-1]*q[-1], u[:-1] * q[:-1])  # update zleva
	res = avgq - dtdx/2 * (Fp - Fl)
	return res, None


def laxwen_cons_1d(q, u, dtdx, bc, dx, vertices):
	"""
	Updater for equation
		q_t + (u(x)q)_x = 0 i.e. q_t + u(x)q_x + u(x)_xq = 0
	Uses Lax-Wndroff scheme

	:param q: solution in previous time step
	:param u: sampled velocity
	:param dtdx: time step to space step ratio
	:return: solution in current time step
	"""

	if bc == 'zeroDirichlet':
		Fp = np.append(u[1:] * q[1:], 0)  # update zprava
		Fl = np.append(0, u[:-1] * q[:-1])  # update zleva
	elif bc == 'periodic':
		Fp = np.append(u[1:] * q[1:], u[0] * q[0])  # update zprava
		Fl = np.append(u[-1] * q[-1], u[:-1] * q[:-1])  # update zleva

	res = q - dtdx/2 * (Fp - Fl) + u*dtdx**2/2*(Fp + Fl - u*2*q)

	return res, None


def quick_noncons_1d(q, u, dtdx, bc, dx, vertices):
	"""

	Quadratic Upstream Interpolation for Convective Kinematics (QUICK) from article

	Leonard, B. P. (1979). A Stable And Accurate Convective Modelling Procedure
		Based on Quadraticu pstream Interpolation.
	Computer Methods in Applied Mechanics and Engineering, 19, 59–98.

	Only theoretical method, has unreasonably strict stability conditions:
	alph + c/4 <1/2  [(38), p. 74]
	c^2 < 2*alph  [(39), p. 74]

	:param q:
	:param u:
	:param dtdx:
	:param dx:
	:return:
	"""
	
	if bc == 'zeroDirichlet':
		q_L = np.append(0,  q[:-1])
		q_R = np.append(q[1:], 0)
		q_FL = np.append(0, np.append(0, q[:-2]))
		q_FR = np.append(np.append(q[2:], 0), 0)
	elif bc == 'periodic':
		q_L = np.append(q[-1], q[:-1])
		q_R = np.append(q[1:], q[0])
		q_FL = np.append(q[-2:], q[:-2])
		q_FR = np.append(q[2:], q[:2])
	else:
		q_L = np.append(q[0], q[:-1])
		q_R = np.append(q[1:], q[-1])
		q_FL = np.append(q[:2], q[:-2])
		q_FR = np.append(q[2:], q[-2:])

	if min(u) >= 0:
		# q_r =  1/2(q_C + q_R) - 1/8(q_L + q_R - 2q_C)  # update zprava [(23), p.71]
		# q_l = 1/2(q_L + q_C) - 1/8(q_FL + q_C - 2q_L)  # update zleva  [(25), p.71]
		# res = q + [u_l * q_l - u_r * q_r + Gamma_r*(dq/dx)_r - Gamma_l*(dq/dx)_l]*dt/dx [(3), p.61]
		Fr = (q + q_R)/2 - (q_L - 2*q + q_R)/8  # update zprava
		Fl = (q_L + q)/2 - (q_FL - 2*q_L + q)/8  # update zleva
		res = q - dtdx * u * (Fr - Fl) + dtdx * diff*(q_L + q_R - 2*q)/dx
	else:
		Fr = (q + q_R) / 2 - (q - 2 * q_R + q_FR) / 8  # update zprava
		Fl = (q_L + q) / 2 - (q_L - 2*q + q_R) / 8  # update zleva
		res = q - dtdx * u * (Fr - Fl) + dtdx * diff*(q_L + q_R - 2*q)/dx
		pass

	return res, None


def leith_noncons_1d(q, u, dtdx, bc, dx, vertices):
	"""
	Based on

	Leonard, B. P. (1979). A Stable And Accurate Convective Modelling Procedure
		Based on Quadraticu pstream Interpolation.
	Computer Methods in Applied Mechanics and Engineering, 19, 59–98.

	Just really basic streaming estimation procedure for comparison with QUISKEST

	Stability condition is given by
	0 <= alpha <= (1-c**2)/2

	:param q:
	:param u:
	:param dtdx:
	:return:
	"""

	if bc == 'zeroDirichlet':
		q_L = np.append(0,  q[:-1])
		q_R = np.append(q[1:], 0)
		Fp = np.append(u[1:] * q[1:], 0)  # update zprava
		Fl = np.append(0, u[:-1] * q[:-1])  # update zleva
	elif bc == 'periodic':
		q_L = np.append(q[-1], q[:-1])
		q_R = np.append(q[1:], q[0])
		Fp = np.append(u[1:] * q[1:], u[0] * q[0])  # update zprava
		Fl = np.append(u[-1] * q[-1], u[:-1] * q[:-1])  # update zleva
	else:
		q_L = np.append(q[0], q[:-1])
		q_R = np.append(q[1:], q[-1])
		Fp = np.append(u[1:] * q[1:], q[-1])  # update zprava
		Fl = np.append(q[0], u[:-1] * q[:-1])  # update zleva

	res = q - dtdx * (Fp - Fl) / 2 + (u * dtdx) ** 2 * (q_L + q_R - 2 * q) / 2  # [(46) p. 77]

	return res, None


def leith_cons_1d(q, u, dtdx, bc, dx, vertices):
	"""
	Based on

	Leonard, B. P. (1979). A Stable And Accurate Convective Modelling Procedure
		Based on Quadraticu pstream Interpolation.
	Computer Methods in Applied Mechanics and Engineering, 19, 59–98.

	Just really basic streaming estimation procedure for comparison with QUISKEST

	Stability condition is given by
	0 <= alpha <= (1-c**2)/2

	:param q:
	:param u:
	:param dtdx:
	:return:
	"""

	if bc == 'zeroDirichlet':
		q_L = np.append(0,  q[:-1])
		q_R = np.append(q[1:], 0)
	elif bc == 'periodic':
		q_L = np.append(q[-1], q[:-1])
		q_R = np.append(q[1:], q[0])
	else:
		q_L = np.append(q[0], q[:-1])
		q_R = np.append(q[1:], q[-1])

	c_L = np.append(u[1], u[:-1]) * dtdx
	c_R = np.append(u[1:], u[-1]) * dtdx

	res = q + c_L / 2 * ((q_L + q) - c_L*(q - q_L)) - c_R / 2 * ((q + q_R) - c_R * (q_R - q))  # [(45) p. 76]

	return res, None


def quickest_noncons_1d(q, u, dtdx, bc, dx, vertices):
	"""

	QUICK with Estimated Streaming Terms

	Leonard, B. P. (1979). A Stable And Accurate Convective Modelling Procedure
		Based on Quadraticu pstream Interpolation.
	Computer Methods in Applied Mechanics and Engineering, 19, 59–98.

	Stability conditions can be found in Fig.18 p.82

	:param q:
	:param u:
	:param dtdx:
	:param dx:
	:return:
	"""

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

	grad_l = (q - q_L)/dx
	curv_l = (q_FL + q - 2*q_L)/dx**2  # [(52), p.79]

	grad_r = (q_R - q)/dx
	curv_r = (q_L + q_R - 2*q)/dx**2  # [(53), p.79]

	c = u * dtdx  # Courant number

	# diff is diffusion coefficient currently defined globally in the module
	res = q - c * ((1/2*(q+q_R) - dx/2*c*grad_r - dx**2/6*(1-c**2-3*diff)*curv_r)
				 - (1/2*(q_L+q) - dx/2*c*grad_l - dx**2/6*(1-c**2-3*diff)*curv_l))\
				 + diff*((dx*grad_r-dx**2*c*curv_r) - (dx*grad_l - dx**2/2*c*curv_l))  # [(57), p.80]

	return res, None


def quickest_cons_1d(q, u, dtdx, bc, dx, vertices):
	"""

	QUICK with Estimated Streaming Terms for conservative problem:
		q_t + (u(x)q)_x = 0 i.e. q_t + u(x)q_x + u(x)_xq = 0

	Leonard, B. P. (1979). A Stable And Accurate Convective Modelling Procedure
		Based on Quadraticu pstream Interpolation.
	Computer Methods in Applied Mechanics and Engineering, 19, 59–98.

	Stability conditions can be found in Fig.18 p.82

	:param q:
	:param u:
	:param dtdx:
	:param dx:
	:return:
	"""

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

	u_L = np.append(u[1], u[:-1])
	u_R = np.append(u[1:], u[-1])
	u_FL = np.append(u[:2], u[:-2])

	grad_l = (u*q - u_L*q_L)/dx
	curv_l = (u_FL*q_FL + u*q - 2*u_L*q_L)/dx**2  # [(52), p.79]

	grad_r = (u_R*q_R - u*q)/dx
	curv_r = (u_L*q_L + u_R*q_R - 2*u*q)/dx**2  # [(53), p.79]

	c = u * dtdx  # Courant number

	# diff is diffusion coefficient currently defined globally in the module
	# TODO is c in the following right? it stays the same as in nonconservative case, that might not be correct
	# as it probably does not align correct speed values to the interfaces
	res = q - dtdx * ((1/2*(u*q + u_R*q_R) - dx/2*c*grad_r - dx**2/6*(1-c**2-3*diff)*curv_r)
				    - (1/2*(u_L*q_L + u*q) - dx/2*c*grad_l - dx**2/6*(1-c**2-3*diff)*curv_l)) \
			+ diff*((dx*grad_r-dx**2*dtdx*curv_r) - (dx*grad_l - dx**2/2*dtdx*curv_l))  # [(57), p.80]

	return res, None


def quickest_cons_1dv2(q, u, dtdx, bc, dx, vertices):
	"""

	QUICK with Estimated Streaming Terms for conservative problem:
		q_t + (u(x)q)_x = 0 i.e. q_t + u(x)q_x + u(x)_xq = 0
	this is second more stremlined version without the diffusion
	coeficient and with nicer for of the update formula

	Leonard, B. P. (1979). A Stable And Accurate Convective Modelling Procedure
		Based on Quadraticu pstream Interpolation.
	Computer Methods in Applied Mechanics and Engineering, 19, 59–98.

	Stability conditions can be found in Fig.18 p.82

	:param q:
	:param u:
	:param dtdx:
	:param dx:
	:return:
	"""

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

	grad_l = (q - q_L)/dx
	curv_l = (q_FL + q - 2*q_L)/dx**2  # [(52), p.79]

	grad_r = (q_R - q)/dx
	curv_r = (q_L + q_R - 2*q)/dx**2  # [(53), p.79]

	c_L = np.append(u[1], u[:-1]) * dtdx
	c_R = np.append(u[1:], u[-1]) * dtdx

	res = q - (c_R * (1/2*(q + q_R) - dx/2*c_R*grad_r - dx**2/6*(1-c_R**2)*curv_r)
			-  c_L * (1/2*(q_L + q) - dx/2*c_L*grad_l - dx**2/6*(1-c_L**2)*curv_l))
	# derived from (51), (54) [p. 79] and (55), (56) [p. 80]

	return res, None


def afs_noncons_1d(q, u, dtdx, bc, dx, vertices):
	"""
	AFS

	Eymann, Timothy Andrew (2013). Active Flux Schemes. Dissertation in the University of Michigan (Aerospace Engineering and Scientific Computing)

	Scheme p. 18

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

	ni = u * dtdx
	if min(u) >= 0:
		res = ni ** 2 * (ni - 1) * q_LL_mid + ni ** 2 * (3 - 2 * ni) * q_L + ni * (1 - ni) * q_L_mid + (1 - ni) ** 2 * (
				1 + 2 * ni) * q - ni * (1 - ni) ** 2 * q_R_mid
		vertices = ni * (3 * ni - 2) * q_L_mid + 6 * ni * (1 - ni) * q + (1 - ni) * (1 - 3 * ni) * q_R_mid
	else:
		res = - ni ** 2 * (ni + 1) * q_RR_mid + ni ** 2 * (3 + 2 * ni) * q_R - ni * (1 + ni) * q_R_mid + (1 + ni) ** 2 * (
				1 - 2 * ni) * q + ni * (1 + ni) ** 2 * q_L_mid
		vertices = ni * (3 * ni + 2) * q_R_mid - 6 * ni * (1 + ni) * q + (1 + ni) * (1 + 3 * ni) * q_L_mid
	return res, vertices


def upwind_cons_1d(q, u, dtdx, bc, dx, vertices):
	"""
	First order iterative upwind
	method for equation
		q_t + (u(x)q)_x = 0  (9.8)
	(conservative form)

	:param q: solution in prevorius time
	:param u: sampled velocity
	:param dtdx: dt/dx
	:return: solution in current time step
	"""
	if min(u) >= 0:
		res = q[1:] - dtdx * (u[1:] * q[1:] - u[:-1] * q[:-1])
		if bc == 'zeroDirichlet':
			res = np.append(0, res)  # zero boundary condition on left side
		elif bc == 'periodic':
			res = np.append(q[-1], res)
	else:
		res = q[:-1] - dtdx * (u[1:] * q[1:] - u[:-1] * q[:-1])
		if bc == 'zeroDirichlet':
			res = np.append(res, 0)  # zero boundary condition on right side
		elif bc == 'periodic':
			res = np.append(res, q[0])
	return res, None


def upwind_noncons_1d(q, u, dtdx, bc, dx, vertices):
	"""
	First order iterative upwind
	method for equation
		q_t + u(x) * q_x = 0  (9.12)
	(nonconservative form)

	:param q: solution in prevorius time
	:param u: sampled velocity
	:param dtdx: dt/dx
	:return:solution in current time step
	"""
	dq = q[1:] - q[:-1]

	if min(u) >= 0:
		res = q[1:] - dtdx * u[1:] * dq
		if bc == 'zeroDirichlet':
			res = np.append(0, res)
		elif bc == 'periodic':
			res = np.append(q[-1], res)
	else:
		res = q[:-1] - dtdx * (u[:-1] * dq)
		if bc == 'zeroDirichlet':
			res = np.append(res, 0)
		elif bc == 'periodic':
			res = np.append(res, q[0])
	return res, None


def upwind_hrc_noncons_1d(q, u, dtdx, bc, dx, vertices):
	"""
	First order iterative upwind
	method with high resolution correction for equation
		q_t + u(x) * q_x = 0  (9.12)
	(nonconservative form)

	:param q: solution in prevorius time step
	:param u: sampled velocity
	:param dtdx: dt/dx
	:return: solution in current time step
	"""
	dp = q[1:] - q[:-1]

	q_R = (u[1:] * (1 - dtdx * u[1:])) * dp / 2
	Fm = (u[:-1] * (1 - dtdx * u[:-1])) * dp / 2
	if min(u) >= 0:
		res = q[1:] - dtdx * u[1:] * dp - dtdx * (-q_R + Fm)
		if bc == 'zeroDirichlet':
			res = np.append(0, res)
		elif bc == 'periodic':
			res = np.append(q[-1], res)
	else:
		res = q[:-1] - dtdx * (u[:-1] * dp) - dtdx * (q_R - Fm)
		if bc == 'zeroDirichlet':
			res = np.append(res, 0)
		elif bc == 'periodic':
			res = np.append(res, q[0])
	return res, None
