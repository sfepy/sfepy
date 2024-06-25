# -*- Mode: Python -*-
# cython: language_level=3
"""
Low level term evaluation functions.
"""
cimport cython

cimport numpy as np
import numpy as np

from sfepy.discrete.common.extmods.cmapping cimport (Mapping, CMapping)
from sfepy.discrete.common.extmods._fmfield cimport (FMField,
                                                     array2fmfield4,
                                                     array2fmfield3,
                                                     array2fmfield2,
                                                     array2fmfield1,
                                                     array2pint1,
                                                     array2pint2)

from sfepy.discrete.common.extmods.types cimport int32, float64, complex128

cdef extern from 'common.h':
    cdef void _errclear 'errclear'()

cdef extern from 'terms.h':
    cdef int32 _dq_state_in_qp \
         'dq_state_in_qp'(FMField *out, FMField *state, int32 offset,
                          FMField *bf,
                          int32 *conn, int32 nEl, int32 nEP)

    cdef int32 _dq_grad \
         'dq_grad'(FMField *out, FMField *state, int32 offset,
                   Mapping *vg, int32 *conn, int32 nEl, int32 nEP)

    cdef int32 _dq_div_vector \
         'dq_div_vector'(FMField *out, FMField *state, int32 offset,
                         Mapping *vg,
                         int32 *conn, int32 nEl, int32 nEP)

    cdef int32 _d_volume_surface \
         'd_volume_surface'(FMField *out, FMField *in_,
                            Mapping *sg,
                            int32 *conn, int32 nEl, int32 nEP)

    cdef int32 _di_surface_moment \
         'di_surface_moment'(FMField *out, FMField *in_,
                             Mapping *sg,
                             int32 *conn, int32 nEl, int32 nEP)

    cdef int32 _dq_finite_strain_tl \
         'dq_finite_strain_tl'(FMField *mtxF, FMField *detF, FMField *vecCS,
                               FMField *trC, FMField *in2C, FMField *vecInvCS,
                               FMField *vecES,
                               FMField *state, int32 offset,
                               Mapping *vg,
                               int32 *conn, int32 nEl, int32 nEP)

    cdef int32 _dq_finite_strain_ul \
         'dq_finite_strain_ul'(FMField *mtxF, FMField *detF, FMField *vecBS,
                               FMField *trB, FMField *in2B, FMField *vecES,
                               FMField *state, int32 offset,
                               Mapping *vg,
                               int32 *conn, int32 nEl, int32 nEP)

    cdef int32 _dq_tl_finite_strain_surface \
         'dq_tl_finite_strain_surface'(FMField *mtxF, FMField *detF,
                                       FMField *mtxFI,
                                       FMField *state, int32 offset,
                                       Mapping *sg,
                                       int32 *fis, int32 nFa, int32 nFP,
                                       int32 *conn, int32 nEl, int32 nE)

    cdef int32 _dq_tl_he_stress_bulk \
         'dq_tl_he_stress_bulk'(FMField *out,FMField *mat,
                                FMField *detF, FMField *vecInvCS)

    cdef int32 _dq_ul_he_stress_bulk \
         'dq_ul_he_stress_bulk'(FMField *out,FMField *mat,
                                FMField *detF)

    cdef int32 _dq_tl_he_stress_bulk_active \
         'dq_tl_he_stress_bulk_active'(FMField *out,FMField *mat,
                                       FMField *detF, FMField *vecInvCS)

    cdef int32 _dq_tl_he_stress_neohook \
         'dq_tl_he_stress_neohook'(FMField *out, FMField *mat,
                                   FMField *detF, FMField *trC,
                                   FMField *vecInvCS)

    cdef int32 _dq_ul_he_stress_neohook \
         'dq_ul_he_stress_neohook'(FMField *out, FMField *mat,
                                   FMField *detF, FMField *trB, FMField *vecBS)

    cdef int32 _dq_tl_he_stress_mooney_rivlin \
         'dq_tl_he_stress_mooney_rivlin'(FMField *out, FMField *mat,
                                         FMField *detF, FMField *trC,
                                         FMField *vecInvCS, FMField *vecCS,
                                         FMField *in2C)

    cdef int32 _dq_ul_he_stress_mooney_rivlin \
         'dq_ul_he_stress_mooney_rivlin'(FMField *out, FMField *mat,
                                         FMField *detF, FMField *trB,
                                         FMField *vecBS, FMField *in2B)

    cdef int32 _dq_tl_he_tan_mod_bulk \
         'dq_tl_he_tan_mod_bulk'(FMField *out, FMField *mat,
                                 FMField *detF, FMField *vecInvCS)

    cdef int32 _dq_ul_he_tan_mod_bulk \
         'dq_ul_he_tan_mod_bulk'(FMField *out, FMField *mat, FMField *detF)

    cdef int32 _dq_tl_he_tan_mod_bulk_active \
         'dq_tl_he_tan_mod_bulk_active'(FMField *out, FMField *mat,
                                        FMField *detF, FMField *vecInvCS)

    cdef int32 _dq_tl_he_tan_mod_neohook \
         'dq_tl_he_tan_mod_neohook'(FMField *out, FMField *mat,
                                    FMField *detF, FMField *trC,
                                    FMField *vecInvCS)

    cdef int32 _dq_ul_he_tan_mod_neohook \
         'dq_ul_he_tan_mod_neohook'(FMField *out, FMField *mat,
                                    FMField *detF, FMField *trB,
                                    FMField *vecBS)

    cdef int32 _dq_tl_he_tan_mod_mooney_rivlin \
         'dq_tl_he_tan_mod_mooney_rivlin'(FMField *out, FMField *mat,
                                          FMField *detF, FMField *trC,
                                          FMField *vecInvCS, FMField *vecCS,
                                          FMField *in2C)

    cdef int32 _dq_ul_he_tan_mod_mooney_rivlin \
         'dq_ul_he_tan_mod_mooney_rivlin'(FMField *out, FMField *mat,
                                          FMField *detF, FMField *trB,
                                          FMField *vecBS, FMField *in2B)

    cdef int32 _dw_he_rtm \
         'dw_he_rtm'(FMField *out,
                     FMField *stress, FMField *tan_mod,
                     FMField *mtxF, FMField *detF,
                     Mapping *vg,
                     int32 isDiff, int32 mode_ul)

    cdef int32 _de_he_rtm \
         'de_he_rtm'(FMField *out,
                     FMField *stress, FMField *detF,
                     Mapping *vg,
                     int32 *elList, int32 elList_nRow,
                     int32 mode_ul)

    cdef int32 _dq_tl_stress_bulk_pressure \
         'dq_tl_stress_bulk_pressure'(FMField *out, FMField *pressure_qp,
                                      FMField *detF, FMField *vecInvCS)
    cdef int32 _dq_ul_stress_bulk_pressure \
         'dq_ul_stress_bulk_pressure'(FMField *out, FMField *pressure_qp,
                                      FMField *det)
    cdef int32 _dq_tl_tan_mod_bulk_pressure_u \
         'dq_tl_tan_mod_bulk_pressure_u'(FMField *out, FMField *pressure_qp,
                                         FMField *detF, FMField *vecInvCS)
    cdef int32 _dq_ul_tan_mod_bulk_pressure_u \
         'dq_ul_tan_mod_bulk_pressure_u'(FMField *out, FMField *pressure_qp,
                                         FMField *detF)

    cdef int32 _dw_tl_volume \
         'dw_tl_volume'(FMField *out, FMField *mtxF,
                        FMField *vecInvCS, FMField *detF,
                        Mapping *vgs, Mapping *vgv,
                        int32 transpose, int32 mode)
    cdef int32 _dw_ul_volume \
         'dw_ul_volume'(FMField *out, FMField *detF,
                        Mapping *vgs, Mapping *vgv,
                        int32 transpose, int32 mode)

    cdef int32 _dw_tl_diffusion \
         'dw_tl_diffusion'(FMField *out, FMField *pressure_grad,
                           FMField *mtxD, FMField *ref_porosity,
                           FMField *mtxF, FMField *detF,
                           Mapping *vg, int32 mode)

    cdef int32 _d_tl_surface_flux \
         "d_tl_surface_flux"( FMField *out, FMField *pressure_grad,
                              FMField *mtxD, FMField *ref_porosity,
                              FMField *mtxFI, FMField *detF,
                              Mapping *sg, int32 mode )

    cdef int32 _dw_tl_surface_traction \
         'dw_tl_surface_traction'(FMField *out, FMField *traction,
                                  FMField *detF, FMField *mtxFI,
                                  FMField *bf, Mapping *sg,
                                  int32 *fis, int32 nFa, int32 nFP,
                                  int32 mode)

    cdef int32 _d_tl_volume_surface \
         'd_tl_volume_surface'(FMField *out, FMField *coors,
                               FMField *detF, FMField *mtxFI,
                               FMField *bf, Mapping *sg,
                               int32 *conn, int32 nFa, int32 nFP)

    cdef int32 _dq_def_grad \
         'dq_def_grad'(FMField *out, FMField *state, Mapping *vg,
                       int32 *conn, int32 nEl, int32 nEP, int32 mode)

    cdef int32 _he_residuum_from_mtx \
         'he_residuum_from_mtx'(FMField *out, FMField *mtxD,
                                FMField *state,
                                int32 *conn, int32 nEl, int32 nEP,
                                int32 *elList, int32 elList_nRo)
    cdef int32 _he_eval_from_mtx \
         'he_eval_from_mtx'(FMField *out, FMField *mtxD,
                            FMField *stateV, FMField *stateU,
                            int32 *conn, int32 nEl, int32 nEP,
                            int32 *elList, int32 elList_nRo)

    cdef int32 _dw_laplace \
         'dw_laplace'(FMField *out, FMField *grad,
                      FMField *coef, Mapping *vg,
                      int32 isDiff)

    cdef int32 _d_laplace \
         'd_laplace'(FMField *out, FMField *gradP1, FMField *gradP2,
                     FMField *coef, Mapping *vg)
    cdef int32 _dw_diffusion \
         'dw_diffusion'(FMField *out, FMField *grad,
                        FMField *mtxD, Mapping *vg,
                        int32 isDiff)
    cdef int32 _d_diffusion \
         'd_diffusion'(FMField *out, FMField *gradP1, FMField *gradP2,
                       FMField *mtxD, Mapping *vg)
    cdef int32 _dw_diffusion_r \
         'dw_diffusion_r'(FMField *out, FMField *mtxD, Mapping *vg)
    cdef int32 _d_surface_flux \
         'd_surface_flux'(FMField *out, FMField *grad,
                          FMField *mtxD, Mapping *sg, int32 mode)
    cdef int32 _dw_surface_flux \
         'dw_surface_flux'(FMField *out, FMField *grad,
                           FMField *mat, FMField *bf, Mapping *sg,
                           int32 *fis, int32 nFa, int32 nFP, int32 mode)
    cdef int32 _dw_convect_v_grad_s \
         'dw_convect_v_grad_s'(FMField *out, FMField *val_v, FMField *grad_s,
                               Mapping *vvg, Mapping *svg,
                               int32 isDiff)

    cdef int32 _dw_lin_elastic \
         'dw_lin_elastic'(FMField *out, float64 coef, FMField *strain,
                          FMField *mtxD, Mapping *vg,
                          int32 isDiff)
    cdef int32 _d_lin_elastic \
         'd_lin_elastic'(FMField *out, float64 coef, FMField *strainV,
                         FMField *strainU, FMField *mtxD, Mapping *vg)

    cdef int32 _d_sd_lin_elastic \
         'd_sd_lin_elastic'(FMField *out, float64 coef, FMField *gradV,
                            FMField *gradU, FMField *gradW, FMField *mtxD,
                            Mapping *vg)

    cdef int32 _dw_lin_prestress \
         'dw_lin_prestress'(FMField *out, FMField *stress, Mapping *vg)

    cdef int32 _dw_lin_strain_fib \
         'dw_lin_strain_fib'(FMField *out, FMField *mtxD, FMField *mat,
                             Mapping *vg)

    cdef int32 _de_cauchy_strain \
         'de_cauchy_strain'(FMField *out, FMField *strain,
                            Mapping *vg, int32 mode)
    cdef int32 _de_cauchy_stress \
         'de_cauchy_stress'(FMField *out, FMField *strain,
                            FMField *mtxD,  Mapping *vg,
                            int32 mode)
    cdef int32 _dq_cauchy_strain \
         'dq_cauchy_strain'(FMField *out, FMField *state, int32 offset,
                            Mapping *vg,
                            int32 *conn, int32 nEl, int32 nEP)

    cdef int32 _dw_nonsym_elastic \
         'dw_nonsym_elastic'(FMField *out, FMField *grad, FMField *mtxD,
                             Mapping *vg, int32 isDiff)

    cdef int32 _dw_surface_ltr \
         'dw_surface_ltr'(FMField *out, FMField *traction, Mapping *sg)

    cdef int32 _dw_volume_lvf \
         'dw_volume_lvf'(FMField *out, FMField *forceQP, Mapping *vg)

    cdef int32 _dw_surface_v_dot_n_s \
         'dw_surface_v_dot_n_s'(FMField *out,
                                FMField *coef, FMField *val_qp,
                                Mapping *rsg,
                                Mapping *csg,
                                int32 isDiff)

    cdef int32 _dw_surface_s_v_dot_n \
         'dw_surface_s_v_dot_n'(FMField *out,
                                FMField *coef, FMField *val_qp,
                                Mapping *rsg,
                                Mapping *csg,
                                int32 isDiff)

    cdef int32 _dw_volume_dot_vector \
         'dw_volume_dot_vector'(FMField *out, FMField *coef, FMField *val_qp,
                                Mapping *rvg, Mapping *cvg,
                                int32 isDiff)

    cdef int32 _dw_volume_dot_scalar \
         'dw_volume_dot_scalar'(FMField *out, FMField *coef, FMField *val_qp,
                                Mapping *rvg, Mapping *cvg,
                                int32 isDiff)

    cdef int32 _dw_v_dot_grad_s_vw \
         'dw_v_dot_grad_s_vw'(FMField *out, FMField *coef, FMField *grad,
                              Mapping *vvg, Mapping *svg,
                              int32 isDiff)

    cdef int32 _dw_v_dot_grad_s_sw \
         'dw_v_dot_grad_s_sw'(FMField *out, FMField *coef, FMField *val_qp,
                              Mapping *vvg, Mapping *svg,
                              int32 isDiff)

    cdef int32 _term_ns_asm_div_grad \
         'term_ns_asm_div_grad'(FMField *out, FMField *grad,
                                FMField *viscosity, Mapping *vgv,
                                Mapping *vgs, int32 isDiff)

    cdef int32 _term_ns_asm_convect \
         'term_ns_asm_convect'(FMField *out, FMField *grad, FMField *state,
                               Mapping *vg, int32 isDiff)

    cdef int32 _dw_lin_convect \
         'dw_lin_convect'(FMField *out, FMField *grad, FMField *stateB,
                          Mapping *vg, int32 isDiff)

    cdef int32 _dw_div \
         'dw_div'(FMField *out, FMField *coef, FMField *div,
                  Mapping *svg, Mapping *vvg,
                  int32 isDiff)

    cdef int32 _dw_grad \
         'dw_grad'(FMField *out, FMField *coef, FMField *state,
                   Mapping *svg, Mapping *vvg,
                   int32 isDiff)

    cdef int32 _dw_st_pspg_c \
         'dw_st_pspg_c'(FMField *out,
                        FMField *stateB, FMField *stateU,
                        FMField *coef,
                        Mapping *vg_p, Mapping *vg_u,
                        int32 *conn, int32 nEl, int32 nEP,
                        int32 isDiff)

    cdef int32 _dw_st_supg_p \
         'dw_st_supg_p'(FMField *out,
                        FMField *stateB, FMField *gradP,
                        FMField *coef,
                        Mapping *vg_u, Mapping *vg_p,
                        int32 isDiff)

    cdef int32 _dw_st_supg_c \
         'dw_st_supg_c'(FMField *out,
                        FMField *stateB, FMField *stateU,
                        FMField *coef, Mapping *vg,
                        int32 *conn, int32 nEl, int32 nEP,
                        int32 isDiff)

    cdef int32 _dw_st_grad_div \
         'dw_st_grad_div'(FMField *out, FMField *div,
                          FMField *coef, Mapping *vg,
                          int32 isDiff)

    cdef int32 _dw_biot_grad \
         'dw_biot_grad'(FMField *out, float64 coef, FMField *pressure_qp,
                        FMField *mtxD, Mapping *svg, Mapping *vvg,
                        int32 isDiff)

    cdef int32 _dw_biot_div \
         'dw_biot_div'(FMField *out, float64 coef, FMField *strain,
                       FMField *mtxD, Mapping *svg, Mapping *vvg,
                       int32 isDiff)

    cdef int32 _d_biot_div \
         'd_biot_div'(FMField *out, float64 coef, FMField *state,
                      FMField *strain, FMField *mtxD, Mapping *vg)

    cdef int32 _dw_piezo_coupling \
         'dw_piezo_coupling'(FMField *out, FMField *strain,
                             FMField *charge_grad,
                             FMField *mtxG, Mapping *vg,
                             int32 mode)

    cdef int32 _d_piezo_coupling \
         'd_piezo_coupling'(FMField *out, FMField *strain,
                            FMField *charge_grad,
                            FMField *mtxG, Mapping *vg)

    cdef int32 _dw_electric_source \
         'dw_electric_source'(FMField *out, FMField *grad, FMField *coef,
                              Mapping *vg)

    cdef int32 _d_sd_diffusion \
         'd_sd_diffusion'(FMField *out,
                          FMField *grad_q, FMField *grad_p,
                          FMField *grad_w, FMField *div_w,
                          FMField *mtxD, Mapping *vg)

    cdef int32 _dw_adj_convect1 \
         'dw_adj_convect1'(FMField *out, FMField *stateW, FMField *gradU,
                           Mapping *vg, int32 isDiff)

    cdef int32 _dw_adj_convect2 \
         'dw_adj_convect2'(FMField *out, FMField *stateW, FMField *stateU,
                           Mapping *vg, int32 isDiff)

    cdef int32 _dw_st_adj_supg_c \
         'dw_st_adj_supg_c'(FMField *out, FMField *stateW,
                            FMField *stateU, FMField *gradU,
                            FMField *coef, Mapping *vg,
                            int32 *conn, int32 nEl, int32 nEP,
                            int32 isDiff)

    cdef int32 _dw_st_adj1_supg_p \
         'dw_st_adj1_supg_p'(FMField *out, FMField *stateW, FMField *gradP,
                             FMField *coef, Mapping *vg_w,
                             int32 *conn_w, int32 nEl_w, int32 nEP_w,
                             int32 isDiff)

    cdef int32 _dw_st_adj2_supg_p \
         'dw_st_adj2_supg_p'(FMField *out, FMField *gradU, FMField *stateR,
                             FMField *coef,
                             Mapping *vg_u, Mapping *vg_r,
                             int32 *conn_r, int32 nEl_r, int32 nEP_r,
                             int32 isDiff)

    cdef int32 _d_of_nsMinGrad \
         'd_of_nsMinGrad'(FMField *out, FMField *grad,
                          FMField *viscosity, Mapping *vg)

    cdef int32 _d_of_nsSurfMinDPress \
         'd_of_nsSurfMinDPress'(FMField *out, FMField *pressure,
                                float64 weight, float64 bpress,
                                Mapping *sg, int32 isDiff)

    cdef int32 _d_sd_div \
         'd_sd_div'(FMField *out, FMField *divU, FMField *gradU,
                    FMField *stateP, FMField *divMV, FMField *gradMV,
                    Mapping *vg_u, int32 mode)

    cdef int32 _d_sd_div_grad \
         'd_sd_div_grad'(FMField *out, FMField *gradU, FMField *gradW,
                         FMField *divMV, FMField *gradMV, FMField *viscosity,
                         Mapping *vg_u, int32 mode)

    cdef int32 _d_sd_convect \
         'd_sd_convect'(FMField *out, FMField *stateU, FMField *gradU,
                        FMField *stateW, FMField *divMV, FMField *gradMV,
                        Mapping *vg_u, int32 mode)

    cdef int32 _d_sd_volume_dot \
         'd_sd_volume_dot'(FMField *out, FMField *stateP, FMField *stateQ,
                           FMField *divMV, Mapping *vg, int32 mode)

    cdef int32 _d_sd_st_grad_div \
         'd_sd_st_grad_div'(FMField *out, FMField *divU, FMField *gradU,
                            FMField *divW, FMField *gradW, FMField *divMV,
                            FMField *gradMV, FMField *coef,
                            Mapping *vg_u, int32 mode)

    cdef int32 _d_sd_st_supg_c \
         'd_sd_st_supg_c'(FMField *out, FMField *stateB, FMField *gradU,
                          FMField *gradW, FMField *divMV, FMField *gradMV,
                          FMField *coef, Mapping *vg_u, int32 mode)

    cdef int32 _d_sd_st_pspg_c \
         'd_sd_st_pspg_c'(FMField *out, FMField *stateB, FMField *gradU,
                          FMField *gradR, FMField *divMV, FMField *gradMV,
                          FMField *coef, Mapping *vg_u, int32 mode)

    cdef int32 _d_sd_st_pspg_p \
         'd_sd_st_pspg_p'(FMField *out, FMField *gradR, FMField *gradP,
                          FMField *divMV, FMField *gradMV, FMField *coef,
                          Mapping *vg_p, int32 mode)

    cdef int32 _mulAB_integrate \
         'mulAB_integrate'(FMField *out,
                           FMField *A, FMField *B,
                           Mapping *vg, int32 mode)

    cdef int32 _actBfT \
         'actBfT'(FMField *out,
                  FMField *bf, FMField *A)

    cdef int32 _sym2nonsym \
         'sym2nonsym'(FMField *out,
                      FMField *A)

def errclear():
    _errclear()

def dq_state_in_qp(np.ndarray out not None,
                   np.ndarray state not None,
                   np.ndarray bf not None,
                   np.ndarray conn not None):
    cdef int32 ret
    cdef FMField[1] _out, _state, _bf
    cdef (int32 *) _conn
    cdef int32 n_el, n_ep

    array2fmfield4(_out, out)
    array2fmfield1(_state, state)
    array2fmfield4(_bf, bf)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    ret = _dq_state_in_qp(_out, _state, 0, _bf, _conn, n_el, n_ep)
    return ret

def dq_grad(np.ndarray out not None,
            np.ndarray state not None,
            CMapping cmap not None,
            np.ndarray conn not None):
    cdef int32 ret
    cdef FMField[1] _out, _state
    cdef (int32 *) _conn
    cdef int32 n_el, n_ep

    array2fmfield4(_out, out)
    array2fmfield1(_state, state)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    ret = _dq_grad(_out, _state, 0, cmap.geo, _conn, n_el, n_ep)
    return ret

def dq_div_vector(np.ndarray out not None,
                  np.ndarray state not None,
                  CMapping cmap not None,
                  np.ndarray conn not None):
    cdef int32 ret
    cdef FMField[1] _out, _state
    cdef (int32 *) _conn
    cdef int32 n_el, n_ep

    array2fmfield4(_out, out)
    array2fmfield1(_state, state)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    ret = _dq_div_vector(_out, _state, 0, cmap.geo, _conn, n_el, n_ep)
    return ret

def d_volume_surface(np.ndarray out not None,
                     np.ndarray in_ not None,
                     CMapping cmap not None,
                     np.ndarray conn not None):
    cdef int32 ret
    cdef FMField[1] _out, _in_
    cdef (int32 *) _conn
    cdef int32 n_el, n_ep

    array2fmfield4(_out, out)
    array2fmfield2(_in_, in_)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    ret = _d_volume_surface(_out, _in_, cmap.geo, _conn, n_el, n_ep)
    return ret

def di_surface_moment(np.ndarray out not None,
                      np.ndarray in_ not None,
                      CMapping cmap not None,
                      np.ndarray conn not None):
    cdef int32 ret
    cdef FMField[1] _out, _in_
    cdef (int32 *) _conn
    cdef int32 n_el, n_ep

    array2fmfield4(_out, out)
    array2fmfield2(_in_, in_)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    ret = _di_surface_moment(_out, _in_, cmap.geo, _conn, n_el, n_ep)
    return ret

def dq_finite_strain_tl(np.ndarray mtx_f not None,
                        np.ndarray det_f not None,
                        np.ndarray vec_cs not None,
                        np.ndarray tr_c not None,
                        np.ndarray in_2c not None,
                        np.ndarray vec_inv_cs not None,
                        np.ndarray vec_es not None,
                        np.ndarray state not None,
                        CMapping cmap not None,
                        np.ndarray conn not None):
    cdef int32 ret
    cdef FMField[1] _mtx_f, _det_f, _vec_cs, _tr_c, _in_2c
    cdef FMField[1] _vec_inv_cs, _vec_es, _state
    cdef (int32 *) _conn
    cdef int32 n_el, n_ep

    array2fmfield4(_mtx_f, mtx_f)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_vec_cs, vec_cs)
    array2fmfield4(_tr_c, tr_c)
    array2fmfield4(_in_2c, in_2c)
    array2fmfield4(_vec_inv_cs, vec_inv_cs)
    array2fmfield4(_vec_es, vec_es)
    array2fmfield1(_state, state)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    ret = _dq_finite_strain_tl(_mtx_f, _det_f, _vec_cs, _tr_c, _in_2c,
                               _vec_inv_cs, _vec_es, _state, 0, cmap.geo,
                               _conn, n_el, n_ep)
    if ret:
        raise ValueError('ccore error (see above)')

    return ret

def dq_finite_strain_ul(np.ndarray mtx_f not None,
                        np.ndarray det_f not None,
                        np.ndarray vec_bs not None,
                        np.ndarray tr_b not None,
                        np.ndarray in_2b not None,
                        np.ndarray vec_es not None,
                        np.ndarray state not None,
                        CMapping cmap not None,
                        np.ndarray conn not None):
    cdef int32 ret
    cdef FMField[1] _mtx_f, _det_f, _vec_bs, _tr_b, _in_2b
    cdef FMField[1] _vec_es, _state
    cdef int32 *_conn
    cdef int32 n_el, n_ep

    array2fmfield4(_mtx_f, mtx_f)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_vec_bs, vec_bs)
    array2fmfield4(_tr_b, tr_b)
    array2fmfield4(_in_2b, in_2b)
    array2fmfield4(_vec_es, vec_es)
    array2fmfield1(_state, state)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    ret = _dq_finite_strain_ul(_mtx_f, _det_f, _vec_bs, _tr_b, _in_2b,
                               _vec_es, _state, 0, cmap.geo,
                               _conn, n_el, n_ep)
    if ret:
        raise ValueError('ccore error (see above)')

    return ret

def dq_tl_finite_strain_surface(np.ndarray mtx_f not None,
                                np.ndarray det_f not None,
                                np.ndarray mtx_fi not None,
                                np.ndarray state not None,
                                CMapping cmap not None,
                                np.ndarray fis not None,
                                np.ndarray conn not None):
    cdef int32 ret
    cdef FMField[1] _mtx_f, _det_f, _mtx_fi, _state
    cdef (int32 *) _conn, _fis
    cdef int32 n_el, n_ep, n_fa, n_fp

    array2fmfield4(_mtx_f, mtx_f)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_mtx_fi, mtx_fi)
    array2fmfield1(_state, state)
    array2pint2(&_fis, &n_fa, &n_fp, fis)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    ret = _dq_tl_finite_strain_surface(_mtx_f, _det_f, _mtx_fi, _state, 0,
                                       cmap.geo,
                                       _fis, n_fa, n_fp, _conn, n_el, n_ep)
    if ret:
        raise ValueError('ccore error (see above)')

    return ret

def dq_tl_he_stress_bulk(np.ndarray out not None,
                         np.ndarray mat not None,
                         np.ndarray det_f not None,
                         np.ndarray vec_inv_cs not None):
    cdef int32 ret
    cdef FMField[1] _out, _mat, _det_f, _vec_inv_cs

    array2fmfield4(_out, out)
    array2fmfield4(_mat, mat)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_vec_inv_cs, vec_inv_cs)

    ret = _dq_tl_he_stress_bulk(_out, _mat, _det_f, _vec_inv_cs)
    return ret

def dq_ul_he_stress_bulk(np.ndarray out not None,
                         np.ndarray mat not None,
                         np.ndarray det_f not None):
    cdef int32 ret
    cdef FMField[1] _out, _mat, _det_f

    array2fmfield4(_out, out)
    array2fmfield4(_mat, mat)
    array2fmfield4(_det_f, det_f)

    ret = _dq_ul_he_stress_bulk(_out, _mat, _det_f)
    return ret

def dq_tl_he_stress_bulk_active(np.ndarray out not None,
                                np.ndarray mat not None,
                                np.ndarray det_f not None,
                                np.ndarray vec_inv_cs not None):
    cdef int32 ret
    cdef FMField[1] _out, _mat, _det_f, _vec_inv_cs

    array2fmfield4(_out, out)
    array2fmfield4(_mat, mat)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_vec_inv_cs, vec_inv_cs)

    ret = _dq_tl_he_stress_bulk_active(_out, _mat, _det_f, _vec_inv_cs)
    return ret

def dq_tl_he_stress_neohook(np.ndarray out not None,
                            np.ndarray mat not None,
                            np.ndarray det_f not None,
                            np.ndarray tr_c not None,
                            np.ndarray vec_inv_cs not None):
    cdef int32 ret
    cdef FMField[1] _out, _mat, _det_f, _tr_c, _vec_inv_cs

    array2fmfield4(_out, out)
    array2fmfield4(_mat, mat)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_tr_c, tr_c)
    array2fmfield4(_vec_inv_cs, vec_inv_cs)

    ret = _dq_tl_he_stress_neohook(_out, _mat, _det_f, _tr_c, _vec_inv_cs)
    return ret

def dq_ul_he_stress_neohook(np.ndarray out not None,
                            np.ndarray mat not None,
                            np.ndarray det_f not None,
                            np.ndarray tr_b not None,
                            np.ndarray vec_bs not None):
    cdef int32 ret
    cdef FMField[1] _out, _mat, _det_f, _tr_b, _vec_bs

    array2fmfield4(_out, out)
    array2fmfield4(_mat, mat)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_tr_b, tr_b)
    array2fmfield4(_vec_bs, vec_bs)

    ret = _dq_ul_he_stress_neohook(_out, _mat, _det_f, _tr_b, _vec_bs)
    return ret

def dq_tl_he_stress_mooney_rivlin(np.ndarray out not None,
                                  np.ndarray mat not None,
                                  np.ndarray det_f not None,
                                  np.ndarray tr_c not None,
                                  np.ndarray vec_inv_cs not None,
                                  np.ndarray vec_cs not None,
                                  np.ndarray in_2c not None):
    cdef int32 ret
    cdef FMField[1] _out, _mat, _det_f, _tr_c, _vec_inv_cs
    cdef FMField[1] _vec_cs, _in_2c

    array2fmfield4(_out, out)
    array2fmfield4(_mat, mat)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_tr_c, tr_c)
    array2fmfield4(_vec_inv_cs, vec_inv_cs)
    array2fmfield4(_vec_cs, vec_cs)
    array2fmfield4(_in_2c, in_2c)

    ret = _dq_tl_he_stress_mooney_rivlin(_out, _mat, _det_f, _tr_c,
                                         _vec_inv_cs, _vec_cs, _in_2c)
    return ret

def dq_ul_he_stress_mooney_rivlin(np.ndarray out not None,
                                  np.ndarray mat not None,
                                  np.ndarray det_f not None,
                                  np.ndarray tr_b not None,
                                  np.ndarray vec_bs not None,
                                  np.ndarray in_2b not None):
    cdef int32 ret
    cdef FMField[1] _out, _mat, _det_f, _tr_b, _vec_bs, _in_2b

    array2fmfield4(_out, out)
    array2fmfield4(_mat, mat)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_tr_b, tr_b)
    array2fmfield4(_vec_bs, vec_bs)
    array2fmfield4(_in_2b, in_2b)

    ret = _dq_ul_he_stress_mooney_rivlin(_out, _mat, _det_f, _tr_b,
                                         _vec_bs, _in_2b)
    return ret

def dq_tl_he_tan_mod_bulk(np.ndarray out not None,
                          np.ndarray mat not None,
                          np.ndarray det_f not None,
                          np.ndarray vec_inv_cs not None):
    cdef int32 ret
    cdef FMField[1] _out, _mat, _det_f, _vec_inv_cs

    array2fmfield4(_out, out)
    array2fmfield4(_mat, mat)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_vec_inv_cs, vec_inv_cs)

    ret = _dq_tl_he_tan_mod_bulk(_out, _mat, _det_f, _vec_inv_cs)
    return ret

def dq_ul_he_tan_mod_bulk(np.ndarray out not None,
                          np.ndarray mat not None,
                          np.ndarray det_f not None):
    cdef int32 ret
    cdef FMField[1] _out, _mat, _det_f

    array2fmfield4(_out, out)
    array2fmfield4(_mat, mat)
    array2fmfield4(_det_f, det_f)

    ret = _dq_ul_he_tan_mod_bulk(_out, _mat, _det_f)
    return ret

def dq_tl_he_tan_mod_bulk_active(np.ndarray out not None,
                                 np.ndarray mat not None,
                                 np.ndarray det_f not None,
                                 np.ndarray vec_inv_cs not None):
    cdef int32 ret
    cdef FMField[1] _out, _mat, _det_f, _vec_inv_cs

    array2fmfield4(_out, out)
    array2fmfield4(_mat, mat)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_vec_inv_cs, vec_inv_cs)

    ret = _dq_tl_he_tan_mod_bulk_active(_out, _mat, _det_f, _vec_inv_cs)
    return ret

def dq_tl_he_tan_mod_neohook(np.ndarray out not None,
                             np.ndarray mat not None,
                             np.ndarray det_f not None,
                             np.ndarray tr_c not None,
                             np.ndarray vec_inv_cs not None):
    cdef int32 ret
    cdef FMField[1] _out, _mat, _det_f, _tr_c, _vec_inv_cs

    array2fmfield4(_out, out)
    array2fmfield4(_mat, mat)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_tr_c, tr_c)
    array2fmfield4(_vec_inv_cs, vec_inv_cs)

    ret = _dq_tl_he_tan_mod_neohook(_out, _mat, _det_f, _tr_c, _vec_inv_cs)
    return ret

def dq_ul_he_tan_mod_neohook(np.ndarray out not None,
                             np.ndarray mat not None,
                             np.ndarray det_f not None,
                             np.ndarray tr_b not None,
                             np.ndarray vec_bs not None):
    cdef int32 ret
    cdef FMField[1] _out, _mat, _det_f, _tr_b, _vec_bs

    array2fmfield4(_out, out)
    array2fmfield4(_mat, mat)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_tr_b, tr_b)
    array2fmfield4(_vec_bs, vec_bs)

    ret = _dq_ul_he_tan_mod_neohook(_out, _mat, _det_f, _tr_b, _vec_bs)
    return ret

def dq_tl_he_tan_mod_mooney_rivlin(np.ndarray out not None,
                                   np.ndarray mat not None,
                                   np.ndarray det_f not None,
                                   np.ndarray tr_c not None,
                                   np.ndarray vec_inv_cs not None,
                                   np.ndarray vec_cs not None,
                                   np.ndarray in_2c not None):
    cdef int32 ret
    cdef FMField[1] _out, _mat, _det_f, _tr_c, _vec_inv_cs
    cdef FMField[1] _vec_cs, _in_2c

    array2fmfield4(_out, out)
    array2fmfield4(_mat, mat)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_tr_c, tr_c)
    array2fmfield4(_vec_inv_cs, vec_inv_cs)
    array2fmfield4(_vec_cs, vec_cs)
    array2fmfield4(_in_2c, in_2c)

    ret = _dq_tl_he_tan_mod_mooney_rivlin(_out, _mat, _det_f, _tr_c,
                                          _vec_inv_cs, _vec_cs, _in_2c)
    return ret

def dq_ul_he_tan_mod_mooney_rivlin(np.ndarray out not None,
                                   np.ndarray mat not None,
                                   np.ndarray det_f not None,
                                   np.ndarray tr_b not None,
                                   np.ndarray vec_bs not None,
                                   np.ndarray in_2b not None):
    cdef int32 ret
    cdef FMField[1] _out, _mat, _det_f, _tr_b, _vec_bs, _in_2b

    array2fmfield4(_out, out)
    array2fmfield4(_mat, mat)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_tr_b, tr_b)
    array2fmfield4(_vec_bs, vec_bs)
    array2fmfield4(_in_2b, in_2b)

    ret = _dq_ul_he_tan_mod_mooney_rivlin(_out, _mat, _det_f, _tr_b,
                                          _vec_bs, _in_2b)
    return ret

def dw_he_rtm(np.ndarray out not None,
              np.ndarray stress not None,
              np.ndarray tan_mod not None,
              np.ndarray mtx_f not None,
              np.ndarray det_f not None,
              CMapping cmap not None,
              int32 is_diff, int32 mode_ul):
    cdef int32 ret
    cdef FMField[1] _out, _stress, _tan_mod, _mtx_f, _det_f

    array2fmfield4(_out, out)
    array2fmfield4(_stress, stress)
    array2fmfield4(_tan_mod, tan_mod)
    array2fmfield4(_mtx_f, mtx_f)
    array2fmfield4(_det_f, det_f)

    ret = _dw_he_rtm(_out, _stress, _tan_mod, _mtx_f, _det_f,
                     cmap.geo, is_diff, mode_ul)
    return ret

def de_he_rtm(np.ndarray out not None,
              np.ndarray stress not None,
              np.ndarray det_f not None,
              CMapping cmap not None,
              np.ndarray el_list not None,
              int32 mode_ul):
    cdef int32 ret
    cdef FMField[1] _out, _stress, _det_f
    cdef int32 *_el_list
    cdef int32 n_el

    array2fmfield4(_out, out)
    array2fmfield4(_stress, stress)
    array2fmfield4(_det_f, det_f)
    array2pint1(&_el_list, &n_el, el_list)

    ret = _de_he_rtm(_out, _stress, _det_f,
                     cmap.geo, _el_list, n_el, mode_ul)
    return ret

def dq_tl_stress_bulk_pressure(np.ndarray out not None,
                               np.ndarray pressure_qp not None,
                               np.ndarray det_f not None,
                               np.ndarray vec_inv_cs not None):
    cdef int32 ret
    cdef FMField[1] _out, _pressure_qp, _det_f, _vec_inv_cs

    array2fmfield4(_out, out)
    array2fmfield4(_pressure_qp, pressure_qp)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_vec_inv_cs, vec_inv_cs)

    ret = _dq_tl_stress_bulk_pressure(_out, _pressure_qp, _det_f, _vec_inv_cs)
    return ret

def dq_ul_stress_bulk_pressure(np.ndarray out not None,
                               np.ndarray pressure_qp not None,
                               np.ndarray det_f not None):
    cdef int32 ret
    cdef FMField[1] _out, _pressure_qp, _det_f

    array2fmfield4(_out, out)
    array2fmfield4(_pressure_qp, pressure_qp)
    array2fmfield4(_det_f, det_f)

    ret = _dq_ul_stress_bulk_pressure(_out, _pressure_qp, _det_f)
    return ret

def dq_tl_tan_mod_bulk_pressure_u(np.ndarray out not None,
                                  np.ndarray pressure_qp not None,
                                  np.ndarray det_f not None,
                                  np.ndarray vec_inv_cs not None):
    cdef int32 ret
    cdef FMField[1] _out, _pressure_qp, _det_f, _vec_inv_cs

    array2fmfield4(_out, out)
    array2fmfield4(_pressure_qp, pressure_qp)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_vec_inv_cs, vec_inv_cs)

    ret = _dq_tl_tan_mod_bulk_pressure_u(_out, _pressure_qp, _det_f,
                                         _vec_inv_cs)
    return ret

def dq_ul_tan_mod_bulk_pressure_u(np.ndarray out not None,
                                  np.ndarray pressure_qp not None,
                                  np.ndarray det_f not None):
    cdef int32 ret
    cdef FMField[1] _out, _pressure_qp, _det_f

    array2fmfield4(_out, out)
    array2fmfield4(_pressure_qp, pressure_qp)
    array2fmfield4(_det_f, det_f)

    ret = _dq_ul_tan_mod_bulk_pressure_u(_out, _pressure_qp, _det_f)
    return ret

def dw_tl_volume(np.ndarray out not None,
                 np.ndarray mtx_f not None,
                 np.ndarray vec_inv_cs not None,
                 np.ndarray det_f not None,
                 CMapping cmap_s not None,
                 CMapping cmap_v not None,
                 int32 transpose, int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _mtx_f, _vec_inv_cs, _det_f

    array2fmfield4(_out, out)
    array2fmfield4(_mtx_f, mtx_f)
    array2fmfield4(_vec_inv_cs, vec_inv_cs)
    array2fmfield4(_det_f, det_f)

    ret = _dw_tl_volume(_out, _mtx_f, _vec_inv_cs, _det_f,
                        cmap_s.geo, cmap_v.geo, transpose, mode)
    return ret

def dw_ul_volume(np.ndarray out not None,
                 np.ndarray det_f not None,
                 CMapping cmap_s not None,
                 CMapping cmap_v not None,
                 int32 transpose, int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _det_f

    array2fmfield4(_out, out)
    array2fmfield4(_det_f, det_f)

    ret = _dw_ul_volume(_out, _det_f, cmap_s.geo, cmap_v.geo, transpose, mode)
    return ret

def dw_tl_diffusion(np.ndarray out not None,
                    np.ndarray pressure_grad not None,
                    np.ndarray mtx_d not None,
                    np.ndarray ref_porosity not None,
                    np.ndarray mtx_f not None,
                    np.ndarray det_f not None,
                    CMapping cmap not None,
                    int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _pressure_grad, _mtx_d, _ref_porosity
    cdef FMField[1] _mtx_f, _det_f

    array2fmfield4(_out, out)
    array2fmfield4(_pressure_grad, pressure_grad)
    array2fmfield4(_mtx_d, mtx_d)
    array2fmfield4(_ref_porosity, ref_porosity)
    array2fmfield4(_mtx_f, mtx_f)
    array2fmfield4(_det_f, det_f)

    ret = _dw_tl_diffusion(_out, _pressure_grad, _mtx_d, _ref_porosity,
                           _mtx_f, _det_f, cmap.geo, mode)
    return ret

def d_tl_surface_flux(np.ndarray out not None,
                      np.ndarray pressure_grad not None,
                      np.ndarray mtx_d not None,
                      np.ndarray ref_porosity not None,
                      np.ndarray mtx_fi not None,
                      np.ndarray det_f not None,
                      CMapping cmap not None,
                      int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _pressure_grad, _mtx_d, _ref_porosity
    cdef FMField[1] _mtx_fi, _det_f

    array2fmfield4(_out, out)
    array2fmfield4(_pressure_grad, pressure_grad)
    array2fmfield4(_mtx_d, mtx_d)
    array2fmfield4(_ref_porosity, ref_porosity)
    array2fmfield4(_mtx_fi, mtx_fi)
    array2fmfield4(_det_f, det_f)

    ret = _d_tl_surface_flux(_out, _pressure_grad, _mtx_d, _ref_porosity,
                             _mtx_fi, _det_f, cmap.geo, mode)
    return ret

def dw_tl_surface_traction(np.ndarray out not None,
                           np.ndarray traction not None,
                           np.ndarray det_f not None,
                           np.ndarray mtx_fi not None,
                           np.ndarray bf not None,
                           CMapping cmap not None,
                           np.ndarray fis not None,
                           int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _traction, _det_f, _mtx_fi, _bf
    cdef int32 *_fis
    cdef int32 n_fa, n_fp

    array2fmfield4(_out, out)
    array2fmfield4(_traction, traction)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_mtx_fi, mtx_fi)
    array2fmfield4(_bf, bf)
    array2pint2(&_fis, &n_fa, &n_fp, fis)

    ret = _dw_tl_surface_traction(_out, _traction, _det_f, _mtx_fi, _bf,
                                       cmap.geo, _fis, n_fa, n_fp, mode)
    return ret

def d_tl_volume_surface(np.ndarray out not None,
                        np.ndarray coors not None,
                        np.ndarray det_f not None,
                        np.ndarray mtx_fi not None,
                        np.ndarray bf not None,
                        CMapping cmap not None,
                        np.ndarray conn not None):
    cdef int32 ret
    cdef FMField[1] _out, _coors, _det_f, _mtx_fi, _bf
    cdef int32 *_conn
    cdef int32 n_fa, n_fp

    array2fmfield4(_out, out)
    array2fmfield2(_coors, coors)
    array2fmfield4(_det_f, det_f)
    array2fmfield4(_mtx_fi, mtx_fi)
    array2fmfield4(_bf, bf)
    array2pint2(&_conn, &n_fa, &n_fp, conn)

    ret = _d_tl_volume_surface(_out, _coors, _det_f, _mtx_fi, _bf,
                               cmap.geo, _conn, n_fa, n_fp)
    return ret

def dq_def_grad(np.ndarray out not None,
                np.ndarray state not None,
                CMapping cmap not None,
                np.ndarray conn not None,
                int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _state
    cdef int32 *_conn
    cdef int32 n_el, n_ep

    array2fmfield4(_out, out)
    array2fmfield1(_state, state)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    ret = _dq_def_grad(_out, _state, cmap.geo, _conn, n_el, n_ep, mode)
    return ret

def he_residuum_from_mtx(np.ndarray out not None,
                         np.ndarray mtx_d not None,
                         np.ndarray state not None,
                         np.ndarray conn not None,
                         np.ndarray el_list not None):
    cdef int32 ret
    cdef FMField[1] _out, _mtx_d, _state
    cdef (int32 *) _conn, _el_list
    cdef int32 n_el, n_ep, n_el2

    array2fmfield4(_out, out)
    array2fmfield4(_mtx_d, mtx_d)
    array2fmfield1(_state, state)
    array2pint2(&_conn, &n_el, &n_ep, conn)
    array2pint1(&_el_list, &n_el2, el_list)

    ret = _he_residuum_from_mtx(_out, _mtx_d, _state,
                                _conn, n_el, n_ep, _el_list, n_el2)
    return ret

def he_eval_from_mtx(np.ndarray out not None,
                     np.ndarray mtx_d not None,
                     np.ndarray state_v not None,
                     np.ndarray state_u not None,
                     np.ndarray conn not None,
                     np.ndarray el_list not None):
    cdef int32 ret
    cdef FMField[1] _out, _mtx_d, _state_v, _state_u
    cdef (int32 *) _conn, _el_list
    cdef int32 n_el, n_ep, n_el2

    array2fmfield4(_out, out)
    array2fmfield4(_mtx_d, mtx_d)
    array2fmfield1(_state_v, state_v)
    array2fmfield1(_state_u, state_u)
    array2pint2(&_conn, &n_el, &n_ep, conn)
    array2pint1(&_el_list, &n_el2, el_list)

    ret = _he_eval_from_mtx(_out, _mtx_d, _state_v, _state_u,
                            _conn, n_el, n_ep, _el_list, n_el2)
    return ret

def dw_laplace(np.ndarray out not None,
               np.ndarray grad not None,
               np.ndarray coef not None,
               CMapping cmap not None,
               int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _grad, _coef

    array2fmfield4(_out, out)
    array2fmfield4(_grad, grad)
    array2fmfield4(_coef, coef)

    ret = _dw_laplace(_out, _grad, _coef, cmap.geo, is_diff)
    return ret

def d_laplace(np.ndarray out not None,
              np.ndarray grad_p1 not None,
              np.ndarray grad_p2 not None,
              np.ndarray coef not None,
              CMapping cmap not None):
    cdef int32 ret
    cdef FMField[1] _out, _grad_p1, _grad_p2, _coef

    array2fmfield4(_out, out)
    array2fmfield4(_grad_p1, grad_p1)
    array2fmfield4(_grad_p2, grad_p2)
    array2fmfield4(_coef, coef)

    ret = _d_laplace(_out, _grad_p1, _grad_p2, _coef, cmap.geo)
    return ret

def dw_diffusion(np.ndarray out not None,
                 np.ndarray grad not None,
                 np.ndarray mtx_d not None,
                 CMapping cmap not None,
                 int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _grad, _mtx_d

    array2fmfield4(_out, out)
    array2fmfield4(_grad, grad)
    array2fmfield4(_mtx_d, mtx_d)

    ret = _dw_diffusion(_out, _grad, _mtx_d, cmap.geo, is_diff)
    return ret

def d_diffusion(np.ndarray out not None,
                np.ndarray grad_p1 not None,
                np.ndarray grad_p2 not None,
                np.ndarray mtx_d not None,
                CMapping cmap not None):
    cdef int32 ret
    cdef FMField[1] _out, _grad_p1, _grad_p2, _mtx_d

    array2fmfield4(_out, out)
    array2fmfield4(_grad_p1, grad_p1)
    array2fmfield4(_grad_p2, grad_p2)
    array2fmfield4(_mtx_d, mtx_d)

    ret = _d_diffusion(_out, _grad_p1, _grad_p2, _mtx_d, cmap.geo)
    return ret

def dw_diffusion_r(np.ndarray out not None,
                   np.ndarray mtx_d not None,
                   CMapping cmap not None):
    cdef int32 ret
    cdef FMField[1] _out, _mtx_d

    array2fmfield4(_out, out)
    array2fmfield4(_mtx_d, mtx_d)

    ret = _dw_diffusion_r(_out, _mtx_d, cmap.geo)
    return ret

def d_surface_flux(np.ndarray out not None,
                   np.ndarray grad not None,
                   np.ndarray mtx_d not None,
                   CMapping cmap not None,
                   int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _grad, _mtx_d

    array2fmfield4(_out, out)
    array2fmfield4(_grad, grad)
    array2fmfield4(_mtx_d, mtx_d)

    ret = _d_surface_flux(_out, _grad, _mtx_d, cmap.geo, mode)
    return ret

def dw_surface_flux(np.ndarray out not None,
                    np.ndarray grad not None,
                    np.ndarray mat not None,
                    np.ndarray bf not None,
                    CMapping cmap not None,
                    np.ndarray fis not None,
                    int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _grad, _mat, _bf
    cdef int32 *_fis
    cdef int32 n_fa, n_fp

    array2fmfield4(_out, out)
    array2fmfield4(_grad, grad)
    array2fmfield4(_mat, mat)
    array2fmfield4(_bf, bf)
    array2pint2(&_fis, &n_fa, &n_fp, fis)

    ret = _dw_surface_flux(_out, _grad, _mat, _bf,
                           cmap.geo, _fis, n_fa, n_fp, mode)
    return ret

def dw_convect_v_grad_s(np.ndarray out not None,
                        np.ndarray val_v not None,
                        np.ndarray grad_s not None,
                        CMapping cmap_v not None,
                        CMapping cmap_s not None,
                        int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _val_v, _grad_s

    array2fmfield4(_out, out)
    array2fmfield4(_val_v, val_v)
    array2fmfield4(_grad_s, grad_s)

    ret = _dw_convect_v_grad_s(_out, _val_v, _grad_s,
                               cmap_v.geo, cmap_s.geo, is_diff)
    return ret

def dw_lin_elastic(np.ndarray out not None,
                   float64 coef,
                   np.ndarray strain not None,
                   np.ndarray mtx_d not None,
                   CMapping cmap not None,
                   int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _strain, _mtx_d

    array2fmfield4(_out, out)
    array2fmfield4(_strain, strain)
    array2fmfield4(_mtx_d, mtx_d)

    ret = _dw_lin_elastic(_out, coef, _strain, _mtx_d, cmap.geo, is_diff)
    return ret

def d_lin_elastic(np.ndarray out not None,
                  float64 coef,
                  np.ndarray strain_v not None,
                  np.ndarray strain_u not None,
                  np.ndarray mtx_d not None,
                  CMapping cmap not None):
    cdef int32 ret
    cdef FMField[1] _out, _strain_u, _strain_v, _mtx_d

    array2fmfield4(_out, out)
    array2fmfield4(_strain_u, strain_u)
    array2fmfield4(_strain_v, strain_v)
    array2fmfield4(_mtx_d, mtx_d)

    ret = _d_lin_elastic(_out, coef, _strain_v, _strain_u, _mtx_d, cmap.geo)
    return ret

def d_sd_lin_elastic(np.ndarray out not None,
                     float64 coef,
                     np.ndarray grad_v not None,
                     np.ndarray grad_u not None,
                     np.ndarray grad_w not None,
                     np.ndarray mtx_d not None,
                     CMapping cmap not None):
    cdef int32 ret
    cdef FMField[1] _out, _grad_u, _grad_v, _grad_w, _mtx_d

    array2fmfield4(_out, out)
    array2fmfield4(_grad_u, grad_u)
    array2fmfield4(_grad_v, grad_v)
    array2fmfield4(_grad_w, grad_w)
    array2fmfield4(_mtx_d, mtx_d)

    ret = _d_sd_lin_elastic(_out, coef, _grad_v, _grad_u, _grad_w,
                            _mtx_d, cmap.geo)
    return ret

def dw_lin_prestress(np.ndarray out not None,
                     np.ndarray stress not None,
                     CMapping cmap not None):
    cdef int32 ret
    cdef FMField[1] _out, _stress

    array2fmfield4(_out, out)
    array2fmfield4(_stress, stress)

    ret = _dw_lin_prestress(_out, _stress, cmap.geo)
    return ret

def dw_lin_strain_fib(np.ndarray out not None,
                      np.ndarray mtx_d not None,
                      np.ndarray mat not None,
                      CMapping cmap not None):
    cdef int32 ret
    cdef FMField[1] _out, _mtx_d, _mat

    array2fmfield4(_out, out)
    array2fmfield4(_mtx_d, mtx_d)
    array2fmfield4(_mat, mat)

    ret = _dw_lin_strain_fib(_out, _mtx_d, _mat, cmap.geo)
    return ret

def de_cauchy_strain(np.ndarray out not None,
                     np.ndarray strain not None,
                     CMapping cmap not None,
                     int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _strain

    array2fmfield4(_out, out)
    array2fmfield4(_strain, strain)

    ret = _de_cauchy_strain(_out, _strain, cmap.geo, mode)
    return ret

def de_cauchy_stress(np.ndarray out not None,
                     np.ndarray strain not None,
                     np.ndarray mtx_d not None,
                     CMapping cmap not None,
                     int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _strain, _mtx_d

    array2fmfield4(_out, out)
    array2fmfield4(_strain, strain)
    array2fmfield4(_mtx_d, mtx_d)

    ret = _de_cauchy_stress(_out, _strain, _mtx_d, cmap.geo, mode)
    return ret

def dq_cauchy_strain(np.ndarray out not None,
                     np.ndarray state not None,
                     CMapping cmap not None,
                     np.ndarray conn not None):
    cdef int32 ret
    cdef FMField[1] _out, _state
    cdef int32 *_conn
    cdef int32 n_el, n_ep

    array2fmfield4(_out, out)
    array2fmfield1(_state, state)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    ret = _dq_cauchy_strain(_out, _state, 0, cmap.geo, _conn, n_el, n_ep)
    return ret

def dw_nonsym_elastic(np.ndarray out not None,
                      np.ndarray grad not None,
                      np.ndarray mtx_d not None,
                      CMapping cmap not None,
                      int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _grad, _mtx_d

    array2fmfield4(_out, out)
    array2fmfield4(_grad, grad)
    array2fmfield4(_mtx_d, mtx_d)

    ret = _dw_nonsym_elastic(_out, _grad, _mtx_d, cmap.geo, is_diff)
    return ret

def dw_surface_ltr(np.ndarray out not None,
                   np.ndarray traction not None,
                   CMapping cmap not None):
    cdef int32 ret
    cdef FMField[1] _out, _traction

    array2fmfield4(_out, out)
    array2fmfield4(_traction, traction)

    ret = _dw_surface_ltr(_out, _traction, cmap.geo)
    return ret

def dw_volume_lvf(np.ndarray out not None,
                  np.ndarray force_qp not None,
                  CMapping cmap not None):
    cdef int32 ret
    cdef FMField[1] _out, _force_qp

    array2fmfield4(_out, out)
    array2fmfield4(_force_qp, force_qp)

    ret = _dw_volume_lvf(_out, _force_qp, cmap.geo)
    return ret

def dw_surface_v_dot_n_s(np.ndarray out not None,
                         np.ndarray coef not None,
                         np.ndarray val_qp not None,
                         CMapping rcmap not None,
                         CMapping ccmap not None,
                         int32 is_diff):

    cdef int32 ret
    cdef FMField[1] _out, _coef, _val_qp

    array2fmfield4(_out, out)
    array2fmfield4(_coef, coef)
    array2fmfield4(_val_qp, val_qp)

    ret = _dw_surface_v_dot_n_s(_out, _coef, _val_qp,
                                rcmap.geo, ccmap.geo, is_diff)
    return ret

def dw_surface_s_v_dot_n(np.ndarray out not None,
                         np.ndarray coef not None,
                         np.ndarray val_qp not None,
                         CMapping rcmap not None,
                         CMapping ccmap not None,
                         int32 is_diff):

    cdef int32 ret
    cdef FMField[1] _out, _coef, _val_qp

    array2fmfield4(_out, out)
    array2fmfield4(_coef, coef)
    array2fmfield4(_val_qp, val_qp)

    ret = _dw_surface_s_v_dot_n(_out, _coef, _val_qp,
                                rcmap.geo, ccmap.geo, is_diff)
    return ret

def dw_volume_dot_vector(np.ndarray out not None,
                         np.ndarray coef not None,
                         np.ndarray val_qp not None,
                         CMapping rcmap not None,
                         CMapping ccmap not None,
                         int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _coef, _val_qp

    array2fmfield4(_out, out)
    array2fmfield4(_coef, coef)
    array2fmfield4(_val_qp, val_qp)

    ret = _dw_volume_dot_vector(_out, _coef, _val_qp,
                                rcmap.geo, ccmap.geo, is_diff)
    return ret

def dw_volume_dot_scalar(np.ndarray out not None,
                         np.ndarray coef not None,
                         np.ndarray val_qp not None,
                         CMapping rcmap not None,
                         CMapping ccmap not None,
                         int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _coef, _val_qp

    array2fmfield4(_out, out)
    array2fmfield4(_coef, coef)
    array2fmfield4(_val_qp, val_qp)

    ret = _dw_volume_dot_scalar(_out, _coef, _val_qp,
                                rcmap.geo, ccmap.geo, is_diff)
    return ret

def dw_v_dot_grad_s_vw(np.ndarray out not None,
                       np.ndarray coef not None,
                       np.ndarray grad not None,
                       CMapping cmap_v not None,
                       CMapping cmap_s not None,
                       int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _coef, _grad

    array2fmfield4(_out, out)
    array2fmfield4(_coef, coef)
    array2fmfield4(_grad, grad)

    ret = _dw_v_dot_grad_s_vw(_out, _coef, _grad,
                              cmap_v.geo, cmap_s.geo, is_diff)
    return ret

def dw_v_dot_grad_s_sw(np.ndarray out not None,
                       np.ndarray coef not None,
                       np.ndarray val_qp not None,
                       CMapping cmap_v not None,
                       CMapping cmap_s not None,
                       int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _coef, _val_qp

    array2fmfield4(_out, out)
    array2fmfield4(_coef, coef)
    array2fmfield4(_val_qp, val_qp)

    ret = _dw_v_dot_grad_s_sw(_out, _coef, _val_qp,
                              cmap_v.geo, cmap_s.geo, is_diff)
    return ret

def term_ns_asm_div_grad(np.ndarray out not None,
                         np.ndarray grad not None,
                         np.ndarray viscosity not None,
                         CMapping cmap_v not None,
                         CMapping cmap_s not None,
                         int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _grad, _viscosity

    array2fmfield4(_out, out)
    array2fmfield4(_grad, grad)
    array2fmfield4(_viscosity, viscosity)

    ret = _term_ns_asm_div_grad(_out, _grad, _viscosity,
                                cmap_v.geo, cmap_s.geo, is_diff)
    return ret

def term_ns_asm_convect(np.ndarray out not None,
                        np.ndarray grad not None,
                        np.ndarray state not None,
                        CMapping cmap not None,
                        int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _grad, _state

    array2fmfield4(_out, out)
    array2fmfield4(_grad, grad)
    array2fmfield4(_state, state)

    ret = _term_ns_asm_convect(_out, _grad, _state, cmap.geo, is_diff)
    return ret

def dw_lin_convect(np.ndarray out not None,
                   np.ndarray grad not None,
                   np.ndarray state_b not None,
                   CMapping cmap not None,
                   int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _grad, _state_b

    array2fmfield4(_out, out)
    array2fmfield4(_grad, grad)
    array2fmfield4(_state_b, state_b)

    ret = _dw_lin_convect(_out, _grad, _state_b, cmap.geo, is_diff)
    return ret

def dw_div(np.ndarray out not None,
           np.ndarray coef not None,
           np.ndarray div not None,
           CMapping cmap_s not None,
           CMapping cmap_v not None,
           int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _coef, _div

    array2fmfield4(_out, out)
    array2fmfield4(_coef, coef)
    array2fmfield4(_div, div)

    ret = _dw_div(_out, _coef, _div, cmap_s.geo, cmap_v.geo, is_diff)
    return ret

def dw_grad(np.ndarray out not None,
            np.ndarray coef not None,
            np.ndarray state not None,
            CMapping cmap_s not None,
            CMapping cmap_v not None,
            int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _coef, _state

    array2fmfield4(_out, out)
    array2fmfield4(_coef, coef)
    array2fmfield4(_state, state)

    ret = _dw_grad(_out, _coef, _state, cmap_s.geo, cmap_v.geo, is_diff)
    return ret

def dw_st_pspg_c(np.ndarray out not None,
                 np.ndarray state_b not None,
                 np.ndarray state_u not None,
                 np.ndarray coef not None,
                 CMapping cmap_p not None,
                 CMapping cmap_u not None,
                 np.ndarray conn not None,
                 int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _state_b, _state_u, _coef
    cdef int32 *_conn
    cdef int32 n_el, n_ep

    array2fmfield4(_out, out)
    array2fmfield4(_state_b, state_b)
    array2fmfield1(_state_u, state_u)
    array2fmfield4(_coef, coef)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    ret = _dw_st_pspg_c(_out, _state_b, _state_u, _coef,
                        cmap_p.geo, cmap_u.geo, _conn, n_el, n_ep, is_diff)
    return ret

def dw_st_supg_p(np.ndarray out not None,
                 np.ndarray state_b not None,
                 np.ndarray grad_p not None,
                 np.ndarray coef not None,
                 CMapping cmap_u not None,
                 CMapping cmap_p not None,
                 int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _state_b, _grad_p, _coef

    array2fmfield4(_out, out)
    array2fmfield4(_state_b, state_b)
    array2fmfield4(_grad_p, grad_p)
    array2fmfield4(_coef, coef)

    ret = _dw_st_supg_p(_out, _state_b, _grad_p, _coef,
                        cmap_u.geo, cmap_p.geo, is_diff)
    return ret

def dw_st_supg_c(np.ndarray out not None,
                 np.ndarray state_b not None,
                 np.ndarray state_u not None,
                 np.ndarray coef not None,
                 CMapping cmap not None,
                 np.ndarray conn not None,
                 int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _state_b, _state_u, _coef
    cdef int32 *_conn
    cdef int32 n_el, n_ep

    array2fmfield4(_out, out)
    array2fmfield4(_state_b, state_b)
    array2fmfield1(_state_u, state_u)
    array2fmfield4(_coef, coef)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    ret = _dw_st_supg_c(_out, _state_b, _state_u, _coef,
                        cmap.geo, _conn, n_el, n_ep, is_diff)
    return ret

def dw_st_grad_div(np.ndarray out not None,
                   np.ndarray div not None,
                   np.ndarray coef not None,
                   CMapping cmap not None,
                   int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _div, _coef

    array2fmfield4(_out, out)
    array2fmfield4(_div, div)
    array2fmfield4(_coef, coef)

    ret = _dw_st_grad_div(_out, _div, _coef, cmap.geo, is_diff)
    return ret

def dw_biot_grad(np.ndarray out not None,
                 float64 coef,
                 np.ndarray pressure_qp not None,
                 np.ndarray mtx_d not None,
                 CMapping cmap_s not None,
                 CMapping cmap_v not None,
                 int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _pressure_qp, _mtx_d

    array2fmfield4(_out, out)
    array2fmfield4(_pressure_qp, pressure_qp)
    array2fmfield4(_mtx_d, mtx_d)

    ret = _dw_biot_grad(_out, coef, _pressure_qp, _mtx_d,
                        cmap_s.geo, cmap_v.geo, is_diff)
    return ret

def dw_biot_div(np.ndarray out not None,
                float64 coef,
                np.ndarray strain not None,
                np.ndarray mtx_d not None,
                CMapping cmap_s not None,
                CMapping cmap_v not None,
                int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _strain, _mtx_d

    array2fmfield4(_out, out)
    array2fmfield4(_strain, strain)
    array2fmfield4(_mtx_d, mtx_d)

    ret = _dw_biot_div(_out, coef, _strain, _mtx_d,
                       cmap_s.geo, cmap_v.geo, is_diff)
    return ret

def d_biot_div(np.ndarray out not None,
               float64 coef,
               np.ndarray state not None,
               np.ndarray strain not None,
               np.ndarray mtx_d not None,
               CMapping cmap not None):
    cdef int32 ret
    cdef FMField[1] _out, _state, _strain, _mtx_d

    array2fmfield4(_out, out)
    array2fmfield4(_state, state)
    array2fmfield4(_strain, strain)
    array2fmfield4(_mtx_d, mtx_d)

    ret = _d_biot_div(_out, coef, _state, _strain, _mtx_d, cmap.geo)
    return ret

def dw_piezo_coupling(np.ndarray out not None,
                      np.ndarray strain not None,
                      np.ndarray charge_grad not None,
                      np.ndarray mtx_g not None,
                      CMapping cmap not None,
                      int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _strain, _charge_grad, _mtx_g

    array2fmfield4(_out, out)
    array2fmfield4(_strain, strain)
    array2fmfield4(_charge_grad, charge_grad)
    array2fmfield4(_mtx_g, mtx_g)

    ret = _dw_piezo_coupling(_out, _strain, _charge_grad, _mtx_g,
                             cmap.geo, mode)
    return ret

def d_piezo_coupling(np.ndarray out not None,
                     np.ndarray strain not None,
                     np.ndarray charge_grad not None,
                     np.ndarray mtx_g not None,
                     CMapping cmap not None):
    cdef int32 ret
    cdef FMField[1] _out, _strain, _charge_grad, _mtx_g

    array2fmfield4(_out, out)
    array2fmfield4(_strain, strain)
    array2fmfield4(_charge_grad, charge_grad)
    array2fmfield4(_mtx_g, mtx_g)

    ret = _d_piezo_coupling(_out, _strain, _charge_grad, _mtx_g, cmap.geo)
    return ret

def dw_electric_source(np.ndarray out not None,
                       np.ndarray grad not None,
                       np.ndarray coef not None,
                       CMapping cmap not None):
    cdef int32 ret
    cdef FMField[1] _out, _grad, _coef

    array2fmfield4(_out, out)
    array2fmfield4(_grad, grad)
    array2fmfield4(_coef, coef)

    ret = _dw_electric_source(_out, _grad, _coef, cmap.geo)
    return ret

def d_sd_diffusion(np.ndarray out not None,
                   np.ndarray grad_q not None,
                   np.ndarray grad_p not None,
                   np.ndarray grad_w not None,
                   np.ndarray div_w not None,
                   np.ndarray mtx_d not None,
                   CMapping cmap not None):
    cdef int32 ret
    cdef FMField[1] _out, _grad_q, _grad_p, _grad_w, _div_w, _mtx_d

    array2fmfield4(_out, out)
    array2fmfield4(_grad_q, grad_q)
    array2fmfield4(_grad_p, grad_p)
    array2fmfield4(_grad_w, grad_w)
    array2fmfield4(_div_w, div_w)
    array2fmfield4(_mtx_d, mtx_d)

    ret = _d_sd_diffusion(_out, _grad_q, _grad_p, _grad_w, _div_w,
                          _mtx_d, cmap.geo)
    return ret

def mulAB_integrate(np.ndarray out not None,
                    np.ndarray A not None,
                    np.ndarray B not None,
                    CMapping cmap not None,
                    mode):
    cdef int32 ret
    cdef FMField[1] _out, _A, _B

    array2fmfield4(_out, out)
    if A.ndim == 4:
        array2fmfield4(_A, A)
    else:
        array2fmfield3(_A, A)

    if B.ndim == 4:
        array2fmfield4(_B, B)
    else:
        array2fmfield3(_B, B)

    if mode == 'ATB':
        imode = 0
    elif mode == 'AB':
        imode = 1
    elif mode == 'ABT':
        imode = 2
    elif mode == 'ATBT':
        imode = 3
    else:
        imode = -1

    ret = _mulAB_integrate(_out, _A, _B, cmap.geo, imode)
    return ret

def actBfT(np.ndarray out not None,
           np.ndarray bf not None,
           np.ndarray A not None):
    cdef int32 ret
    cdef FMField[1] _out, _bf, _A

    array2fmfield4(_out, out)
    array2fmfield4(_bf, bf)
    array2fmfield4(_A, A)

    ret = _actBfT(_out, _bf, _A)
    return ret

def sym2nonsym(np.ndarray out not None,
               np.ndarray A not None):
    cdef int32 ret
    cdef FMField[1] _out, _A

    array2fmfield4(_out, out)
    array2fmfield4(_A, A)

    ret = _sym2nonsym(_out, _A)
    return ret

def dw_adj_convect1(np.ndarray out not None,
                    np.ndarray state_w not None,
                    np.ndarray grad_u not None,
                    CMapping cmap not None,
                    int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _state_w, _grad_u

    array2fmfield4(_out, out)
    array2fmfield4(_state_w, state_w)
    array2fmfield4(_grad_u, grad_u)

    ret = _dw_adj_convect1(_out, _state_w, _grad_u, cmap.geo, is_diff)
    return ret

def dw_adj_convect2(np.ndarray out not None,
                    np.ndarray state_w not None,
                    np.ndarray state_u not None,
                    CMapping cmap not None,
                    int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _state_w, _state_u

    array2fmfield4(_out, out)
    array2fmfield4(_state_w, state_w)
    array2fmfield4(_state_u, state_u)

    ret = _dw_adj_convect2(_out, _state_w, _state_u, cmap.geo, is_diff)
    return ret

def dw_st_adj_supg_c(np.ndarray out not None,
                     np.ndarray state_w not None,
                     np.ndarray state_u not None,
                     np.ndarray grad_u not None,
                     np.ndarray coef not None,
                     CMapping cmap not None,
                     np.ndarray conn not None,
                     int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _state_w, _state_u, _grad_u, _coef
    cdef int32 *_conn
    cdef int32 n_el, n_ep

    array2fmfield4(_out, out)
    array2fmfield1(_state_w, state_w)
    array2fmfield4(_state_u, state_u)
    array2fmfield4(_grad_u, grad_u)
    array2fmfield4(_coef, coef)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    ret = _dw_st_adj_supg_c(_out, _state_w, _state_u, _grad_u, _coef,
                            cmap.geo, _conn, n_el, n_ep, is_diff)
    return ret

def dw_st_adj1_supg_p(np.ndarray out not None,
                      np.ndarray state_w not None,
                      np.ndarray grad_p not None,
                      np.ndarray coef not None,
                      CMapping cmap_w not None,
                      np.ndarray conn_w not None,
                      int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _state_w, _grad_p, _coef
    cdef int32 *_conn_w
    cdef int32 n_el_w, n_ep_w

    array2fmfield4(_out, out)
    array2fmfield1(_state_w, state_w)
    array2fmfield4(_grad_p, grad_p)
    array2fmfield4(_coef, coef)
    array2pint2(&_conn_w, &n_el_w, &n_ep_w, conn_w)

    ret = _dw_st_adj1_supg_p(_out, _state_w, _grad_p, _coef,
                             cmap_w.geo, _conn_w, n_el_w, n_ep_w, is_diff)
    return ret

def dw_st_adj2_supg_p(np.ndarray out not None,
                      np.ndarray grad_u not None,
                      np.ndarray state_r not None,
                      np.ndarray coef not None,
                      CMapping cmap_u not None,
                      CMapping cmap_r not None,
                      np.ndarray conn_r not None,
                      int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _grad_u, _state_r, _coef
    cdef int32 *_conn_r
    cdef int32 n_el_r, n_ep_r

    array2fmfield4(_out, out)
    array2fmfield4(_grad_u, grad_u)
    array2fmfield1(_state_r, state_r)
    array2fmfield4(_coef, coef)
    array2pint2(&_conn_r, &n_el_r, &n_ep_r, conn_r)

    ret = _dw_st_adj2_supg_p(_out, _grad_u, _state_r, _coef,
                             cmap_u.geo, cmap_r.geo, _conn_r, n_el_r, n_ep_r,
                             is_diff)
    return ret

def d_of_nsMinGrad(np.ndarray out not None,
                   np.ndarray grad not None,
                   np.ndarray viscosity not None,
                   CMapping cmap not None):
    cdef int32 ret
    cdef FMField[1] _out, _grad, _viscosity

    array2fmfield4(_out, out)
    array2fmfield4(_grad, grad)
    array2fmfield4(_viscosity, viscosity)

    ret = _d_of_nsMinGrad(_out, _grad, _viscosity, cmap.geo)
    return ret

def d_of_nsSurfMinDPress(np.ndarray out not None,
                         np.ndarray pressure not None,
                         float64 weight,
                         float64 bpress,
                         CMapping cmap not None,
                         int32 is_diff):
    cdef int32 ret
    cdef FMField[1] _out, _pressure

    array2fmfield4(_out, out)
    array2fmfield4(_pressure, pressure)

    ret = _d_of_nsSurfMinDPress(_out, _pressure, weight, bpress,
                                cmap.geo, is_diff)
    return ret

def d_sd_div(np.ndarray out not None,
             np.ndarray div_u not None,
             np.ndarray grad_u not None,
             np.ndarray state_p not None,
             np.ndarray div_mv not None,
             np.ndarray grad_mv not None,
             CMapping cmap_u not None,
             int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _state_p, _div_u, _grad_u
    cdef FMField[1] _div_mv, _grad_mv

    array2fmfield4(_out, out)
    array2fmfield4(_div_u, div_u)
    array2fmfield4(_grad_u, grad_u)
    array2fmfield4(_state_p, state_p)
    array2fmfield4(_div_mv, div_mv)
    array2fmfield4(_grad_mv, grad_mv)

    ret = _d_sd_div(_out, _div_u, _grad_u, _state_p, _div_mv, _grad_mv,
                    cmap_u.geo, mode)
    return ret

def d_sd_div_grad(np.ndarray out not None,
                  np.ndarray grad_u not None,
                  np.ndarray grad_w not None,
                  np.ndarray div_mv not None,
                  np.ndarray grad_mv not None,
                  np.ndarray viscosity not None,
                  CMapping cmap_u not None,
                  int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _grad_u, _grad_w, _div_mv, _grad_mv
    cdef FMField[1] _viscosity

    array2fmfield4(_out, out)
    array2fmfield4(_grad_u, grad_u)
    array2fmfield4(_grad_w, grad_w)
    array2fmfield4(_div_mv, div_mv)
    array2fmfield4(_grad_mv, grad_mv)
    array2fmfield4(_viscosity, viscosity)

    ret = _d_sd_div_grad(_out, _grad_u, _grad_w, _div_mv, _grad_mv, _viscosity,
                         cmap_u.geo, mode)
    return ret

def d_sd_convect(np.ndarray out not None,
                 np.ndarray state_u not None,
                 np.ndarray grad_u not None,
                 np.ndarray state_w not None,
                 np.ndarray div_mv not None,
                 np.ndarray grad_mv not None,
                 CMapping cmap_u not None,
                 int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _state_u, _grad_u, _state_w
    cdef FMField[1] _div_mv, _grad_mv

    array2fmfield4(_out, out)
    array2fmfield4(_state_u, state_u)
    array2fmfield4(_grad_u, grad_u)
    array2fmfield4(_state_w, state_w)
    array2fmfield4(_div_mv, div_mv)
    array2fmfield4(_grad_mv, grad_mv)

    ret = _d_sd_convect(_out, _state_u, _grad_u, _state_w, _div_mv, _grad_mv,
                        cmap_u.geo, mode)
    return ret

def d_sd_volume_dot(np.ndarray out not None,
                    np.ndarray state_p not None,
                    np.ndarray state_q not None,
                    np.ndarray div_mv not None,
                    CMapping cmap not None,
                    int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _state_p, _state_q, _div_mv

    array2fmfield4(_out, out)
    array2fmfield4(_state_p, state_p)
    array2fmfield4(_state_q, state_q)
    array2fmfield4(_div_mv, div_mv)

    ret = _d_sd_volume_dot(_out, _state_p, _state_q, _div_mv, cmap.geo, mode)
    return ret

def d_sd_st_grad_div(np.ndarray out not None,
                     np.ndarray div_u not None,
                     np.ndarray grad_u not None,
                     np.ndarray div_w not None,
                     np.ndarray grad_w not None,
                     np.ndarray div_mv not None,
                     np.ndarray grad_mv not None,
                     np.ndarray coef not None,
                     CMapping cmap_u not None,
                     int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _div_u, _grad_u, _div_w, _grad_w
    cdef FMField[1] _div_mv, _grad_mv, _coef

    array2fmfield4(_out, out)
    array2fmfield4(_div_u, div_u)
    array2fmfield4(_grad_u, grad_u)
    array2fmfield4(_div_w, div_w)
    array2fmfield4(_grad_w, grad_w)
    array2fmfield4(_div_mv, div_mv)
    array2fmfield4(_grad_mv, grad_mv)
    array2fmfield4(_coef, coef)

    ret = _d_sd_st_grad_div(_out, _div_u, _grad_u, _div_w, _grad_w,
                            _div_mv, _grad_mv, _coef, cmap_u.geo, mode)
    return ret

def d_sd_st_supg_c(np.ndarray out not None,
                   np.ndarray state_b not None,
                   np.ndarray grad_u not None,
                   np.ndarray grad_w not None,
                   np.ndarray div_mv not None,
                   np.ndarray grad_mv not None,
                   np.ndarray coef not None,
                   CMapping cmap_u not None,
                   int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _state_b, _grad_u, _grad_w
    cdef FMField[1] _div_mv, _grad_mv, _coef

    array2fmfield4(_out, out)
    array2fmfield4(_state_b, state_b)
    array2fmfield4(_grad_u, grad_u)
    array2fmfield4(_grad_w, grad_w)
    array2fmfield4(_div_mv, div_mv)
    array2fmfield4(_grad_mv, grad_mv)
    array2fmfield4(_coef, coef)

    ret = _d_sd_st_supg_c(_out, _state_b, _grad_u, _grad_w, _div_mv, _grad_mv,
                          _coef, cmap_u.geo, mode)
    return ret

def d_sd_st_pspg_c(np.ndarray out not None,
                   np.ndarray state_b not None,
                   np.ndarray grad_u not None,
                   np.ndarray grad_r not None,
                   np.ndarray div_mv not None,
                   np.ndarray grad_mv not None,
                   np.ndarray coef not None,
                   CMapping cmap_u not None,
                   int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _state_b, _grad_u, _grad_r
    cdef FMField[1] _div_mv, _grad_mv, _coef

    array2fmfield4(_out, out)
    array2fmfield4(_state_b, state_b)
    array2fmfield4(_grad_u, grad_u)
    array2fmfield4(_grad_r, grad_r)
    array2fmfield4(_div_mv, div_mv)
    array2fmfield4(_grad_mv, grad_mv)
    array2fmfield4(_coef, coef)

    ret = _d_sd_st_pspg_c(_out, _state_b, _grad_u, _grad_r, _div_mv, _grad_mv,
                          _coef, cmap_u.geo, mode)
    return ret

def d_sd_st_pspg_p(np.ndarray out not None,
                   np.ndarray grad_r not None,
                   np.ndarray grad_p not None,
                   np.ndarray div_mv not None,
                   np.ndarray grad_mv not None,
                   np.ndarray coef not None,
                   CMapping cmap_p not None,
                   int32 mode):
    cdef int32 ret
    cdef FMField[1] _out, _grad_r, _grad_p
    cdef FMField[1] _div_mv, _grad_mv, _coef

    array2fmfield4(_out, out)
    array2fmfield4(_grad_r, grad_r)
    array2fmfield4(_grad_p, grad_p)
    array2fmfield4(_div_mv, div_mv)
    array2fmfield4(_grad_mv, grad_mv)
    array2fmfield4(_coef, coef)

    ret = _d_sd_st_pspg_p(_out, _grad_r, _grad_p, _div_mv, _grad_mv, _coef,
                          cmap_p.geo, mode)
    return ret
