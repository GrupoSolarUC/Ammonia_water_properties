# -*- coding: utf-8 -*-
"""
Created on Tue May 30 12:14:03 2023

@author: adria
"""

from math import log, exp
from inverse import inversefunc
from bisect import bisect

class AmmoniaFunctionError(Exception):
    def __init__(self, message):
        super().__init__(message)


constants_1 = {1: {'a': 0.322302e1, 'm': 0, 'n': 0},
               2: {'a': -0.384206e0, 'm': 0, 'n': 1},
               3: {'a': 0.460965e-1, 'm': 0, 'n': 2},
               4: {'a': -0.378945e-2, 'm': 0, 'n': 3},
               5: {'a': 0.135610e-3, 'm': 0, 'n': 4},
               6: {'a': 0.487755e0, 'm': 1, 'n': 0},
               7: {'a': -0.120108e0, 'm': 1, 'n': 1},
               8: {'a': 0.106154e-1, 'm': 1, 'n': 2},
               9: {'a': -0.533589e-3, 'm': 2, 'n': 3},
               10: {'a': 0.785041e1, 'm': 4, 'n': 0},
               11: {'a': -0.115941e2, 'm': 5, 'n': 0},
               12: {'a': -0.523150e-1, 'm': 5, 'n': 1},
               13: {'a': 0.489596e1, 'm': 6, 'n': 0},
               14: {'a': 0.421059e-1, 'm': 13, 'n': 1},
               'p0': 2, 'T0': 100}

constants_2 = {1: {'a': 0.324004e1, 'm': 0, 'n': 0},
               2: {'a': -0.395920e0, 'm': 0, 'n': 1},
               3: {'a': 0.435624e-1, 'm': 0, 'n': 2},
               4: {'a': -0.218943e-2, 'm': 0, 'n': 3},
               5: {'a': -0.143526e1, 'm': 1, 'n': 0},
               6: {'a': 0.105256e1, 'm': 1, 'n': 1},
               7: {'a': -0.719281e-1, 'm': 1, 'n': 2},
               8: {'a': 0.122362e2, 'm': 2, 'n': 0},
               9: {'a': -0.224368e1, 'm': 2, 'n': 1},
               10: {'a': -0.201780e2, 'm': 3, 'n': 0},
               11: {'a': 0.110834e1, 'm': 3, 'n': 1},
               12: {'a': 0.145399e2, 'm': 4, 'n': 0},
               13: {'a': 0.644312e0, 'm': 4, 'n': 2},
               14: {'a': -0.221246e1, 'm': 5, 'n': 0},
               15: {'a': -0.756266e0, 'm': 5, 'n': 2},
               16: {'a': -0.135529e1, 'm': 6, 'n': 0},
               17: {'a': 0.183541e0, 'm': 7, 'n': 2},
               'p0': 2, 'T0': 100}

constants_3 = {1: {'a': 1.98022017e1, 'm': 0, 'n': 0},
               2: {'a': -1.18092669e1, 'm': 0, 'n': 1},
               3: {'a': 2.77479980e1, 'm': 0, 'n': 6},
               4: {'a': -2.88634277e1, 'm': 0, 'n': 7},
               5: {'a': -5.91616608e1, 'm': 1, 'n': 0},
               6: {'a': 5.78091305e2, 'm': 2, 'n': 1},
               7: {'a': -6.21736743e0, 'm': 2, 'n': 2},
               8: {'a': -3.42198402e3, 'm': 3, 'n': 2},
               9: {'a': 1.19403127e4, 'm': 4, 'n': 3},
               10: {'a': -2.45413777e4, 'm': 5, 'n': 4},
               11: {'a': 2.91591865e4, 'm': 6, 'n': 5},
               12: {'a': -1.84782290e4, 'm': 7, 'n': 6},
               13: {'a': 2.34819434e1, 'm': 7, 'n': 7},
               14: {'a': 4.80310617e3, 'm': 8, 'n': 7},
               'p0': 2}

constants_4 = {1: {'a': -0.761080e1, 'm': 0, 'n': 1},
               2: {'a': 0.256905e2, 'm': 0, 'n': 4},
               3: {'a': -0.247092e3, 'm': 0, 'n': 8},
               4: {'a': 0.325952e3, 'm': 0, 'n': 9},
               5: {'a': -0.158854e3, 'm': 0, 'n': 12},
               6: {'a': 0.619084e2, 'm': 0, 'n': 14},
               7: {'a': 0.114314e2, 'm': 1, 'n': 0},
               8: {'a': 0.118157e1, 'm': 1, 'n': 1},
               9: {'a': 0.284179e1, 'm': 2, 'n': 1},
               10: {'a': 0.741609e1, 'm': 3, 'n': 3},
               11: {'a': 0.891844e3, 'm': 5, 'n': 3},
               12: {'a': -0.161309e4, 'm': 5, 'n': 4},
               13: {'a': 0.622106e3, 'm': 5, 'n': 5},
               14: {'a': -0.207588e3, 'm': 6, 'n': 2},
               15: {'a': -0.687393e1, 'm': 6, 'n': 4},
               16: {'a': 0.350716e1, 'm': 8, 'n': 0},
               'h0': 100, 'T0': 273.16}

constants_5 = {1: {'a': 0.128827e1, 'm': 0, 'n': 0},
               2: {'a': 0.125247e0, 'm': 1, 'n': 0},
               3: {'a': -0.208748e1, 'm': 2, 'n': 0},
               4: {'a': 0.217696e1, 'm': 3, 'n': 0},
               5: {'a': 0.235687e1, 'm': 0, 'n': 2},
               6: {'a': -0.886987e1, 'm': 1, 'n': 2},
               7: {'a': 0.102635e2, 'm': 2, 'n': 2},
               8: {'a': -0.237440e1, 'm': 3, 'n': 2},
               9: {'a': -0.670515e1, 'm': 0, 'n': 3},
               10: {'a': 0.164508e2, 'm': 1, 'n': 3},
               11: {'a': -0.936849e1, 'm': 2, 'n': 3},
               12: {'a': 0.842254e1, 'm': 0, 'n': 4},
               13: {'a': -0.858807e1, 'm': 1, 'n': 4},
               14: {'a': -0.277049e1, 'm': 0, 'n': 5},
               15: {'a': -0.961248e0, 'm': 4, 'n': 6},
               16: {'a': 0.988009e0, 'm': 2, 'n': 7},
               17: {'a': 0.308482e0, 'm': 1, 'n': 10},
               'h0': 1000, 'T0': 324}


def T_from_px(p,x):
    """
    Función para obtener la temperatura (en K) como función de la presión (p, en MPa) y la fracción molar de amoniaco en la fase líquida (x)
    """
    return constants_1['T0']*sum([ constants_1[i]['a']*(1-x)**constants_1[i]['m']*log(constants_1['p0']/p)**constants_1[i]['n'] for i in range(1,15) ])

def T_from_py(p,y):
    """
    Función para obtener la temperatura (en K) como función de la presión (p, en MPa) y la fracción molar de amoniaco en la fase gaseosa (y)
    """
    return constants_2['T0']*sum([ constants_2[i]['a']*(1-y)**(constants_2[i]['m']/4)*log(constants_2['p0']/p)**constants_2[i]['n'] for i in range(1,18) ])

def y_from_px(p,x):
    """
    Función para obtener la fracción molar de amoniaco en la fase gaseosa como función de la presión (p, en MPa) y la fracción molar de amoniaco en la fase líquida (x)
    """
    summation = sum([ constants_3[i]['a']*(p/constants_3['p0'])**constants_3[i]['m']*x**(constants_3[i]['n']/3) for i in range(1,15) ])
    return 1 - exp(log(1-x)*summation)

def Hl_from_Tx(T,x):
    """
    Función para obtener la entalpía de la fase líquida (en kJ/kg) como función de la presión (p, en MPa) y la fracción molar de amoniaco en la fase líquida (x)
    """
    return constants_4['h0']*sum([ constants_4[i]['a']*(T/constants_4['T0'] - 1)**constants_4[i]['m']*x**constants_4[i]['n'] for i in range(1,17) ])

def Hg_from_Ty(T,y):
    """
    Función para obtener la entalpía de la fase gaseosa (en kJ/kg) como función de la presión (p, en MPa) y la fracción molar de amoniaco en la fase gaseosa (y)
    """
    return constants_5['h0']*sum([ constants_5[i]['a']*(1 - T/constants_5['T0'])**constants_5[i]['m']*(1 - y)**(constants_5[i]['n']/4) for i in range(1,18) ])

def molar_to_massic_quality(q):
    """
    Función para convertir la fracción molar de amoniaco en fracción másica
    """
    MA = 17.031
    MW = 18.015
    return q*MA/(q*MA + (1-q)*MW)
    
def massic_to_molar_quality(q):
    """
    Función para convertir la fracción másica de amoniaco en fracción molar
    """
    MA = 17.031
    MW = 18.015
    return q*MW/(q*MW + (1-q)*MA)

def x_from_py(p,y):
    """
    Función para obtener la fracción molar de amoniaco en la fase líquida como función de la presión (p, en MPa) y la fracción molar de amoniaco en la fase gaseosa (y)
    """
    def y_from_x(x):
        return y_from_px(p,x)
    return inversefunc(y_from_x, y_values = [y], domain = [0,0.99999])[0]

def p_from_Tx(T, x, domain):
    """
    Función para obtener la presión (en MPa) como función de la temperatura (T, en K) y la fracción molar de amoniaco en la fase líquida (x)
    """
    def T_from_p(p):
        return T_from_px(p, x)
    return inversefunc(T_from_p, y_values = [T], domain = domain)[0]

def p_from_Ty(T, y, domain):
    """
    Función para obtener la presión (en MPa) como función de la temperatura (T, en K) y la fracción molar de amoniaco en la fase gaseosa (y)
    Se debe especificar un dominio (domain) para la presión en forma de lista: [p_min, p_max]. En este caso, el dominio es entre 0.002 MPa y 2 MPa
    """
    def T_from_p(p):
        return T_from_py(p, y)
    return inversefunc(T_from_p, y_values = [T], domain = domain)[0]

def p_from_xy(x,y):
    """
    Función para obtener la presión (en MPa) como función de la fracción molar de amoniaco en la fase líquida (x) y la fracción molar de amoniaco en la fase gaseosa (y)
    """
    T_min = T_from_px(1e-3, x)
    T_max = T_from_px(3,x)
    def T_from_p(p):
        return T_from_px(p, x)
    p_from_Tx = inversefunc(T_from_p, domain = [1e-3,3])
    T = 325
    it = 0
    while True:
        if T > T_max:
            T = T_max
        if T < T_min:
            T = T_min
        p = p_from_Tx(T)
        T_new = T_from_py(p, y)
        it = it + 1
        if abs(T - T_new) < 1e-5:
            break
        if it == 200:
            break
        T = T_new
    return p

def x_from_pT(p,T):
    """
    Función para obtener la fracción molar de amoniaco en la fase líquida como función de la presión (p, en MPa) y la temperatura (T, en K)
    """
    def T_from_x(x):
        return T_from_px(p, x)
    return inversefunc(T_from_x, y_values = [T], domain = [0,1])[0]

def Hl_from_px(p, x):
    """
    Función para obtener la entalpía de la fase líquida como función de la presión (p, en MPa) y la fracción molar de amoniaco en la fase líquida (x)
    """
    T = T_from_px(p, x)
    return Hl_from_Tx(T, x)

def p_from_xHl(x, Hl, domain):
    """
    Función para obtener la la presión (en MPa) como función de la fracción molar de amoniaco en la fase líquida (x) y la entalpía de la fase líquida (Hl, en kJ/kg)
    Se debe especificar un dominio (domain) para la presión en forma de lista: [p_min, p_max]. En este caso, el dominio es entre 0.002 MPa y 2 MPa
    """
    def Hl_from_p(p):
        return Hl_from_px(p, x)
    return inversefunc(Hl_from_p, y_values = [Hl], domain = domain)[0]

def Hg_from_py(p, y):
    """
    Función para obtener la entalpía de la fase gaseosa (en kJ/kg) como función de la presión (p, en MPa) y la fracción molar de amoniaco en la fase gaseosa
    """
    T = T_from_py(p, y)
    return Hg_from_Ty(T, y)

def p_from_yHg(y, Hg, domain):
    """
    Función para obtener la presión (en MPa) como función de la fracción molar de amoniaco en la fase gaseosa (y) y la entalpía de la fase gaseosa (Hg, en kJ/kg)
    Se debe especificar un dominio (domain) para la presión en forma de lista: [p_min, p_max]. En este caso, el dominio es entre 0.002 MPa y 2 MPa
    """
    def Hg_from_p(p):
        return Hg_from_py(p, y)
    return inversefunc(Hg_from_p, y_values = [Hg], domain = domain)[0]

def y_from_pT(p, T):
    """
    Función para obtener la fracción molar de amoniaco en la fase gaseosa como función de la presión (p, en MPa) y la temperatura (T, en K)
    """
    def T_from_y(y):
        return T_from_py(p, y)
    return inversefunc(T_from_y, y_values = [T], domain = [0,1])[0]

def y_from_HgT(Hg, T):
    """
    Función para obtener la fracción molar de amoniaco en la fase gaseosa como función de la entalpía de la fase gaseosa (Hg, en kJ/kg) y la temperatura (T, en K)
    """
    def Hg_from_y(y):
        return Hg_from_Ty(T, y)
    return inversefunc(Hg_from_y, y_values = [Hg], domain = [0,1])[0]

def y_from_pHg(p, Hg):
    """
    Función para obtener la fracción molar de amoniaco en la fase gaseosa como función de la presión (p, en MPa) y la entalpía de la fase gaseosa (Hg, en kJ/kg)
    """
    def Hg_from_y(y):
        return Hg_from_py(p, y)
    return inversefunc(Hg_from_y, y_values = [Hg], domain = [0,1])[0]

def Hl_from_py(p, y):
    """
    Función para obtener la entalpía de la fase líquida (en kJ/kg) como función de la presión (p, en MPa) y la fracción molar de amoniaco en la fase gaseosa (y)
    """
    x = x_from_py(p, y)
    return Hl_from_px(p, x)

def p_from_yHl(y, Hl):
    """
    Función para obtener la presión (en MPa) como función de la fracción molar de amoniaco en la fase gaseosa (y) y la entalpía de la fase líquida (Hl, en kJ/kg)
    """
    def Hl_from_p(p):
        return Hl_from_py(p, y)
    return inversefunc(Hl_from_p, y_values = [Hl], domain = [0.002, 2])[0]

def Hg_from_px(p, x):
    """
    Función para obtener la entalpía de la fase gaseosa (en kJ/kg) como función de la presión (p, en MPa) y la fracción molar de amoniaco en la fase líquida (x)
    """
    y = y_from_px(p, x)
    return Hg_from_py(p, y)

def p_from_xHg(x, Hg):
    """
    Función para obtener la presión (en MPa) como función de la fracción molar de amoniaco en la fase líquida (x) y la entalpía de la fase gaseosa (Hg, en kJ/kg)
    """
    def Hg_from_p(p):
        return Hg_from_px(p, x)
    return inversefunc(Hg_from_p, y_values = [Hg], domain = [0.002,2])


#########################################################
#########################################################
#########################################################
#########################################################
#########################################################

## Función principal

def PropsSI_NH3_H2O(variable_required, variable_1_name, variable_1_value, variable_2_name, variable_2_value):
    """
    Función que entrega las propiedades de la mezcla de agua con amoniaco.

    Parámetros:
        - variable_required: Variable que será entregada por la función como resultado.
        - variable_1_name: Variable entregada como primer dato para obtener un resultado.
        - variable_1_value: Valor de la primera variable entregada como dato.
        - variable_2_name: Variable entregada como segundo dato para obtener un resultado.
        - variable_2_value: Valor de la segunda variable entregada como dato.

    Retorna:
        - Valor de la variable solicitada en el primer parámetro, como resultado de las dos variables entregadas como datos.
        
    Los parámetros 'variable_required', 'variable_1_name', y 'variable_2_name' pueden ser los siguientes:
        - 'p': Presión [Pa]
        - 'T': Temperatura [K]
        - 'Ql': Fracción másica de amoniaco en la fase líquida (entre 0 y 1)
        - 'Qg': Fracción másica de amoniaco en la fase gaseosa (entre 0 y 1)
        - 'Hl': Entalpía de la fase líquida [J/kg]
        - 'Hg': Entalpía de la fase gaseosa [J/kg]
        
        Los únicos pares de variables que no se pueden usar como variables de entrada son los siguientes:
            - 'p', 'Hl'
            - 'T', 'Hl'
            - 'Hl', 'Hg'
    """
    if not variable_required in ['T', 'p', 'Ql', 'Qg', 'Hl', 'Hg']:
        raise AmmoniaFunctionError("Los nombres de las variables deben ser 'T', 'p', 'Ql', 'Qg', 'Hl' o 'Hg'")
    if not variable_1_name in ['T', 'p', 'Ql', 'Qg', 'Hl', 'Hg']:
        raise AmmoniaFunctionError("Los nombres de las variables deben ser 'T', 'p', 'Ql', 'Qg', 'Hl' o 'Hg'")
    if not variable_2_name in ['T', 'p', 'Ql', 'Qg', 'Hl', 'Hg']:
        raise AmmoniaFunctionError("Los nombres de las variables deben ser 'T', 'p', 'Ql', 'Qg', 'Hl' o 'Hg'")
    if variable_1_name == variable_2_name or variable_1_name == variable_required or variable_2_name == variable_required:
        raise AmmoniaFunctionError("Las tres variables (la solicitada y las dos entregadas) deben ser diferentes")
    
    if (variable_1_name == 'T' and variable_2_name == 'Hl') or (variable_1_name == 'Hl' and variable_2_name == 'T'):
        raise AmmoniaFunctionError("No es posible combinar las variables de entrada 'T' y 'Hl'")
    if (variable_1_name == 'p' and variable_2_name == 'Hl') or (variable_1_name == 'Hl' and variable_2_name == 'p'):
        raise AmmoniaFunctionError("No es posible combinar las variables de entrada 'p' y 'Hl'")
    if (variable_1_name == 'Hl' and variable_2_name == 'Hg') or (variable_1_name == 'Hg' and variable_2_name == 'Hl'):
        raise AmmoniaFunctionError("No es posible combinar las variables de entrada 'Hl' y 'Hg'")

    ## Caso 1: p, Ql
    if (variable_1_name == 'p' and variable_2_name == 'Ql') or (variable_1_name == 'Ql' and variable_2_name == 'p'):
        if variable_1_name == 'p':
            p_orig = variable_1_value
            x = massic_to_molar_quality(variable_2_value)
        else:
            x = massic_to_molar_quality(variable_1_value)
            p_orig = variable_2_value
        p = p_orig/1e6
        if variable_required == 'Qg' or variable_required == 'Hg':
            if x < 0.05 or x >= 1:
                raise AmmoniaFunctionError('La fracción de amoniaco en la fase líquida debe ser mayor o igual a ' + str(molar_to_massic_quality(0.05)) + ' y menor a uno.')
            if p_orig < 20000 or p_orig > 2e6:
                raise AmmoniaFunctionError('La presión debe ser mayor o igual a 20 kPa y menor o igual a 2 MPa (especificar presión en Pa).')
            y = y_from_px(p, x)
            if variable_required == 'Qg':
                return molar_to_massic_quality(y)
            if variable_required == 'Hg':
                return 1000*Hg_from_py(p, y)
        if variable_required == 'T' or variable_required == 'Hl':
            if p_orig < 2000 or p_orig > 2e6:
                raise AmmoniaFunctionError('La presión debe ser mayor o igual a 2 kPa y menor o igual a 2 MPa (especificar presión en Pa).')
            if x <= 0 or x >= 1:
                raise AmmoniaFunctionError('La fracción de amoniaco en la fase líquida debe ser mayor a cero y menor a uno.')
            if variable_required == 'T':
                return T_from_px(p, x)
            if variable_required == 'Hl':
                return 1000*Hl_from_px(p, x)
    
    ## Caso 2: p, Qg
    if (variable_1_name == 'p' and variable_2_name == 'Qg') or (variable_1_name == 'Qg' and variable_2_name == 'p'):
        if variable_1_name == 'p':
            p_orig = variable_1_value
            gas_phase_q = variable_2_value
        else:
            gas_phase_q = variable_1_value
            p_orig = variable_2_value
        y = massic_to_molar_quality(gas_phase_q)
        p = p_orig/1e6
        if p_orig < 20000 or p_orig > 2e6:
            raise AmmoniaFunctionError('La presión debe ser mayor o igual a 20 kPa y menor o igual a 2 MPa (especificar presión en Pa).')
        if variable_required == 'Hg' or variable_required == 'T':
            if y <= 0 or y >= 1:
                raise AmmoniaFunctionError('La fracción de amoniaco en la fase gaseosa debe ser mayor a cero y menor a uno.')
            if variable_required == 'Hg':
                return 1000*Hg_from_py(p, y)
            if variable_required == 'T':
                return T_from_py(p, y)
        if variable_required == 'Ql' or variable_required == 'Hl':
            y_min = y_from_px(p, 0.05)
            if y < y_min or y >= 1:
                raise AmmoniaFunctionError('Para una presión de ' + str(p_orig) + ' Pa, la fracción de amoniaco en la fase gaseosa debe ser mayor o igual a ' + str(molar_to_massic_quality(y_min)) + ' y menor a uno.')
            x = x_from_py(p, y)
            if variable_required == 'Ql':
                return molar_to_massic_quality(x)
            if variable_required == 'Hl':
                return 1000*Hl_from_px(p, x)
    
    ## Caso 3: T, Ql
    if (variable_1_name == 'T' and variable_2_name == 'Ql') or (variable_1_name == 'Ql' and variable_2_name == 'T'):
        if variable_1_name == 'T':
            T = variable_1_value
            liquid_phase_q = variable_2_value
        else:
            liquid_phase_q = variable_1_value
            T = variable_2_value
        x = massic_to_molar_quality(liquid_phase_q)
        if variable_required == 'p' or variable_required == 'Hl':
            if x <= 0 or x >= 1:
                raise AmmoniaFunctionError('La fracción de amoniaco en la fase líquida debe ser mayor a cero y menor a uno.')
            T_min = T_from_px(0.002, x)
            T_max = T_from_px(2, x)
            if T < T_min or T > T_max:
                raise AmmoniaFunctionError('Para una fracción de amoniaco en la fase líquida igual a ' + str(liquid_phase_q) + ', la temperatura debe ser mayor o igual a ' + str (T_min) + ' K y menor o igual a ' + str(T_max) + ' K.')
            if variable_required == 'Hl':
                return 1000*Hl_from_Tx(T, x)
            if variable_required == 'p':
                p = p_from_Tx(T, x, [0.002,2])
                return 1e6*p
        if variable_required == 'Qg' or variable_required == 'Hg':
            if x < 0.05 or x >= 1:
                raise AmmoniaFunctionError('La fracción de amoniaco en la fase líquida debe ser mayor o igual a ' + str(molar_to_massic_quality(0.05)) + ' y menor a uno.')
            T_min = T_from_px(0.02, x)
            T_max = T_from_px(2, x)
            if T < T_min or T > T_max:
                raise AmmoniaFunctionError('Para una fracción de amoniaco en la fase líquida igual a ' + str(liquid_phase_q) + ', la temperatura debe ser mayor o igual a ' + str (T_min) + ' K y menor o igual a ' + str(T_max) + ' K.')
            p = p_from_Tx(T, x, [0.02,2])
            y = y_from_px(p, x)
            if variable_required == 'Qg':
                return molar_to_massic_quality(y)
            if variable_required == 'Hg':
                return 1000*Hg_from_Ty(T, y)
    
    ## Caso 4: T, Qg  (FALTA RESTRICCIÓN PARA QUE x > 0.05)
    if (variable_1_name == 'T' and variable_2_name == 'Qg') or (variable_1_name == 'Qg' and variable_2_name == 'T'):
        if variable_1_name == 'T':
            T = variable_1_value
            gas_phase_q = variable_2_value
        else:
            gas_phase_q = variable_1_value
            T = variable_2_value
        y = massic_to_molar_quality(gas_phase_q)
        if y <= 0 or y >= 1:
            raise AmmoniaFunctionError('La fracción de amoniaco en la fase gaseosa debe ser mayor a cero y menor a uno.')
        T_min = T_from_py(0.02, y)
        T_max = T_from_py(2, y)
        if T < T_min or T > T_max:
            raise AmmoniaFunctionError('Para una fracción de amoniaco en la fase gaseosa igual a ' + str(gas_phase_q) + ', la temperatura debe ser mayor o igual a ' + str (T_min) + ' K y menor o igual a ' + str(T_max) + 'K.')
        if variable_required == 'Hg':
            return 1000*Hg_from_Ty(T, y)
        p = p_from_Ty(T, y, [0.02,2])
        if variable_required == 'p':
            return 1e6*p
        x = x_from_py(p, y)
        if variable_required == 'Ql':
            return molar_to_massic_quality(x)
        if variable_required == 'Hl':
            return 1000*Hl_from_Tx(T, x)
    
    ## Caso 5: T, p
    if (variable_1_name == 'T' and variable_2_name == 'p') or (variable_1_name == 'p' and variable_2_name == 'T'):
        if variable_1_name == 'T':
            T = variable_1_value
            p_orig = variable_2_value
        else:
            p_orig = variable_1_value
            T = variable_2_value
        p = p_orig/1e6
        if variable_required == 'Ql' or variable_required == 'Hl':
            if p_orig < 2000 or p_orig > 2e6:
                raise AmmoniaFunctionError('La presión debe ser mayor o igual a 2 kPa y menor o igual a 2 MPa (especificar presión en Pa).')
            T_min = T_from_px(p, 1)
            T_max = T_from_px(p, 0)
            if T <= T_min or T >= T_max:
                raise AmmoniaFunctionError('Para una presión de '+ str(p_orig)+' Pa, la temperatura debe ser mayor a ' + str (T_min) + ' K y menor a ' + str(T_max) + 'K.')
            x = x_from_pT(p, T)
            if variable_required == 'Ql':
                return molar_to_massic_quality(x)
            if variable_required == 'Hl':
                return 1000*Hl_from_Tx(T, x)
        if variable_required == 'Qg' or variable_required == 'Hg':
            if p_orig < 20000 or p_orig > 2e6:
                raise AmmoniaFunctionError('La presión debe ser mayor o igual a 20 kPa y menor o igual a 2 MPa (especificar presión en Pa).')
            T_min = T_from_py(p, 1)
            T_max = T_from_py(p, 0)
            if T <= T_min or T >= T_max:
                raise AmmoniaFunctionError('Para una presión de '+ str(p_orig)+' Pa, la temperatura debe ser mayor a ' + str (T_min) + ' K y menor a ' + str(T_max) + 'K.')
            y = y_from_pT(p, T)
            if variable_required == 'Qg':
                return molar_to_massic_quality(y)
            if variable_required == 'Hg':
                return 1000*Hg_from_Ty(T, y)
        
    ## Caso 6: Ql, Qg
    if (variable_1_name == 'Ql' and variable_2_name == 'Qg') or (variable_1_name == 'Qg' and variable_2_name == 'Ql'):
        
        y_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
                  0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6,
                  0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8,
                  0.8125, 0.825, 0.8375, 0.85, 0.86, 0.87, 0.88, 0.89,
                  0.9, 0.90833, 0.9166, 0.925, 0.93125, 0.9375, 0.94375, 0.95,
                  0.955, 0.96, 0.965, 0.97, 0.975, 0.9775, 0.98, 0.9825,
                  0.985, 0.9875, 0.99, 0.991, 0.992, 0.993, 0.994, 0.995,
                  0.9955, 0.996, 0.9965, 0.997, 0.9975, 0.998, 0.9985, 0.99875,
                  0.999, 0.9991, 0.9992, 0.9993, 0.9994, 0.9995, 0.999684, 0.9998,
                  0.9999, 0.99993, 0.99996, 0.999975, 0.99999, 0.999993, 0.999996, 0.999999]
        
        x_min_list = [0.00415828, 0.00725252, 0.01033296, 0.01347837, 0.01669839, 0.02002784,
                      0.02354702, 0.02736104, 0.02941827, 0.03160071, 0.03393017, 0.0364314,
                      0.03913267, 0.04206647, 0.04527043, 0.04878845, 0.05267221, 0.05698309,
                      0.06179474, 0.06719645, 0.07329774, 0.08023472, 0.08817903, 0.0973508,
                      0.10248374, 0.10803867, 0.11406765, 0.12063236, 0.12631911, 0.13244112,
                      0.13905236, 0.14621829, 0.15402011, 0.16107822, 0.16866791, 0.17707067,
                      0.18385457, 0.19117267, 0.19912155, 0.20783184, 0.21546905, 0.22384047,
                      0.23313058, 0.24361071, 0.25570335, 0.26255851, 0.2701192,  0.27857266,
                      0.28819303, 0.29940467, 0.31291589, 0.31922467, 0.32622868, 0.33411111,
                      0.34313863, 0.35372146, 0.35979203, 0.36653998, 0.37414155, 0.38285265,
                      0.393066,   0.40543138, 0.42114761, 0.43097219, 0.44284865, 0.44839853,
                      0.45455822, 0.46148391, 0.46940206, 0.47865961, 0.50143601, 0.52339064,
                      0.55524431, 0.57098802, 0.5948407,  0.61411055, 0.64982971, 0.66312481,
                      0.68335726, 0.73053364]
        
        x_max_list = [0.01460431, 0.02132113, 0.02877275, 0.03701659, 0.04611661, 0.05614862,
                      0.06719805, 0.07936689, 0.08590799, 0.09277598, 0.0999915,  0.10757524,
                      0.11555182, 0.12395138, 0.13280512, 0.14214381, 0.15203874, 0.16251921,
                      0.17366353, 0.18555366, 0.19829932, 0.21202813, 0.22696042, 0.24331502,
                      0.25213723, 0.26146452, 0.27136816, 0.2819426,  0.29095686, 0.30054636,
                      0.31080306, 0.32184032, 0.3338095,  0.34463222, 0.35630052, 0.36929357,
                      0.37987036, 0.39138867, 0.40405196, 0.41813196, 0.43067147, 0.44463996,
                      0.46042853, 0.47861361, 0.50009731, 0.5125144,  0.52640873, 0.54218829,
                      0.56046596, 0.58218946, 0.60899899, 0.62176234, 0.63611891, 0.65252412,
                      0.67165051, 0.69455489, 0.70793468, 0.72302475, 0.74028157, 0.76039398,
                      0.78435571, 0.81372392, 0.85099377, 0.87375683, 0.90015824, 0.91193203,
                      0.92446168, 0.93791392, 0.95239827, 0.96797299, 0.99996709, 0.99995936,
                      0.99996557, 0.99964257, 0.99966114, 0.99996852, 0.99998729, 0.99997642,
                      0.99997124, 0.99997188]
        
        y_mass_lim_inf = 0.04739839027485557
        y_mass_lim_sup = 0.9999989422230661
        
        if variable_1_name == 'Ql':
            liquid_phase_q = variable_1_value
            gas_phase_q = variable_2_value
        else:
            gas_phase_q = variable_1_value
            liquid_phase_q = variable_2_value
        if gas_phase_q < y_mass_lim_inf or gas_phase_q > y_mass_lim_sup:
            raise AmmoniaFunctionError('La fracción de amoniaco en la fase gaseosa debe ser mayor o igual a ' + str(y_mass_lim_inf) + ' y menor o igual a ' + str(y_mass_lim_sup) + '.')
        y = massic_to_molar_quality(gas_phase_q)
        x = massic_to_molar_quality(liquid_phase_q)
        index = bisect(y_list, y)
        x_min = x_min_list[index - 1] + ((y - y_list[index - 1])/(y_list[index] - y_list[index - 1]))*(x_min_list[index] - x_min_list[index - 1])
        x_max = x_max_list[index - 1] + ((y - y_list[index - 1])/(y_list[index] - y_list[index - 1]))*(x_max_list[index] - x_max_list[index - 1])
        if x < x_min or x > x_max:
            raise AmmoniaFunctionError('Para una fracción de amoniaco en la fase gaseosa igual a ' + str(gas_phase_q) + ', la fracción de amoniaco en la fase líquida debe ser mayor o igual a '+str(molar_to_massic_quality(x_min))+' y menor o igual a '+ str(molar_to_massic_quality(x_max))+'.')
        p = p_from_xy(x, y)
        if variable_required == 'p':
            return 1e6*p
        if variable_required == 'T':
            return T_from_px(p, x)
        if variable_required == 'Hl':
            return 1000*Hl_from_px(p, x)
        if variable_required == 'Hg':
            return 1000*Hg_from_py(p, y)
        
    ## Caso 7: Hl, Ql
    if (variable_1_name == 'Hl' and variable_2_name == 'Ql') or (variable_1_name == 'Ql' and variable_2_name == 'Hl'):
        if variable_1_name == 'Hl':
            Hl = variable_1_value/1000
            liquid_phase_q = variable_2_value
        else:
            liquid_phase_q = variable_1_value
            Hl = variable_2_value/1000
        x = massic_to_molar_quality(liquid_phase_q)
        if variable_required == 'p' or variable_required == 'T':
            if x <= 0 or x >= 1:
                raise AmmoniaFunctionError('La fracción de amoniaco en la fase líquida debe ser mayor a cero y menor a uno.')
            Hl_min = Hl_from_px(0.002, x)
            Hl_max = Hl_from_px(2, x)
            if Hl < Hl_min or Hl > Hl_max:
                raise AmmoniaFunctionError('Para una fracción de amoniaco en la fase líquida igual a ' + str(liquid_phase_q) + ', la entalpía de la fase líquida debe ser mayor o igual a ' + str(1000*Hl_min) + ' J/kg y menor o igual a ' + str(1000*Hl_max) + '  J/kg.')
            p = p_from_xHl(x, Hl, [0.002,2])
            if variable_required == 'p':
                return 1e6*p
            if variable_required == 'T':
                return T_from_px(p, x)
        if variable_required == 'Qg' or variable_required == 'Hg':
            if x < 0.05 or x >= 1:
                raise AmmoniaFunctionError('La fracción de amoniaco en la fase líquida debe ser mayor o igual a ' + str(molar_to_massic_quality(0.05)) + ' y menor a uno.')
            Hl_min = Hl_from_px(0.02, x)
            Hl_max = Hl_from_px(2, x)
            if Hl < Hl_min or Hl > Hl_max:
                raise AmmoniaFunctionError('Para una fracción de amoniaco en la fase líquida igual a ' + str(liquid_phase_q) + ', la entalpía de la fase líquida debe ser mayor o igual a ' + str(1000*Hl_min) + ' J/kg y menor o igual a ' + str(1000*Hl_max) + '  J/kg.')
            p = p_from_xHl(x, Hl, [0.02,2])
            y = y_from_px(p, x)
            if variable_required == 'Qg':
                return molar_to_massic_quality(y)
            if variable_required == 'Hg':
                return 1000*Hg_from_py(p, y)
        
    ## Caso 8: Hg, Qg (FALTA RESTRICCIÓN PARA QUE x > 0.05)
    if (variable_1_name == 'Hg' and variable_2_name == 'Qg') or (variable_1_name == 'Qg' and variable_2_name == 'Hg'):
        if variable_1_name == 'Hg':
            Hg = variable_1_value/1000
            gas_phase_q = variable_2_value
        else:
            gas_phase_q = variable_1_value
            Hg = variable_2_value/1000
        y = massic_to_molar_quality(gas_phase_q)
        if y <= 0 or y >= 1:
            raise AmmoniaFunctionError('La fracción de amoniaco en la fase gaseosa debe ser mayor a cero y menor a uno.')
        Hg_min = Hg_from_py(0.02, y)
        Hg_max = Hg_from_py(2, y)
        if Hg < Hg_min or Hg > Hg_max:
            raise AmmoniaFunctionError('Para una fracción de amoniaco en la fase gaseosa igual a ' + str(gas_phase_q) + ', la entalpía de la fase gaseosa debe ser mayor o igual a ' + str(1000*Hg_min) + ' J/kg y menor o igual a ' + str(1000*Hg_max) + '  J/kg.')
        p = p_from_yHg(y, Hg, [0.02, 2])
        if variable_required == 'p':
            return 1e6*p
        if variable_required == 'T':
            return T_from_py(p, y)
        x = x_from_py(p, y)
        if variable_required == 'Ql':
            return molar_to_massic_quality(x)
        if variable_required == 'Hl':
            return 1000*Hl_from_px(p, x)
    
    ## Caso 9: T, Hl. Prohibido, función no monótona
    
    ## Caso 10: T, Hg
    if (variable_1_name == 'T' and variable_2_name == 'Hg') or (variable_1_name == 'Hg' and variable_2_name == 'T'):
        if variable_1_name == 'T':
            T = variable_1_value
            Hg = variable_2_value/1000
        else:
            Hg = variable_1_value/1000
            T = variable_2_value
        T_max_0_02 = T_from_py(0.02, 0)
        T_min_0_02 = T_from_py(0.02, 1)
        T_max_2 = T_from_py(2, 0)
        T_min_2 = T_from_py(2, 1)
        if T >= T_max_2:
            raise AmmoniaFunctionError('La temperatura debe ser menor a ' + str(T_max_2) + ' K.')
        elif T <= T_min_0_02:
            raise AmmoniaFunctionError('La temperatura debe ser mayor a ' + str(T_min_0_02) + ' K.')
        elif T > T_max_0_02:
            y_min = 0
            y_max = y_from_pT(2, T)
        elif T < T_min_2:
            y_min = y_from_pT(0.02, T)
            y_max = 1
        else:
            y_min = y_from_pT(0.02, T)
            y_max = y_from_pT(2, T)
        Hg_max = Hg_from_Ty(T, y_min)
        Hg_min = Hg_from_Ty(T, y_max)
        if y_min == 0:
            condition_1 = Hg < Hg_max
            text1 = ' menor a '
        else:
            condition_1 = Hg <= Hg_max
            text1 = ' menor o igual a '
        if y_max == 1:
            condition_2 = Hg > Hg_min
            text2 = ' mayor a '
        else:
            condition_2 = Hg >= Hg_min
            text2 = ' mayor o igual a '
        if not (condition_1 and condition_2):
            raise AmmoniaFunctionError('Para una temperatura de ' + str(T) + ' K, la entalpía de la fase gaseosa debe ser' + text2 + str(1000*Hg_min) + ' K y' + text1 + str(1000*Hg_max) + ' K.')
        y = y_from_HgT(Hg, T)
        if variable_required == 'Qg':
            return molar_to_massic_quality(y)
        p = p_from_Ty(T, y, [0.01,2.5])
        if variable_required == 'p':
            return 1e6*p
        x = x_from_py(p, y)
        if variable_required == 'Ql':
            return molar_to_massic_quality(x)
        if variable_required == 'Hl':
            return 1000*Hl_from_Tx(T, x)
        
    ## Caso 11: p, Hl. Prohibido, función no monótona

    ## Caso 12: p, Hg
    if (variable_1_name == 'p' and variable_2_name == 'Hg') or (variable_1_name == 'Hg' and variable_2_name == 'p'):
        if variable_1_name == 'p':
            p_orig = variable_1_value
            Hg = variable_2_value/1000
        else:
            Hg = variable_1_value/1000
            p_orig = variable_2_value
        if p_orig < 2000 or p_orig > 2e6:
            raise AmmoniaFunctionError('La presión debe ser mayor o igual a 2 kPa y menor o igual a 2 MPa (especificar presión en Pa).')
        p = p_orig/1e6
        Hg_max = Hg_from_py(p, 0)
        Hg_min = Hg_from_py(p, 1)
        if Hg <= Hg_min or Hg >= Hg_max:
            raise AmmoniaFunctionError('Para una presión de ' + str(p_orig) + ' Pa, la entalpía de la fase gaseosa debe ser mayor a ' + str(1000*Hg_min) + ' J/kg y menor a ' + str(1000*Hg_max) + ' J/kg.')
        y = y_from_pHg(p, Hg)
        if variable_required == 'Qg':
            return molar_to_massic_quality(y)
        if variable_required == 'T':
            return T_from_py(p, y)
        x = x_from_py(p, y)
        if variable_required == 'Ql':
            return molar_to_massic_quality(x)
        if variable_required == 'Hl':
            return 1000*Hl_from_px(p, x)
        
    ## Caso 13: Qg, Hl
    if (variable_1_name == 'Qg' and variable_2_name == 'Hl') or (variable_1_name == 'Hl' and variable_2_name == 'Qg'):
        if variable_1_name == 'Qg':
            gas_phase_q = variable_1_value
            Hl = variable_2_value/1000
        else:
            Hl = variable_1_value/1000
            gas_phase_q = variable_2_value
        y = massic_to_molar_quality(gas_phase_q)
        if y <= 0 or y >= 96:
            raise AmmoniaFunctionError('La fracción de amoniaco en la fase gaseosa debe ser mayor a cero y menor o igual a ' + str(molar_to_massic_quality(0.96)) + '.')
        Hl_min = Hl_from_py(0.02, y)
        Hl_max = Hl_from_py(2, y)
        if Hl < Hl_min or Hl > Hl_max:
            raise AmmoniaFunctionError('Para una fracción de amoniaco en la fase gaseosa igual a ' + str(gas_phase_q) + ', la entalpía de la fase líquida debe ser mayor o igual a ' + str(Hl_min) + ' J/kg y menor o igual a ' + str(Hl_max) + ' J/kg.')
        p = p_from_yHl(y, Hl)
        if variable_required == 'p':
            return 1e6*p
        if variable_required == 'T':
            return T_from_py(p, y)
        if variable_required == 'Hg':
            return 1000*Hg_from_py(p, y)
        x = x_from_py(p, y)
        if variable_required == 'Ql':
            return molar_to_massic_quality(x)
        
    ## Caso 14: Ql, Hg
    if (variable_1_name == 'Ql' and variable_2_name == 'Hg') or (variable_1_name == 'Hg' and variable_2_name == 'Ql'):
        if variable_1_name == 'Ql':
            liquid_phase_q = variable_1_value
            Hg = variable_2_value/1000
        else:
            Hg = variable_1_value/1000
            liquid_phase_q = variable_2_value
        x = massic_to_molar_quality(liquid_phase_q)
        if x <= 0 or x > 0.5:
            raise AmmoniaFunctionError("Para usar 'Ql' y 'Hg' como variables de entrada, la fracción de amoniaco en la fase líquida debe ser mayor a cero y menor o igual a " + str(molar_to_massic_quality(0.5)) + ".")
        Hg_min = Hg_from_px(0.02, x)
        Hg_max = Hg_from_px(2, x)
        if Hg < Hg_min or Hg > Hg_max:
            raise AmmoniaFunctionError('Para una fracción de amoniaco en la fase líquida igual a ' + str(liquid_phase_q) + ', la entalpía de la fase gaseosa debe ser mayor o igual a ' + str(1000*Hg_min) + ' J/kg y menor o igual a ' + str(1000*Hg_max) + ' kJ/kg.')
        p = p_from_xHg(x, Hg)
        if variable_required == 'p':
            return 1e6*p
        if variable_required == 'T':
            return T_from_px(p, x)
        if variable_required == 'Hl':
            return 1000*Hl_from_px(p, x)
        if variable_required == 'Qg':
            y = y_from_px(p, x)
            return molar_to_massic_quality(y)