import streamlit as st
#import pandas as pd
import numpy as np
from micropile_buckling.src.micropile_buckling import (get_cross_section_parameters_circular, get_w_f_elastoplastic_eq1, get_p_f_elastoplastic_eq2, 
                                    get_p_f_elastoplastic_eq3, get_p_f_elastoplastic_eq4, get_Ncr_by_iteration, get_Nb_Rd_EC3,
                                    display_micropile_buckling_resistance)

st.title('Micro-pile buckling check')
st.subheader('Vogt, N., & Vogt, S. (2013). Biegeknickwiderstand von Mikropfählen gemäß den Eurocodes. Bautechnik, 90(9), 550-558.')

# Cross-sectional parameters
st.header('Cross-sectional parameters')
col1, col2 = st.columns(2)
d = float(col1.text_input('Steel bar diameter [mm]', 50.0))
D = float(col2.text_input('Pile shaft diameter [mm]', 270.0)) 
E = float(col1.text_input("Young's modulus [kN/m^2]", 200000.0))
fy = float(col2.text_input('Characteristic yield stress [kN/m^2]', 500.0))
A, I, EI, fyA = get_cross_section_parameters_circular(d, D, E, fy)
st.write('A = {0:.2f} [mm^2]; I = {1:.2f} [cm^4]; EI = {2:.2f} [kNm^2]; fyA = {3:.2f} [kN]\n\n'.format(A, I, EI, fyA))

# Selection for p-y cur2e
st.header('Maximal mobilized lateral displacement and stress (soil-pile interaction)')
col1, col2, col3 = st.columns(3)
c_u = float(col1.text_input('Undrained shear strength [kN/m^2]', 25.0))
w_f = get_w_f_elastoplastic_eq1(D, c_u)

eq_p_f = col2.selectbox('Pile contact surface (with soil)', ('Ideally smooth', 'Ideally rough', 'Pile not completely encircled by soil', 'User defined'), index=1)
if eq_p_f == 'Ideally smooth':
    p_f = get_p_f_elastoplastic_eq2(c_u)
elif eq_p_f == 'Ideally rough':
    p_f = get_p_f_elastoplastic_eq3(c_u)
elif eq_p_f == 'Pile not completely encircled by soil':
    p_f = get_p_f_elastoplastic_eq4(c_u)
else: # user defined
    a = float(col3.text_input('User-defined constant a, for p_f = a*c_u', 10.5))
    p_f = a*c_u

col1.write('w_f = {0:.2f} [mm]'.format(w_f))
col2.write('p_f = {0:.2f} [kN/m^2]\n\n'.format(p_f))

# Micropile length and buckling parameters
st.header('Micropile length and buckling parameters')
col1, col2 = st.columns(2)
L = float(col1.text_input('Pile length [m]', 5.0))
k_imp = float(col2.text_input('Imperfection parameter kappa [m]', 1/200))        # [1/m]
buckling_curve = col1.selectbox('Buckling curve (Figure 6.4, EC3)', ('a0', 'a', 'b', 'c', 'd'), index=3) 
gamma_M = float(col2.text_input('gamma_M', 1.1))

# Buckling resistance check for the given c_u value
st.header('Design resistance against buckling')
col1, col2, col3 = st.columns(3)
Lcr, Ncr = get_Ncr_by_iteration(w_f*0.001, EI, p_f, D*0.001, k_imp, L)
Nb_Rd = get_Nb_Rd_EC3(fyA, Ncr, buckling_curve, gamma_M)
col1.write('Lcr = {0:.2f} [m]'.format(Lcr))
col2.write('Ncr = {0:.2f} [kN]'.format(Ncr))
col3.write('Nb_Rd = {0:.2f} [kN]'.format(Nb_Rd))

# Data generation and plotting
st.header('Plot of bifurcation load vs. undrained shear strength')
fig = display_micropile_buckling_resistance(D, EI, fyA, k_imp, L, p_f, c_u, Nb_Rd, buckling_curve, gamma_M)
st.pyplot(fig)