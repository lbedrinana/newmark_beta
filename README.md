# newmark_beta
#=========================================================          
#     newmark.py
#=========================================================
# Coded by Luis Bedriñana,
#          Universidad de Ingenieria y Tecnologia - UTEC    
#          Jun, 2020
# Version: 1.1
# 
# Script to analyze the dynamic linear response of 1-DOF subjected to GM record by the
# Newmark beta method
#
# INPUT:
#   per: period of the analyzed system
#   damp: Equivalent damping ratio of the analyzed system
#   GM_acc: Ground motion record read form a text file 
#
# OUTPUT:
#   Out_acce: Time history of the acceleration response
#   Out_disp: Time history of the displacement response 
#   peak_acc: Peak acceleration
#   peak_disp: peak displacement
