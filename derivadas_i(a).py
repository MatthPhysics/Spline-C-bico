# -*- coding: utf-8 -*-
"""derivadas_i(a).ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1uz-ndKvfrc3A0cwFOQoeMlDDK4S4XBHX
"""

from sympy import *
import math

init_printing(pretty_print=True)

x = Symbol('x')

cos(pi * x).diff(x, 1)

cos(pi * x).diff(x, 1).subs({x: 0.5})

cos(pi * x).diff(x, 2).subs({x: 0.5})

D_0 = (1.0 - (0.7573593128807148 * (x - 0.0)) + (0.0 * (x - 0.0)**2) - (6.62741699796952*(x-0.0)**3)).diff(x,1).subs({x:0.5})
D_0_2 = (1.0 - (0.7573593128807148 * (x - 0.0)) + (0.0 * (x - 0.0)**2) - (6.62741699796952*(x-0.0)**3)).diff(x,2).subs({x:0.5})
print("Resultado da derivada primeira: ", D_0)
print("Resultado da derivada segunda: ", D_0_2)

D_1 = (0.7071067811865476 - (1.9999999999999998 * (x-0.25)) - (4.97056274847714*(x -0.25)**2)+ (6.627416997969518*(x-0.25)**3)).diff(x,1).subs({x:0.5})
D_1_2 = (0.7071067811865476 - (1.9999999999999998 * (x-0.25)) - (4.97056274847714*(x -0.25)**2)+ (6.627416997969518*(x-0.25)**3)).diff(x,2).subs({x:0.5})
print("Resultado da derivada primeira: ", D_1)
print("Resultado da derivada segunda: ", D_1_2)

D_2 = (6.123233995736766e-17 - (3.242640687119285 * (x-0.5)) - (1.8947806286936005e-15 * (x -0.5)**2)+ (6.627416997969523*(x-0.5)**3)).diff(x,1).subs({x:0.5})
D_2_2 = (6.123233995736766e-17 - (3.242640687119285 * (x-0.5)) - (1.8947806286936005e-15 * (x -0.5)**2)+ (6.627416997969523*(x-0.5)**3)).diff(x,2).subs({x:0.5})
print("Resultado da derivada primeira: ", D_2)
print("Resultado da derivada segunda: ", D_2_2)

D_3 = (-0.7071067811865475 - (2.0 * (x-0.75)) + (4.970562748477141*(x -0.75)**2) - (6.627416997969521*(x-0.75)**3)).diff(x,1).subs({x:0.5})
D_3_2 = (-0.7071067811865475 - (2.0 * (x-0.75)) + (4.970562748477141*(x -0.75)**2) - (6.627416997969521*(x-0.75)**3)).diff(x,2).subs({x:0.5})
print("Resultado da derivada primeira: ", D_3)
print("Resultado da derivada segunda: ", D_3_2)

DF_0 = (1.0 + (0.0 * (x-0.0)) - (5.193321002610615*(x -0.0)**2)+(2.028118006381504*(x-0.0)**3)).diff(x,1).subs({x:0.5})
DF_0_2 = (1.0 + (0.0 * (x-0.0)) - (5.193321002610615*(x -0.0)**2)+(2.028118006381504*(x-0.0)**3)).diff(x,2).subs({x:0.5})
print("Resultado da derivada primeira: ", DF_0)
print("Resultado da derivada segunda: ", DF_0_2)

DF_1 = (0.7071067811865476-(2.2163883751087754 * (x-0.25)) - (3.672232497824487*(x -0.25)**2)+ (4.896309997099313*(x-0.25)**3)).diff(x,1).subs({x:0.5})
DF_1_2 = (0.7071067811865476-(2.2163883751087754 * (x-0.25)) - (3.672232497824487*(x -0.25)**2)+ (4.896309997099313*(x-0.25)**3)).diff(x,2).subs({x:0.5})
print("Resultado da derivada primeira: ", DF_1)
print("Resultado da derivada segunda: ", DF_1_2)

DF_2 = (6.123233995736766e-17 - (3.134446499564897 * (x-0.5)) - (2.0325621527752865e-15*(x -0.5)**2) + (4.89630999709932*(x-0.5)**3)).diff(x,1).subs({x:0.5})
DF_2_2 = (6.123233995736766e-17 - (3.134446499564897 * (x-0.5)) - (2.0325621527752865e-15*(x -0.5)**2) + (4.89630999709932*(x-0.5)**3)).diff(x,2).subs({x:0.5})
print("Resultado da derivada primeira: ", DF_2)
print("Resultado da derivada segunda: ", DF_2_2)

DF_3 = (-0.7071067811865475 - (2.216388375108776 * (x-0.75))+ (3.6722324978244876*(x -0.75)**2) + (2.028118006381503*(x-0.75)**3)).diff(x,1).subs({x:0.5})
DF_3_2 = (-0.7071067811865475 - (2.216388375108776 * (x-0.75))+ (3.6722324978244876*(x -0.75)**2) + (2.028118006381503*(x-0.75)**3)).diff(x,2).subs({x:0.5})
print("Resultado da derivada primeira: ", DF_3)
print("Resultado da derivada segunda: ", DF_3_2)