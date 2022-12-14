# -*- coding: utf-8 -*-
"""derivadas_i(b).ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1rT269kqRYzSDCPToczjnVMkKNxcjO4-c
"""

from sympy import *
import math

init_printing(pretty_print=True)

x = Symbol('x')

cos(pi * x).diff(x, 1)

cos(pi * x).diff(x, 1).subs({x: 0.5})

cos(pi * x).diff(x, 2).subs({x: 0.5})

D_0 = (1.0 - (0.36075774689582807 * (x-0.0)) + (0.0*(x -0.0)**2)- (15.885183552888194*(x-0.0)**3)).diff(x, 1).subs({x: 0.5})
D_0_2 = (1.0 - (0.36075774689582807 * (x-0.0)) + (0.0*(x -0.0)**2)- (15.885183552888194*(x-0.0)**3)).diff(x, 2).subs({x: 0.5})
print("Resultado da derivada primeira: ", D_0)
print("Resultado da derivada segunda: ", D_0_2)

D_1 = (0.9238795325112867-(1.1053757259374621 * (x-0.125)) - (5.956943832333073*(x -0.125)**2) + (7.411948440395709*(x-0.125)**3)).diff(x, 1).subs({x: 0.5})
D_1_2 = (0.9238795325112867-(1.1053757259374621 * (x-0.125)) - (5.956943832333073*(x -0.125)**2) + (7.411948440395709*(x-0.125)**3)).diff(x, 2).subs({x: 0.5}) 
print("Resultado da derivada primeira: ", D_1)
print("Resultado da derivada segunda: ", D_1_2)

D_2 = (0.7071067811865476- (2.2471766008771814 * (x-0.25)) - (3.1774631671846816*(x -0.25)**2) + (3.134253197030707*(x-0.25)**3)).diff(x, 1).subs({x: 0.5})
D_2_2 = (0.7071067811865476- (2.2471766008771814 * (x-0.25)) - (3.1774631671846816*(x -0.25)**2) + (3.134253197030707*(x-0.25)**3)).diff(x, 2).subs({x: 0.5})
print("Resultado da derivada primeira: ", D_2)
print("Resultado da derivada segunda: ", D_2_2)

D_3 = (0.38268343236508984-(2.894624274062538 * (x-0.375))- (2.0021182182981665*(x -0.375)**2) + (5.338981915461776*(x-0.375)**3)).diff(x, 1).subs({x: 0.5})
D_3_2 = (0.38268343236508984-(2.894624274062538 * (x-0.375))- (2.0021182182981665*(x -0.375)**2) + (5.338981915461776*(x-0.375)**3)).diff(x, 2).subs({x: 0.5})
print("Resultado da derivada primeira: ", D_3)
print("Resultado da derivada segunda: ", D_3_2)

D_4 = (6.123233995736766e-17 - (3.1448890513498085 * (x-0.5)) - (5.94952051473768e-16*(x -0.5)**2) + (5.338981915461786*(x-0.5)**3)).diff(x, 1).subs({x: 0.5})
D_4_2 = (6.123233995736766e-17 - (3.1448890513498085 * (x-0.5)) - (5.94952051473768e-16*(x -0.5)**2) + (5.338981915461786*(x-0.5)**3)).diff(x, 2).subs({x: 0.5})
print("Resultado da derivada primeira: ", D_4)
print("Resultado da derivada segunda: ", D_4_2)

D_5 = (-0.3826834323650897- (2.8946242740625374 * (x-0.625)) + (2.002118218298169*(x -0.625)**2) + (3.1342531970306786*(x-0.625)**3)).diff(x, 1).subs({x: 0.5})
D_5_2 = (-0.3826834323650897- (2.8946242740625374 * (x-0.625)) + (2.002118218298169*(x -0.625)**2) + (3.1342531970306786*(x-0.625)**3)).diff(x, 2).subs({x: 0.5})
print("Resultado da derivada primeira: ", D_5)
print("Resultado da derivada segunda: ", D_5_2)

D_6 = (-0.7071067811865475-(2.247176600877182 * (x-0.75)) + (3.1774631671846736*(x -0.75)**2) + (7.411948440395752*(x-0.75)**3)).diff(x, 1).subs({x: 0.5})
D_6_2 = (-0.7071067811865475-(2.247176600877182 * (x-0.75)) + (3.1774631671846736*(x -0.75)**2) + (7.411948440395752*(x-0.75)**3)).diff(x, 2).subs({x: 0.5})
print("Resultado da derivada primeira: ", D_6)
print("Resultado da derivada segunda: ", D_6_2)

D_7 = (-0.9238795325112867- (1.1053757259374628 * (x-0.875)) + (5.956943832333081*(x -0.875)**2) - (15.885183552888215*(x-0.875)**3)).diff(x, 1).subs({x: 0.5})
D_7_2 = (-0.9238795325112867- (1.1053757259374628 * (x-0.875)) + (5.956943832333081*(x -0.875)**2) - (15.885183552888215*(x-0.875)**3)).diff(x, 2).subs({x: 0.5})
print("Resultado da derivada primeira: ", D_7)
print("Resultado da derivada segunda: ", D_7_2)

DF_0 = (1.0+ (0.0 * (x-0.0)) - (4.998540328123639*(x -0.0)**2) + (1.0146432707679172*(x-0.0)**3)).diff(x, 1).subs({x: 0.5})
DF_0_2 = (1.0+ (0.0 * (x-0.0)) - (4.998540328123639*(x -0.0)**2) + (1.0146432707679172*(x-0.0)**3)).diff(x, 2).subs({x: 0.5})
print("Resultado da derivada primeira: ", DF_0)
print("Resultado da derivada segunda: ", DF_0_2)

DF_1 = (0.9238795325112867- (1.2020736787136634 * (x-0.125)) - (4.61804910158567*(x -0.125)**2) + (2.889459572093371*(x-0.125)**3)).diff(x, 1).subs({x: 0.5})
DF_1_2 = (0.9238795325112867- (1.2020736787136634 * (x-0.125)) - (4.61804910158567*(x -0.125)**2) + (2.889459572093371*(x-0.125)**3)).diff(x, 2).subs({x: 0.5}) 
print("Resultado da derivada primeira: ", DF_1)
print("Resultado da derivada segunda: ", DF_1_2)

DF_2 = (0.7071067811865476 - (2.221142536668204 * (x-0.25)) - (3.5345017620506556*(x -0.25)**2) + (4.324381846583953*(x-0.25)**3)).diff(x, 1).subs({x: 0.5})
DF_2_2 = (0.7071067811865476 - (2.221142536668204 * (x-0.25)) - (3.5345017620506556*(x -0.25)**2) + (4.324381846583953*(x-0.25)**3)).diff(x, 2).subs({x: 0.5})
print("Resultado da derivada primeira: ", DF_2)
print("Resultado da derivada segunda: ", DF_2_2)

DF_3 = (0.38268343236508984 - (2.9020625781222456 * (x-0.375)) - (1.912858569581673*(x -0.375)**2) + (5.1009561855511265*(x-0.375)**3)).diff(x, 1).subs({x: 0.5})
DF_3_2 = (0.38268343236508984 - (2.9020625781222456 * (x-0.375)) - (1.912858569581673*(x -0.375)**2) + (5.1009561855511265*(x-0.375)**3)).diff(x, 2).subs({x: 0.5})
print("Resultado da derivada primeira: ", DF_3)
print("Resultado da derivada segunda: ", DF_3_2)

DF_4 = (-0.3826834323650897- (2.902062578122245 * (x-0.625)) + (1.9128585695816755*(x -0.625)**2) + (4.3243818465839245*(x-0.625)**3)).diff(x, 1).subs({x: 0.5})
DF_4_2 = (-0.3826834323650897- (2.902062578122245 * (x-0.625)) + (1.9128585695816755*(x -0.625)**2) + (4.3243818465839245*(x-0.625)**3)).diff(x, 2).subs({x: 0.5})
print("Resultado da derivada primeira: ", DF_4)
print("Resultado da derivada segunda: ", DF_4_2)

DF_5 = (-0.3826834323650897- (2.902062578122245 * (x-0.625)) + (1.9128585695816755*(x -0.625)**2) + (4.3243818465839245*(x-0.625)**3)).diff(x, 1).subs({x: 0.5})
DF_5_2 = (-0.3826834323650897- (2.902062578122245 * (x-0.625)) + (1.9128585695816755*(x -0.625)**2) + (4.3243818465839245*(x-0.625)**3)).diff(x, 2).subs({x: 0.5})
print("Resultado da derivada primeira: ", DF_5)
print("Resultado da derivada segunda: ", DF_5_2)

DF_6 = (-0.7071067811865475-(2.221142536668205 * (x-0.75)) + (3.534501762050647*(x -0.75)**2) + (2.8894595720934197*(x-0.75)**3)).diff(x, 1).subs({x: 0.5})
DF_6_2 = (-0.7071067811865475-(2.221142536668205 * (x-0.75)) + (3.534501762050647*(x -0.75)**2) + (2.8894595720934197*(x-0.75)**3)).diff(x, 2).subs({x: 0.5})
print("Resultado da derivada primeira: ", DF_6)
print("Resultado da derivada segunda: ", DF_6_2)

DF_7 = (-0.9238795325112867- (1.2020736787136639 * (x-0.875)) + (4.6180491015856795*(x -0.875)**2) + (1.0146432707678652*(x-0.875)**3)).diff(x, 1).subs({x: 0.5})
DF_7_2 = (-0.9238795325112867- (1.2020736787136639 * (x-0.875)) + (4.6180491015856795*(x -0.875)**2) + (1.0146432707678652*(x-0.875)**3)).diff(x, 2).subs({x: 0.5})
print("Resultado da derivada primeira: ", DF_7)
print("Resultado da derivada segunda: ", DF_7_2)