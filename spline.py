import numpy as np
import math
import matplotlib.pyplot as plt

def spline(x, y):
	
	
	n = len(x)
	a = {k: v for k, v in enumerate(y)}
	h = {k: x[k+1] - x[k] for k in range(n - 1)}

	A = [[1] + [0] * (n-1)]
	for i in range(1, n - 1):
		linha = [0] * n
		linha[i-1] = h[i-1]
		linha[i] = 2*(h[i-1] + h[i])
		linha[i + 1] = h[i]	
		A.append(linha)

	A.append([0] * (n - 1) + [1])


	B = [0]

	for k in range(1, n-1):
		valor = 3*(a[k+1] -a[k]) / h[k] - 3*(a[k] - a[k-1]) / h[k-1]
		B.append(valor)
	B.append(0)
	print("Os valores de B: ", B)
	c = dict(zip(range(n), np.linalg.solve(A, B)))

	b = {}
	d = {}
	for k in range(n - 1):
		b[k] = (1/h[k]) * (a[k + 1] - a[k]) - (h[k]/3) * (2*c[k]+c[k+1])
		d[k] = (c[k+1] - c[k])/ (3 * h[k])	


	s = {}
	#print("Os coeficientes calculados são: \n")
	#for k in range(n - 1):
	#	print(a[k],"\n",b[k],"\n",c[k],"\n",d[k])	
	
	for k in range(n - 1):
		eq = f'{a[k]}{b[k]:+} * (x{-x[k]:+}){c[k]:+}*(x {- x[k]:+})**2{d[k]:+}*(x{-x[k]:+})**3'
		s[k] = {'eqs': eq, 'dominio': [x[k], x[k+1]]}
	
	print("Os splines para este conjunto de pontos são: \n")
	for k, v in s.items():
		print(f'S{k}', v)

	
	return s

x_1 = [0., 0.25, 0.5, 0.75, 1.]
x_2 = [0., 0.125, 0.250, 0.375, 0.5, 0.625, 0.75, 0.875, 1.]
y_1 = [math.cos(math.pi*i) for i in x_1]
y_2 = [math.cos(math.pi*i) for i in x_2]


eqs = spline(x_2, y_2)
print("\n")

#gera os gráficos dos splines 
for key , value in eqs.items():
	def p(x):
		return eval(value['eqs'])
	t = np.linspace(*value['dominio'], 100)
	plt.plot(t, p(t), label=f"$S_{key}(x)$")

plt.title('Interpolação: Spline natural (9 pontos)')
plt.scatter(x_2, y_2)
plt.legend()
plt.show()






 
