""" 

Código que ordena las soluciones encontradas por los solucionadores de Python
y COMSOL, la idea es darle un orden general a los nodos, primero de menor a 
mayor en x, y posteriormente de menor a mayor en y 

"""

import numpy as np

# 'simulaciones_comsol/datos_slnes/comsol9.txt'
# 'datos_python/python2.txt'
filename_in = 'datos_python/python14.txt'
filename_out = filename_in
datos = np.loadtxt(filename_in)
datosb = np.copy(datos)
idx_sort = np.argsort(datos[:,0])

datos = datos[idx_sort]
datosc = np.copy(datos)

for i in range(np.shape(datos)[0] - 1):
    if abs(datos[i,0]- datos[i+1,0]) < 0.00005:
        datos[i+1,0] = datos[i,0]


vals, idx_start, cont = np.unique(datos[:,0], return_index = True, return_counts= True)

for i in range(np.size(cont)):
    if cont[i] != 1:
        inicio = idx_start[i]
        fin = idx_start[i] + cont[i]
        b = datos[inicio:fin, :]
        idx_sort = np.argsort(b[:,1])
        b = b[idx_sort]
        datos[inicio:fin, :] = b 

np.savetxt(filename_out, datos)
