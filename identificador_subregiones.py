""" Este código busca identificar las regiones del dominio con diferentes índices de refracción, de tal forma que
crea dos vectores (llamados ref_indi_puntos y ref_indi_elem) en cuyas posiciones se encuentran los índices asociados a
 los nodos y a los elementos. Además, grafica los diferentes subdominios para confirmarle al usario las definiciones
 hechas


 El usuario debe ingresar el path del archivo .msh y debe ingresar en el array 'indices' los diferentes índicdes de
 refracción para las diferentes regiones que haya definidio durante la creación de la malla. Por ejemplo, para un
 dominio con 2 regiones, se ingresan dos índices de refracción. El orden de ingreso de estos índices debe verificarse
 con la gráfica que entrega el programa


 Finalmente, el programa crea un archivo .msh con los dos vectores previamente explicados añadidos como point_data y
 cell_data, este archivo lo guarda en otra ubicación que debe ingresar el usuario"""

import meshio
import numpy as np
import matplotlib.pyplot as plt


path_entrada = 'mallas_msh/circulo1.msh'
refs = np.array([1, 2])


mesh = meshio.read(path_entrada)
puntos = mesh.points
elementos = mesh.cells['triangle']

regiones = mesh.cell_data['triangle']['gmsh:physical']
regiones = np.array(regiones)

npuntos = np.shape(puntos)[0]
nelementos = np.shape(elementos)[0]

subdoms_key = np.unique(regiones)

ref_puntos = np.zeros(npuntos)

for i in range(nelementos):
    for j in range(np.size(subdoms_key)):
        if regiones[i] == subdoms_key[j]:

            p0 = elementos[i, 0]
            p1 = elementos[i, 1]
            p2 = elementos[i, 2]
            ref_puntos[p0]= j+1
            ref_puntos[p1] = j+1
            ref_puntos[p2] = j+1

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(puntos[:, 0], puntos[:, 1], ref_puntos, zdir='z', s=20, c=None, depthshade=True)
