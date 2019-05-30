""" Modulo con funciones de apoyo para el programa FEM2D_main.py

Contiene las funciones:
    + leer_msh_imponer_materiales
    
    + hallar_det_jaco_tria
    
    + hallar_jaco_tria
    
    + hallar_matriz_rigidez_local
    
    + hallar_matriz_masa_local
    
    + hallar_vector_cargas_local
    
    + ensamble

"""

import numpy as np
import meshio

D = np.array([[-1, 1, 0],[-1, 0 ,1]])
Dt = np.transpose(D)

#%%

def leer_msh_imponer_materiales(filename, mr1, er1, mr2, er2):
    """ funcion que lee el archivo .vtk y lo separa el 4 arrays de numpy, 
    uno con los puntos, otro con elementos y otros 2  con los parametros del material
    asociados a los puntos"""
    
    mesh = meshio.read(filename)  
    puntos = np.array(mesh.points) 
    npuntos = np.shape(puntos)[0]
    cells = mesh.cells
    elementos = np.array(cells['triangle'])
    nelementos = np.shape(elementos)[0]

    # en esta parte se definen los parámetros de las regiones a partir de las relativas
    permit_idx = np.array([er1, er2])
    permea_idx = np.array([mr1, mr2])
    
    regiones = mesh.cell_data['triangle']['gmsh:physical']
    regiones = np.array(regiones)
    
    subdoms_key = np.unique(regiones)

    permit_elem = np.zeros(nelementos)
    permea_elem = np.zeros(nelementos)
    
    for i in range(nelementos):
        for j in range(np.size(subdoms_key)):
            if regiones[i] == subdoms_key[j]:
                permit_elem[i] = permit_idx[j]
                permea_elem[i] = permea_idx[j]

    return puntos, elementos, permit_elem, permea_elem, npuntos, nelementos

#%%
    
def hallar_det_jaco_tria(elemento, puntos): # revisada
    """ funcion que devuelve el determinante del jacobiano para un 
    elemento triangular de 3 nodos"""
    
    p0 = puntos[elemento[0]]
    p1 = puntos[elemento[1]]
    p2 = puntos[elemento[2]]
    det_jaco = (p1[0] - p0[0])*(p2[1] - p0[1]) - (p2[0] - p0[0])*(p1[1] - p0[1])
    return det_jaco

#%%
    
def hallar_jaco_tria(elemento, puntos):
    """ funcion que devuelve el jacobiano para un 
    elemento triangular de 3 nodos"""
    
    p0 = puntos[elemento[0]]
    p1 = puntos[elemento[1]]
    p2 = puntos[elemento[2]]
    jaco = np.array([[p1[0] - p0[0], p1[1] - p0[1] ],[p2[0] - p0[0], p2[1] - p0[1] ]])
    return jaco

#%%
    
def hallar_matriz_rigidez_local(elemento, puntos):
    """ funcion que devuelve la matriz de rigidez local para un elemento"""
    
    invJ = np.linalg.inv(hallar_jaco_tria(elemento,puntos))
    invJ_t = np.transpose(invJ)
    a_local = (hallar_det_jaco_tria(elemento, puntos)/2) * (Dt @ invJ_t @ invJ @ D)
    return a_local

#%%
    
def hallar_matriz_masa_local(elemento, puntos):
    """ funcion que devuelve la matriz de masa local para un elemento"""
    
    m_local = (hallar_det_jaco_tria(elemento, puntos)/24) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]])
    return m_local

#%%
    
def hallar_vector_cargas_local(elemento, puntos, fuente):
    """ funcion que devuelve vector de cargas local para un elemento"""
    
    p0 = puntos[elemento[0]]
    p1 = puntos[elemento[1]]
    p2 = puntos[elemento[2]]
    det = hallar_det_jaco_tria(elemento, puntos)
    centroide_x = (p0[0] + p1[0] +p2[0])/(3)
    centroide_y = (p0[1] + p1[1] + p2[1]) /( 3)
    b_local = (fuente(centroide_x, centroide_y) * (1/6) * det) * np.ones(3)
    return b_local

#%%
    
def ensamble(elemento, a_local, m_local, b_local, a_global, m_global, b_global):
    """ funcion que ensambla las matrices locales en las posiciones pertinentes
    de las matrices globales"""
    
    for i in range(3):
        for j in range(3):
            a_global[elemento[i], elemento[j]] = a_local[i, j] + a_global[elemento[i], elemento[j]]
            m_global[elemento[i], elemento[j]] = m_local[i, j] + m_global[elemento[i], elemento[j]]
        b_global[elemento[i]] = b_local[i] + b_global[elemento[i]]
    return a_global, m_global, b_global

