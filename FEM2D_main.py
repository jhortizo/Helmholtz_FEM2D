"""
Archivo definitivo del proyecto, soluciona la ecuacion de Helmholtz para 2 dimensiones
Desarrollado por Sofia Oando Vasquez, Jose Hernan Ortiz y Jose Manuel Rendon
para la materia Metodos numericos, a cargo del profesor Nicolas Guarin-Zapata


El usuario debe ingresar:

    + Archivo con malla, habiendo identificado las regiones con identificador_subregiones.py
    + Función fuente
    + Frecuencia angular de la onda
    + permitividades y permeabilidades relativas de las 2 subregiones
    + Condiciones de frontera
    + Opciones para los archivos de salida


El programa entrega:
    
    Gráfico de los resultados
    archivo .txt con resultados
    archivo .vtk con resultadoS
"""

import FEM2D_modulo as fcn
import numpy as np
import matplotlib.pyplot as plt
import meshio
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

    
#%% ENTRADAS

mesh_path = 'mallas_msh/planos5.msh' # path del archivo .msh
lado_hor = 5.0 #Longitud horizontal del dominio rectangular
lado_ver = 5.0 #Longitud vertical del dominio rectangular

sol_id = 'NA' #cadena al final del nombre del archivo que guarda el programa

w = 2975752870.9460554 # frecuencia angular


save_slns = False # activar para generar archivos .txt y .vtu

# permeabilidad y permitividad relativa de 2 regiones

mr1 = 1
er1 = 1

mr2 = 1
er2 = 4

# activadores de las condiciones de dirichlet, si inactivas, se soluciona
# para condiciones de Neumann homogéneas

Dirichlet_izq = True 
Dirichlet_der = True
Dirichlet_arr = True
Dirichlet_aba = True

# condiciones de dirichlet para los 4 bordes

diri_condicion_izq = 1
diri_condicion_der = 0
diri_condicion_arriba = 0
diri_condicion_abajo = 0


def fuente(x, y): #función fuente
#    f = np.exp(-30*((x-2.4)**2 + (y-2.5)**2)) # guasiana en 2.5, 2.5
    f = 0
    return f

#%% FUNCION DE SOLUCION

def main(mesh_path,sol_id,w,er1,mr1,er2,mr2, fuente, diri_condicion_izq,
         diri_condicion_der, diri_condicion_arriba, diri_condicion_abajo,
         lado_hor, lado_ver , Dirichlet_izq = True, Dirichlet_der = True,
         Dirichlet_arr = True, Dirichlet_aba = True, save_slns=False):
    
    c_const = 299792458 # velocidad de la luz en el vacío
    
    wavelength = 2*np.pi*c_const/w #longitud de onda en el vacío
    
    # LECTURA ARCHIVO .MSH
    
    puntos, elementos, permit_elem,permea_elem, npuntos, nelementos = fcn.leer_msh_imponer_materiales(mesh_path, mr1, er1, mr2, er2)
    
    k_onda = (2*np.pi/wavelength)
    # HALLANDO MATRICES LOCALES Y ENSAMBLAJE
    
    A_g = lil_matrix((npuntos, npuntos))
    M_g = lil_matrix((npuntos, npuntos))
    b_g = np.zeros(npuntos)
    
    for p in range(nelementos):
        A_l = (1 / permea_elem[p] ) * fcn.hallar_matriz_rigidez_local(elementos[p], puntos)
        M_l = ((k_onda**2) * permit_elem[p])* fcn.hallar_matriz_masa_local(elementos[p], puntos)
        b_l = -fcn.hallar_vector_cargas_local(elementos[p], puntos, fuente)
        A_g, M_g, b_g = fcn.ensamble(elementos[p], A_l, M_l, b_l, A_g, M_g, b_g)
    
    K = lil_matrix((npuntos, npuntos))
    
    K = A_g - M_g
    
    # APLICANDO CONDICIONES DE DIRICHLET FORZADAS EN LOS BORDES
    
    lado_izq = np.where(puntos[:,0] == 0.0)[0]
    lado_der = np.where(puntos[:,0] == lado_hor)[0]
    lado_arriba = np.where(puntos[:,1] == lado_ver)[0]
    lado_abajo = np.where(puntos[:,1] == 0.0)[0]
    
    if Dirichlet_der :
        for index in lado_der:
            b_g[index] = diri_condicion_der
            K[index] = np.zeros(npuntos)
            K[index, index] = 1
    
    if Dirichlet_arr:
        for index in lado_arriba:
            b_g[index] = diri_condicion_arriba
            K[index] = np.zeros(npuntos)
            K[index, index] = 1
        
    if Dirichlet_aba:
        for index in lado_abajo:
            b_g[index] = diri_condicion_abajo
            K[index] = np.zeros(npuntos)
            K[index, index] = 1
        
    if Dirichlet_izq:
        for index in lado_izq:
            b_g[index] = diri_condicion_izq
            K[index] = np.zeros(npuntos)
            K[index, index] = 1
    
    
    # RESOLVIENDO LOS VALORES EN LOS NODOS
    
    U = spsolve(K, b_g)
    
    # GRAFICACION Y ALMACENAMIENTO DE RESULTADOS
    
    datos = np.copy(puntos)
    datos[:,2] = U
    
    fig1, ax1 = plt.subplots()
    tcf = ax1.tricontourf(datos[:, 0], datos[:,1], datos[:,2], 8, 
                          alpha = 0.75, cmap = "seismic")
    fig1.colorbar(tcf)
    
    if save_slns:
    
        np.savetxt('datos_python/python' + sol_id + '.txt' , datos)
        
        points = puntos
        cells = {'triangle': elementos}
        mesh2 = meshio.Mesh(points, cells)
        mesh2.point_data["Ez"] = U
        mesh2.cell_data['triangle'] = {"permit_elem": permit_elem}
        mesh2.cell_data['triangle'] = {"permea_elem": permea_elem}
        meshio.write("archivos_vtk/python" + sol_id + ".vtu", mesh2)
      
    return U

#%% PRINCIPAL

if __name__ == '__main__':
    import timeit
    tic = timeit.default_timer()
    
    U = main(mesh_path,sol_id,w,er1,mr1,er2,mr2, fuente, diri_condicion_izq, 
             diri_condicion_der, diri_condicion_arriba, diri_condicion_abajo,
             lado_hor, lado_ver , Dirichlet_izq, Dirichlet_der, Dirichlet_arr, 
             Dirichlet_aba, save_slns)
    
    toc = timeit.default_timer()
    print('tiempo en segundos que tardó la simulación: ', toc - tic)  
