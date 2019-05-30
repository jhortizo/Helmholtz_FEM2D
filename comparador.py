"""
Archivo para comparar numéricamente las soluciones de COMSOL y Python

Se ingresan los paths de los archivos .txt de python y comsol, previamente ordenados



"""
import numpy as np

python_path = 'datos_python/python14.txt'

comsol_path = 'simulaciones_comsol/datos_slnes/comsol14.txt'

python_data = np.loadtxt(python_path)

comsol_data = np.loadtxt(comsol_path)

error_x =python_data[:,0]-comsol_data[:,0]
error_x = np.linalg.norm(error_x) / np.size(error_x)

if error_x > 1:
    print(' hay errores de ordenamiento en x ')
    print(python_data[:,0]-comsol_data[:,0])
    
else:
    error_y =python_data[:,1]-comsol_data[:,1]
    error_y = np.linalg.norm(error_y) / np.size(error_y)
    if error_y > 1:
        print(' hay errores de ordenamiento en y ')
        print(python_data[:,1]-comsol_data[:,1])
    else: 
        error_Ez =python_data[:,2]-comsol_data[:,2]
        error_Ez = np.linalg.norm(error_Ez) / np.size(error_Ez)
        print (' error de los campos : ' , error_Ez)
        print(python_data[:,2]-comsol_data[:,2])
    
