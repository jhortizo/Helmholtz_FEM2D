Implementación en Python de la solución de la ecuación de Helmholtz con permeabilidades
y permitividades inhomogéneas para 2 dimensiones,

Hecha por : Sofia Obando Vasquez, Jose Hernan Ortiz, Jose Manuel Rendón, Mayo 2019
Proyecto final de la materia Metodos numericos, a cargo del profesor Nicolas Guarin-Zapata
2019-1, Universidad EAFIT 

CONTIENE LOS PROGRAMAS :

	+FEM2D_main.py :archivo principal, soluciona la ecuación de Helmholtz para parámetros del medio
	variables. Se alimenta de las funciones de FEM2D_modulo.py. Arroja la solución, una gráfica
	de la misma, y archivos .vtu y .txt si se quiere.

		entrada : path al archivo .msh
			  frecuencia angular de la onda
			  opciones de guardado de archivo de salida
			  opciones de condiciones de frontera
			  función fuente
			  permeabilidades y permitividades relativas de 2 medios

		salida : gráfica de la solución
			 archivos .vtu y .txt con al solución

	+FEM2D_modulo.py : programa que contiene las funciones necesarias para la solución de la ecuación por 
	el método de elementos finitos. 

	+identificador_subregiones.py : programa que permite la identificación de las subregiones definidas
	en un archivo .msh , dando el orden de identificación de dichas subregiones. Debe usarse en caso 
	de dudas en la definición de los grupos físicos del archivo .msh

	+nodos_ordenador.py : programa que reorganiza los archivos .txt arrojados por los solucionadores,
	para ordenar los nodos de forma análoga en estos y permitir la comparación nodo a nodo

	+ comparador.py : programa en Python que compara los archivos .txt de las soluciones en Python y en 
	COMSOL, y arroja el error entre las dos. Antes de usarlo, asegúrese de ordenar los datos de los
	archivos .txt con nodos_ordenador.py

		entrada : path al archivo .txt de python
			  path al archivo .txt de COMSOL

		salida : error entre las soluciones

RESTRICCIONES DEL PROGRAMA:

	+actualmente solo trabaja para 2 subregiones distinguibles en el dominio a solucionar.
	estas dos subregiones deben tener parámetros de material constantes. dominios con una cantidad de particiones
	mayor generarán errores en la solución

	+las condiciones de dirichlet se imponen correctamente solo cuando se trabaja con dominios
	rectangulares. Dentro de las entradas del programa se encuentran dos variables llamadas 'lado_hor' y 'lado_ver'
	para que el usuario especifique las dimensiones del rectángulo. Sin embargo, no hay restricciones en la forma
	de las subregiones en las que se parte el dominio.

	+Solo se implementaron condiciones de frontera de dirichlet (con cualquier valor) y las de neumann
	(homogeneas), las condiciones de Robin no se pueden simular con el programa.

CONTIENE LAS CARPETAS :

	+ archivos_vtk = contiene los archivos .vtk / .vtu con los datos de las simulaciones hechas en python,
	a este directorio se guardan los archivos .vtu de salida del programa principal. Para ver en Paraview

	+ datos_python = contiene los archivos .txt con los datos de las simulaciones hechas en python,
	a este directorio se guardan los archivos .txt de salida del programa principal}

	+mallas_msh = contiene los archivos .msh con las diferentes mallas usadas para las simulaciones que se
	presentan en python. A continuación se enumeran las mallas presentes, donde planos#, hace referencia a
	mallas con subdominios rectangulares, circulo# hace referencia a mallas con un circulo en medio, y 
	corazon# hace referencia a mallas con un corazón en medio.

			planos1: contiene 558 nodos
			planos2: contiene 13 nodos
			planos3: contiene 2214 nodos
			planos4: contiene 3421 nodos
			planos5: contiene 6987 nodos
			planos6: contiene 20884 nodos
			circulo1: contiene 424 nodos
			circulo2: contiene 44879 nodos
			circulo3: contiene 3419 nodos
			corazon1: contiene 2947 nodos

	Todas las mallas definen dominios cuadrados entre 0 y 5. Todas elaboradas con el programa Gmsh

	+simulaciones_comsol = contiene diversas carpetas con los archivos relacionados con las simulaciones 
	que se hicieron en COMSOL Multiphysics para cotejar las simulaciones en Python. Contiene:

		-archivos_comsol: contiene los archivos de comsol con las simulaciones, guardados de acuerdo
		al numero de la simulación

		-archivos_vtu : contiene los archivos .vtu exportados desde comsol. Para ver el Paraview

		-datos_slnes : contiene los archivos .txt exportados desde comsol

		-mallas_nastran : contiene las mismas mallas que mallas_msh pero en formato .nas

SIMULACIONES LLEVADAS A CABO: (a menos que se diga lo contrario, el índice se controla con la permitividad relativa,
y la fuente es f = 0)

	1 : hecha con geometría planos2, índices de 1 y 1, condiciones de Dirichlet 1 a izquierda y 0 en el resto. 
Da un error de 0.09289521286686218. Tiempo= 0.0497 s tiempo comsol 4s

	2: hecha con geometría planos1, índices de 1 y 1, condiciones de Dirichlet 1 a izquierda y 0 en el resto. 
Da un error de 0.07977691790519442. Tiempo = 1.513 tiempo comsol 4s
	
	3: hecha con geometría planos3, índices de 1 y 1, condiciones de Dirichlet 1 a izquierda y 0 en el resto. 
Da un error de 0.06716363759121191 tiempo 7.268 tiempo comsol 5s

	4: hecha con geometría planos2, índices 1 y 1, condiciones de Dirichlet 1 a izquierda y 0 en el resto. 
Da un error de 0.06499157926054021. En los vértices donde se superponen las condiciones de dirichlet se hace
un promedio (osea que quedan el 0.5).

	5: hecha con geometría planos1, índices de 1 y 2, condiciones de Dirichlet 1 a izquierda y 0 en el resto. 
Da un error de 0.48901269096155037.

	6: hecha con geometría planos3, índices de 1 y 2, condiciones de Dirichlet 1 a izquierda y 0 en el resto. 
Da un error de 0.1714157640861453.

	7: hecha con geometría planos4, índices de 1 y 1, condiciones de Dirichlet 1 a izquierda y 0 en el resto. 
Da un error de 0.0421555302746676. tiempo = 12.668 tiempo comsol 5s

	8: hecha con geometría planos5, índices de 1 y 1, condiciones de Dirichlet 1 a izquierda y 0 en el resto. 
Da un error de 0.034807815068731776. tiempo = 35.22 tiempo comsol 9s

	9: hecha con geometría planos6, índices de 1 y 1, condiciones de Dirichlet 1 a izquierda y 0 en el resto. 
Da un error de 0.02299845894309353. tiempo = 230.18 tiempo comsol 20s

	10: hecha con geometria circulo1, indices 1 y 2, condiciones de dirichlet 1 a isq y 0 en el resto. 
Da un error 0.7143367552787476

	11: hecha con geometria circulo3, indices 1 y 2, condiciones de dirichlet 1 a isq y 0 en el resto. 
Da un error 0.03697256652029045

	12: hecha con geometria circulo2, indices 1 y 2, condiciones de dirichlet 1 a isq y 0 en el resto. 
Da un error 0.008261156758867112.

	13: hecha con geometria corazon1, indices 1 y 2, condiciones de dirichlet 1 a izq y 0 en el resto. 
Da un error 0.2734359232883579

	14. hecha con geometría planos5, índices de 1 y 2, condiciones de Dirichlet 1 a izquierda y 0 en el resto. 
Da un error de 0.009100016669280885

	especial1: hecha con geometría planos5, índices 1 y 2, condiciones de Dirichlet homogéneas, función fuente gaussiana

	especial2: hecha con geometría planos5, región 1 con permitividad = 4 y permeabilidad = 1, región 2 con
	permitividad = 1 y permeabilidad = 4
