#               UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
#       POSGRADO EN CIENCIAS DE LA TIERRA - INSTITUTO DE GEOFÍSICA
#  CAMPO 2 : EXPLORACIÓN, AGUAS SUBTERRANEAS, MODELACIÓN Y PERCEPCIÓN REMOTA
#
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#      'VTKfile.py' . Programa para generar un archivo VTK "Paraview"
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                                               'Última modificación: 31/07/23'
#                                                Rodrigo Negrete Juárez
#/||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

import numpy as np
import pandas as pd

filename='filename' # 1. Nombre referencia del archivo final

 # Selecciona un caso [N,U,D]
command = 'N'

# Case [N] -> Se genera el archivo VTK de todo el ensamble, continen contrastes 
#             de densidad positivos y negativos.  
# Case [U] -> Se genera el archivo VTK de todo el ensamble, con los contrastes de
#             densidad positivos, mayores al límite superior LS.  
# Case [D] -> Se genera el archivo VTK de todo el ensamble, con los contrastes de
#             densidad negativos, menores al límite inferior LI.  

#Lectura del archivo txt
data=pd.read_csv('file.dat',header=None) #2. Tu vector de densidades 'data.txt' o 'data.dat'

# 3. (filtrado, límite superior e inferior) Modificar
LI=-1500 
LS=1200

# 4. Número de prismas del modelo
Nx=60    #No. prismas eje x
Ny=40    #No. prismas eje y
Nz=45    #No. prismas eje z

N=Nx*Ny*Nz   #Total de prismas

#NÚMERO DE NODOS EN CADA EJE
Nnx=Nx+1; #No. puntos eje x (vertice - vertice) vista perpendicular al eje
Nny=Ny+1; #No. puntos eje y
Nnz=Nz+1; #No. puntos eje z

# 5. Longitud [m]
Lx=45
Ly=30
Lz=40

NN=Nnx*Nny*Nnz #Total de nodos 
NNstr=str(NN)

x=np.zeros([NN,1])
y=np.zeros([NN,1])
z=np.zeros([NN,1])

# PRIMERA PARTE del archivo VTK (coordenadas de los vértices de los cubos)

#Contador para el ciclo
k=0

for cz in range(Nnz):
    for cy in range(Nny):
        for cx in range(Nnx):
            
            x[k]=cx
            y[k]=cy
            z[k]=cz
            k+=1
            

#Para sacar el delta d = (distancia real [km] )/(No. Cubos)
dx=Lx/Nx
dy=Ly/Ny
dz=Lz/Nz

#Convertimos a las unidades reales
x=x*dx
y=y*dy
z=-1*z*dz # Es a profundidad (-1)

XYZ=np.concatenate((x,y,z), axis=1)

print(XYZ[0]) #Coordenadas de los vértices XYZ

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Vértices de los cubos

# Salto de línea
shx=Nnx
shy=Nny
# Salto de plano
sv=Nnx*Nny

# Posición de los vértices dentro de las coordenadas
a=np.zeros([N,1])
b=np.zeros([N,1])
c=np.zeros([N,1])
d=np.zeros([N,1])
e=np.zeros([N,1])
f=np.zeros([N,1])
g=np.zeros([N,1])
h=np.zeros([N,1])

k2=0

for k in range(Nz):
    for j in range(Ny):
        for i in range(Nx):
            
            a[k2]=(i)+(shx*(j))+(sv*(k))             
            b[k2]=(i+1)+(shx*(j))+(sv*(k))           
            c[k2]=((shx-1)+(i+1))+(shx*(j))+(sv*(k)) 
            d[k2]=((shx)+(i+1))+(shx*(j))+(sv*(k))
            e[k2]=((sv-1)+(i+1))+(shx*(j))+(sv*(k))      
            f[k2]=((sv)+(i+1))+(shx*(j))+(sv*(k))
            g[k2]=((sv+shx-1)+(i+1))+(shx*(j))+(sv*(k)) 
            h[k2]=((sv+shx)+(i+1))+(shx*(j))+(sv*(k))  
            k2+=1

print(k2)

XYZ=np.concatenate((x,y,z), axis=1)

C1=np.ones([N,1])
Cp=C1*8
Ct=C1*11

CELL=np.concatenate((Cp,a,b,c,d,e,f,g,h), axis=1)

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# ARCHIVO VTK dependiendo el caso seleccionado

match command:
    case 'N':
        print('Hello to you too!')

        N9=N*9         #N.prismas x 9 columnas
        Nstr=str(N)    #Convertir a string
        N9str=str(N9)

        # INICIA ESCRITURA ARCHIVO VTK
        lines = ['# vtk DataFile Version 2.0', 'Unstructured Grid Gz', 'ASCII', \
        'DATASET UNSTRUCTURED_GRID']
        with open(filename +'.vtk', 'w') as fi:
            for line in lines:
                fi.write(line)
                fi.write('\n')
            fi.write('\n')
            fi.write('POINTS ' + NNstr + ' float\n')
            np.savetxt(fi,XYZ,fmt='%1.4f',delimiter='\t')
            fi.write('\n')
            fi.write('CELLS ' + Nstr +' '+ N9str+'\n')
            np.savetxt(fi,CELL,fmt='%d',delimiter='\t')
            fi.write('\n')
            fi.write('CELL_TYPES ' + Nstr +'\n')
            np.savetxt(fi,Ct,fmt='%2d',delimiter='\t')
            fi.write('\n')
            fi.write('CELL_DATA ' + Nstr +'\n')
            fi.write('SCALARS CONTRASTE_DENSIDAD float 1\n')
            fi.write('LOOKUP_TABLE default\n')
            np.savetxt(fi,data,fmt='%1.10f',delimiter='\t')


    case 'U':
        print('See you later')
        data.columns=['Rho']
    
        datau=data[data['Rho'] > LS]
        Ctu=Ct[data['Rho'] > LS]
        CELLu=CELL[data['Rho'] > LS]

        Nm=datau.size
        Nm9=Nm*9         #N.prismas condición x 9 columnas
        Nmstr=str(Nm)    #Convertir a string
        Nm9str=str(Nm9)

        # INICIA ESCRITURA ARCHIVO VTK
        lines = ['# vtk DataFile Version 2.0', 'Unstructured Grid Gz', 'ASCII', \
        'DATASET UNSTRUCTURED_GRID']
        with open(filename +'U.vtk', 'w') as fi:
            for line in lines:
                fi.write(line)
                fi.write('\n')
            fi.write('\n')
            fi.write('POINTS ' + NNstr + ' float\n')
            np.savetxt(fi,XYZ,fmt='%1.2f',delimiter='\t')
            fi.write('\n')
            fi.write('CELLS ' + Nmstr +' '+ Nm9str+'\n')
            np.savetxt(fi,CELLu,fmt='%d',delimiter='\t')
            fi.write('\n')
            fi.write('CELL_TYPES ' + Nmstr +'\n')
            np.savetxt(fi,Ctu,fmt='%2d',delimiter='\t')
            fi.write('\n')
            fi.write('CELL_DATA ' + Nmstr +'\n')
            fi.write('SCALARS CONTRASTE_DENSIDAD float 1\n')
            fi.write('LOOKUP_TABLE default\n')
            np.savetxt(fi,datau,fmt='%1.10f',delimiter='\t')

    case 'D':
        print('ajuaaaa')
        data.columns=['Rho']
    
        datau=data[data['Rho'] < LI]
        Ctu=Ct[data['Rho'] < LI]
        CELLu=CELL[data['Rho'] < LI]

        Nm=datau.size
        Nm9=Nm*9         #N.prismas condición x 9 columnas
        Nmstr=str(Nm)    #Convertir a string
        Nm9str=str(Nm9)

        # INICIA ESCRITURA ARCHIVO VTK
        lines = ['# vtk DataFile Version 2.0', 'Unstructured Grid Gz', 'ASCII', \
        'DATASET UNSTRUCTURED_GRID']
        with open(filename +'D.vtk', 'w') as fi:
            for line in lines:
                fi.write(line)
                fi.write('\n')
            fi.write('\n')
            fi.write('POINTS ' + NNstr + ' float\n')
            np.savetxt(fi,XYZ,fmt='%1.2f',delimiter='\t')
            fi.write('\n')
            fi.write('CELLS ' + Nmstr +' '+ Nm9str+'\n')
            np.savetxt(fi,CELLu,fmt='%d',delimiter='\t')
            fi.write('\n')
            fi.write('CELL_TYPES ' + Nmstr +'\n')
            np.savetxt(fi,Ctu,fmt='%2d',delimiter='\t')
            fi.write('\n')
            fi.write('CELL_DATA ' + Nmstr +'\n')
            fi.write('SCALARS CONTRASTE_DENSIDAD float 1\n')
            fi.write('LOOKUP_TABLE default\n')
            np.savetxt(fi,datau,fmt='%1.10f',delimiter='\t')

    case other:
        print('No match found')

print('Fin del programa!')
