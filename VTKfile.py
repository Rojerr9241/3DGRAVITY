#               UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
#       POSGRADO EN CIENCIAS DE LA TIERRA - INSTITUTO DE GEOFÍSICA
#  CAMPO 2 : EXPLORACIÓN, AGUAS SUBTERRANEAS, MODELACIÓN Y PERCEPCIÓN REMOTA
#
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#      'VTKfile.py' . Programa para generar un archivo VTK "Paraview"
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                                               'Última modificación: 05/11/22'
#
#////////////////////////////////////////////////////////////////////////////////

import numpy as np
import pandas as pd

#Lectura del archivo txt
data=pd.read_csv('m_LBP2.dat',header=None)

#NÚMERO DE PRISMAS EN LOS EJES
Nx=45    #No. cubos eje x
Ny=30    #No. cubos eje y
Nz=40    #No. cubos eje z

N=Nx*Ny*Nz     #Número de prismas
N9=N*9         #N.prismas x 9 columnas
Nstr=str(N)    #Convertir a string
N9str=str(N9)

#NÚMERO DE NODOS EN CADA EJE
Nnx=Nx+1; #No. puntos eje x (vertice - vertice) vista perpendicular al eje
Nny=Ny+1; #No. puntos eje y
Nnz=Nz+1; #No. puntos eje z

#DISTANCIA CUBIERTA POR LOS PRISMAS
# Distancia en [m]
Lx=45
Ly=30
Lz=40

NN=Nnx*Nny*Nnz #Número de nodos totales
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
            

# x1 = np.array(G)
# G = np.array(G)
# G = np.array(G)

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

# INICIA ESCRITURA ARCHIVO VTK

lines = ['# vtk DataFile Version 2.0', 'Unstructured Grid Gz', 'ASCII', \
         'DATASET UNSTRUCTURED_GRID']
with open('file.vtk', 'w') as fi:
    for line in lines:
        fi.write(line)
        fi.write('\n')
    fi.write('\n')
    fi.write('POINTS ' + NNstr + ' float\n')

# with open('file.vtk', 'a') as f:
#     f.write('\n'.join(map(str, XYZ)))


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

with open('file.vtk', 'a') as fi:
    np.savetxt(fi,XYZ,fmt='%1.5f',delimiter='\t')
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



