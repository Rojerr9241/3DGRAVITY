!                  UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
!          FACULTAD DE INGENIERÍA-DIVISIÓN DE CIENCIAS DE LA TIERRA
!                       DEPARTAMENTO DE GEOFÍSICA
!_______________________________________________________________________________
! TESIS "INVERSIÓN POR RECRISTALIZACIÓN SIMULADA DEL CAMPO VECTORIAL-TENSORIAL
!            GRAVITACIONAL PARA LA EXPLORACIÓN DE YACIMIENTOS SUBSALINOS"
!_______________________________________________________________________________
!
! >> Programa que calcula el kernel de sensitividad de un ensamble de prismas.
!
!  ** El programa hace uso de la función gbox (Blakely,R.J.,1995; Nava-Flores,M.,2018)
!     donde se tiene que especificar el kernel que se quiere obtener. (Gx,Gy,Gz)
!
!  ** El resultado de este programa, es una matriz G que se guarda en el archivo
!     "KNLGx,KNLGy,KNLGz".dat. Los datos observados de la distribución de se
!      guarda en el archivo dobs"Gx,Gy,Gz".dat.
!
!                                        -> Última Revisión: 05 de noviembre del 2022
!                                              por Rodrigo Negrete Juárez
!
!_______________________________________________________________________________

program KNL_G

  use, intrinsic :: iso_fortran_env, dp=>real64, il=>int32

integer(IL)::Nobsx,Nobsy,i,j,k,r,s,k0,kd,Col,Ren,UnLec,UnEsc,Nx,Ny,Nz,Nnx,Nny,Nnz
real(DP):: z0,gk,gr,dx,dy,dPx,dPy,dPz,Lx,Ly,Lz,t1,t2
real(DP), allocatable:: x0(:),y0(:),xP(:),yP(:),zP(:),rho(:),G(:,:),Gt(:,:),&
                         f(:),gc(:,:)
character (len=3):: V,ref_sal
character(len=100):: DirRes


!*******************************************************************************
!*                        >> Interfaz con la función <<
!*******************************************************************************
!   >> Interfaz para una correcta comunicación con la función.

  interface
    function gboxVEC(x0,y0,z0,x1,y1,z1,x2,y2,z2,rho,V)

          use, intrinsic :: iso_fortran_env, dp=>real64, il=>int32

      character (len=3), intent(in):: V
      real(DP),intent(in):: x0,y0,z0,x1,y1,z1,x2,y2,z2,rho
      real(DP):: gboxVEC

    end function
  end interface

!*******************************************************************************

  !Limpieza de terminal:
    call system('clear')

  print*, 'Programa que calcula el Kernel de sensitividad (Gx, Gy, Gz)'

  ! >> 1.Plano de observación (x0,y0,z0) [Nobsy,Nobsx]

    ref_sal='G1' !Referencia del archivo de salida
    V='Gz' ! Seleccionar la componente que se requiere (Gx, Gy, Gz)

    Nx=10 !No. de prismas en X
    Ny=10 !No. de prismas en Y
    Nz=10 !No. de prismas en Z

    Nobsx = 101  !No. de observaciones en el eje X
    Nobsy = 101  !No. de observaciones en el eje Y
    Ren=Nobsx*Nobsy !No. de renglones de la matriz G
    Col=Nx*Ny*Nz    !No. de columnas de la matriz G

    Lx = 0.1_DP !Longitud del dominio en x
    Ly = 0.1_DP !Longitud del dominio en y
    Lz = 0.1_DP !Longitud del dominio en z

    allocate(x0(Nobsx),y0(Nobsy),G(Ren,Col))

    dx = Lx/DBLE(Nobsx-1._IL) !Delta x malla de observación
    dy = Ly/DBLE(Nobsy-1._IL) !Delta y malla de observación

    !PLANO DE OBSERVACIÓN
    x0=[(i-1._DP, i=1,Nobsx)]*dx
    y0=[(i-1._DP, i=1,Nobsy)]*dy
    z0=-0.0001_DP ! Plano de observación a 1 centímetro de altura evita un ERROR

    ! >> 2.Ensamble de prismas

    Nnx = Nx+1 !No. de nodos en x
    Nny = Ny+1 !No. de nodos en y
    Nnz = Nz+1 !No. de nodos en z

    dPx = Lx/DBLE(Nx) !Delta x prismas
    dPy = Ly/DBLE(Ny) !Delta y prismas
    dPz = Lz/DBLE(Nz) !Delta z prismas

    allocate(xP(Nnx),yP(Nny),zP(Nnz))

    xP=[(i-1._IL, i=1,Nnx)]*dPx      !Vector de prismas en el eje x
    yP=[(i-1._IL, i=1,Nny)]*dPy      !Vector de prismas en el eje y
    zP=[(i-1._IL, i=1,Nnz)]*dPz      !Vector de prismas en el eje z


  ! ! >> 3.Lectura del vector de densidad [MODELO SINTÉTICO], paso opcional
  !     allocate(rho(Col))
  !
  !   DirRes='/home/rodrigo/Escritorio/TESIS2/4to/RHO/'
  !     open(NewUnit=UnLec, file=trim(DirRes)//'G1.dat', status='old', action='read')
  !       do i=1,Col
  !          read(UnLec,*) rho(i)
  !       end do
  !     close(UnLec)

  !Inicia cronómetro:
    call cpu_time(t1)

    print*, 'INICIA CÁLCULO DEL KERNEL G'

!*******************************************************************************
!*                         >>  KERNEL DE SENSITIVIDAD  <<
!*******************************************************************************

   k0=0 !Contador externo

      do r=1,Nobsx
        do s=1,Nobsy

           k0=k0+1
           kd=0 !Contador interno
           gr=0._DP

             do k=1,Nz
               do j=1,Ny
                 do i=1,Nx

                   kd=kd+1

                   gk=gboxVEC(x0(r),y0(s),z0,xP(i),yP(j),zP(k),xP(i+1),yP(j+1),zP(k+1),1._DP,V)
                   G(k0,kd)=gk


                 end do
               end do
             end do

        end do
      end do


  !Inicia cronómetro:
    call cpu_time(t2)

    print *, ' Tiempo de cómputo cpu_time: ', t2-t1


  !!!! ** Multiplicación del kernel con el vector de densidad
  !   f = matmul(G,rho)

    Gt = transpose(G) !Matriz transpuesta

!*******************************************************************************
!*                          >> ALMACENAMIENTO DATOS  <<                        *
!*******************************************************************************
! >>  Almacenamiento del Kernel de sensitividades en archivos ASCII con, formato
!     de tabla (matriz).

!!! KERNEL Transpuesto
    open(NewUnit=UnEsc,file='Knlt'//trim(ref_sal)//''//trim(V)//'.dat',status='replace',action='write')
      do i=1,Col
         do j=1,Ren
            write(UnEsc,'(f15.9,2x)',advance='no') Gt(i,j)
         end do
         write(UnEsc,*)
      end do
    close(UnEsc)

!!! KERNEL
    open(NewUnit=UnEsc,file='Knl'//trim(ref_sal)//''//trim(V)//'.dat',status='replace',action='write')
      do i=1,Ren
         do j=1,Col
            write(UnEsc,'(f15.9,2x)',advance='no') G(i,j)
         end do
         write(UnEsc,*)
      end do
    close(UnEsc)

! !.............................................................................
! ! >>  Almacenamiento de datos en archivos ASCII en forma de columna:
! !.............................................................................
!
!     open(NewUnit=UnEsc,file='dobs'//trim(ref_sal)//''//trim(V)//'.dat',status='replace',action='write')
!       do i=1,Ren
!          write(UnEsc,'(5(f15.9,2x))') f(i)
!       end do
!     close(UnEsc)
!
!
! !.............................................................................
! ! >>  Almacenamiento de resultados en archivos ASCII con formato
! !       de tabla (matriz):
! !.............................................................................
!
!   ! Ciclo para convertir la anomalía que esta en forma vector, a una matriz (Nobsx,Nobsy)
!   allocate(gc(Nobsx,Nobsy))
!   k0=0
!
!    do i=1,Nobsx
!      do j=1,Nobsy
!
!        k0=k0+1
!        gc(j,i)=f(k0)
!
!      end do
!    end do
!
!   ! Almacenamiento de datos en forma de Matriz
!    open(NewUnit=UnEsc,file=trim(ref_sal)//''//trim(V)//'.dat',status='replace',action='write')
!    do i=1,Nobsx
!       do j=1,Nobsy
!          write(UnEsc,'(f15.9,2x)',advance='no') gc(i,j)
!       end do
!       write(UnEsc,*)
!    end do
!    close(UnEsc)

   print*, 'PROGRAMA FINALIZADO CON ÉXITO'

end program KNL_G
