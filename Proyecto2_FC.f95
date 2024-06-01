! Programa que encuentra las asíntotas, raíces y puntos críticos dentro de el intervalo de una función.

Program proyecto2
    implicit none
! Declaramos las variables de entrada 
! Declaramos la función
    real :: F
! variable de distancia total, variable de distancia recorrida
    real :: A, B, DT
    real :: y1, y2, y, DR, DX
! variables dummies para los ciclos 
    integer :: i, j
! Enteros para guardar el número de asintotas 
    integer :: nv, nh
! variable del número de pasos 
    integer :: N 

! variable logica para asintotas verticales
    logical :: condition

! Leemos los límites sobre los que buscaremos 
    write(*,*) 'Escriba los limites del dominio entre los que se buscara'
    read(*,*) A, B

! Obtenemos la distancia y la dividimos sobre el número de pasos para obtener la distancia que recorrerá cada paso 
    nv = 1
    nh = 1
    DT = Abs(B-A)
    N = 100000 
    DX = DT/N 
    write(*,*) 'Cada paso mide ', DX
    
! Con cada paso se calculará el valor de la función y el valor pasado de esta 
    DR = A
    open(1, file='puntos.txt')
    open(2, file='asintv.txt')
    open(3, file='asinth.txt') 

    Do while (i<=N) 
        y = F(DR) 
        y2 = F(DR + DX)
        y1 = F(DR - DX)

     ! Se establecen una serie de criterios para detectar asíntotas verticales
        IF (abs(y2-y1) > 1000) then
          write(*,*) 'Encontramos asintota vertical en', DR 
          nv = nv+1
          write(2,*) DR 
          j = 1
          condition = abs(f(DR+DX)-f(DR))>1000
          Do while (condition)
            DR=DR+10*DX 
            condition = abs(f(DR+DX)-f(DR))>1000
            j = j+1
            if (j > 50) then 
                exit 
            end if
          End do
        
        else if (abs(y2-y1)<=1.0E-02) then 
            ! se llama a una subrutina que aproxime las asíntotas horizontales
            call asin(DR, DX, i, N, nh) 
        End If

        IF (abs(y) < 1.0E-03 ) then
            ! Se escribe la raíz
            call raiz(DR, DX)
        End if

        DR = DR + DX
        i = i+1

    end do

    close(1)
    close(2)
    close(3)
! Graficamos la función con los valores encontrados de raiz
    call graph(nv, nh)

    call system('start gnuplot -p grafica.gpl')

    print*, 'espere a que se genere la grafica y de enter para abrirla'
    read(*,*)

    call system('Funcion.png')

end program 

subroutine Graph(nv, nh)
    integer :: nv, nh 
    real :: x_coord, x_ci, y_ci, x_cf, y_cf

    open(2, file='asintv.txt')
    open(3, file='asinth.txt') 

    open (439, file='grafica.gpl')
    write(439,*) "set terminal pngcairo enhanced font 'Verdana,12'"
    write(439,*) "set output 'Funcion.png'"
    write(439,*) "set title 'Función con sus puntos críticos, raíces y asíntotas'"
    write(439,*) "set xlabel 'x'"
    write(439,*) "set ylabel 'y'"
    write(439,*) "set xrange [-2:200]"
    write(439,*) "set yrange [-2:1]"
    write(439,*) "set border linewidth 1.5"
    write(439,*) 'set xzeroaxis linetype 1 linewidth 2.5'
    write(439,*) "set yzeroaxis linetype 1 linewidth 2.5"
    write(439,*) "set grid"

    if (nv > 1) then 
        do i=1, nv-1 
            read(2,*) x_coord
            write(439,*) "set arrow from  ",x_coord,",graph 0 to ",x_coord,",graph 1 lt 1 lc 'red' linewidth 5 nohead"
        end do
    end if
    if (nh > 1) then 
        do i=1, nh-1 
            read(3,*) x_ci, y_ci
            read(3,*) x_cf, y_cf
            write(439,*) "set arrow from ",x_ci,",",y_ci," to ",x_cf,",",y_cf," lt 1 lc 'red' linewidth 5 nohead"
        end do
    end if

    write(439,*) "plot x**2+40*x+1 with lines ls 1 lw 3 lc 'black', \" 
    write(439,*) "      'puntos.txt' with points ps 1 pointtype 7 pointcolor 'red' title 'puntos críticos y raíces'"

    close(439)
    close(2)
    close(3)
    RETURN
end subroutine Graph


! Subrutina para encontrar asíntotas horizontales o puntos críticos por medio de derivación numérica
Subroutine Refin (DR, DX)
    implicit none 
    real :: y, h, DR, deriv1, F, step, DX
    integer :: j
    h = 1.0E-005
    step = DX/1000

    Do j=1,1000
        deriv1 = (F(DR+h)-F(DR-h))/(2.0*h)
        write(*,*) 'usamos refin'

        If (abs(deriv1) <= h) then
            y = F(DR)
            ! Escribe las coordenadas del punto crítico
            write(1,*) DR, y
            write(*,*) 'Encontramos un punto critico en ', DR
            DR = DR+5*DX
            exit

        end if
        DR = DR + step
    end do

End SUBROUTINE

subroutine asin(DR, DX, i, N, nh)

    implicit none 
    real :: h, DR, DX, F, x1, x2, y, y2, y1
    
    logical :: asin_condition1, asin_condition2
    integer :: i, N, nh
    
    h = 1.0E-003
    x1 = DR + 1000 
    x2 = DR - 1000

    asin_condition1 = abs(f(DR)-f(x2))<=h !Detecta asintotas a la izquierda
    asin_condition2 = abs(f(DR)-f(x1))<=h !Detecta asintotas a la derecha

    If (asin_condition1) then
        y = f(DR)
        Do while (asin_condition1)
            DR=DR+100*DX 
            asin_condition1 = abs(f(DR)-f(x2))<=h
            write(*,*) asin_condition1, DR
        End do
        write(3,*) -100, y
        write(3,*) DR, y
        write(*,*) 'encontramos asintota horizontal a la izquierda desde el infinito hasta', DR 
        nh = nh+1
        ! Vamos a encontrar desde que coordenada x hay una asintota a la izquierda
        
    else if (asin_condition2) then 
        
        ! Encontramos el punto donde comienza una asintota a la derecha
        write(3,*) DR, f(DR)
        write(3,*) 100, f(DR)
        write(*,*) 'encontramos asintota horizontal a la derecha desde', DR, 'hasta infinito'
        nh=nh+1
        i=N+1
    end if

    ! Si no cumple ninguna de esas condiciones, es punto crítico
    y2 = F(DR + DX)
    y1 = F(DR - DX)
    IF (abs((y2-y1)/(2*DX)) < 1.0E-04) then
        ! Se llama a una subrutina q aproxime el punto crítico
        call Refin(DR, DX) 
    End if

end subroutine

Subroutine raiz(DR, DX)
    implicit none 
    real :: h, DR, F, DX, en, x0, xn, xs
    integer :: j
    h = 1.0E-06
    j = 0
    en = h+1
    x0=DR-3*DX
    xn=DR+3*DX
    write(*,*) 'usamos raiz'

    do while (en>h)
        xs=(x0+xn)/2
        if (f(x0)*f(xs)<0) then 
            xn=xs
        else 
            x0=xs
        end if
        en=abs((xs-xn)/xs)
        j=j+1
        If (j>=100) then
            exit 
        end if
    end do

    if (xn>0) then 
        write(1,*) xn, f(xn)
        write(*,*) 'Encontramos raiz en ', xn
        DR = DR + DX 
    else 
        write(1,*) x0, f(x0)
        write(*,*) 'Encontramos raiz en ', x0
        DR = DR + DX
    end if 
end subroutine

FUNCTION F(X)
    pi = 4 * atan(1.0)
    !F = sin(x)-x**2 
    F = x**2+40*x+1
    !F = log(x**2+1)-exp(x/2)*cos(pi*x)
    !F = tanh(x)
    !F = tan(x)
  RETURN
END FUNCTION