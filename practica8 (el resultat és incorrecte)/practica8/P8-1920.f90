! ----------------------------------- Pràctica 8 --------------------------------------- !
! Autor: Javier Rozalén Sarmiento
! Grup: B1B
! Data: 03/12/2019
!
! Funcionalitat: es resol l'equació de Scdxrödinger independent del temps per trobar els
! autovalors i autovectors d'una partícula en una caixa 1D de mida L emprant el mètode 
! artiller amb Runge-Kutta 3
!
! Comentaris: la probabilitat s'hauria de calcular cridant la subrutina trapezoids amb un
! N adequat, però no s'ha fet per falta de temps. Les figures també s'haurien d'haver
! representat amb un format dels nombres de l'eix Y adequat.
!
! Nota: el programa principal i les subrutines tenen indentacions tals que, en alguns
! editors de text (SublimeText), es poden plegar i desplegar per facilitar la lectura.

program practica8
    implicit none
    double precision V_0,E,alfa,delta,a,dx,x
    double precision integral,L,E1,E2,beta
    double precision, allocatable :: vectphi(:),yyin(:),yyout(:)
    integer i,j,N,nequs
    common/cts/V_0,E,alfa,delta,a
    common/b/beta

    nequs=2
    allocate(yyin(2))
    allocate(yyout(2))
    L=8.d0
    beta=0.d0
    alfa=2.d0
    delta=0.05d0
    V_0=-20.d0
    a=7.6199d0

    ! -------------------------------- Apartat 1 --------------------------------------- !
    N=400 ! Nombre de passos per Ralston
    dx=L/dble(N)
    open(12,file="aux.dat") ! Funcions d'ona apartat 1

    ! ---------------- E1  ---------------- !
    E=-21.d0
    yyin(1)=0.d0
    yyin(2)=0.000002d0
    do i=1,N
        x=-L/2.d0+i*dx
        call RLSTN3(x,dx,nequs,yyin,yyout)
        write(12,*) x,yyout(1)
        yyin=yyout
    enddo
    call write(12)

    ! ---------------- E2  ---------------- !
    E=-20.5d0
    yyin(1)=0.d0
    yyin(2)=0.000002d0
    do i=1,N
        x=-L/2.d0+i*dx
        call RLSTN3(x,dx,nequs,yyin,yyout)
        write(12,*) x,yyout(1)
        yyin=yyout
    enddo
    call write(12)

    ! ---------------- E3  ---------------- !
    E=-14.d0
    yyin(1)=0.d0
    yyin(2)=0.000002d0
    do i=1,N
        x=-L/2.d0+i*dx
        call RLSTN3(x,dx,nequs,yyin,yyout)
        write(12,*) x,yyout(1)
        yyin=yyout
    enddo
    call write(12)

    ! ---------------- E4  ---------------- !
    E=-13.d0
    yyin(1)=0.d0
    yyin(2)=0.000002d0
    do i=1,N
        x=-L/2.d0+i*dx
        call RLSTN3(x,dx,nequs,yyin,yyout)
        write(12,*) x,yyout(1)
        yyin=yyout
    enddo
    call write(12)

    ! -------------------------------- Apartat 2 --------------------------------------- !
    N=400
    allocate(vectphi(N))
    open(13,file="aux2.dat") ! Convergència del mètode
    open(14,file="aux3.dat") ! Autovectors normalitzats

    ! Càlcul del primer autovalor
    E1=-21.d0
    E2=-20.5d0
    call artiller(E1,E2,nequs,N,vectphi,13) 
    call write(13)
    call trapezoids(-L/2.d0,L/2.d0,N,abs(vectphi)**2.d0,integral)
    vectphi=vectphi/dsqrt(integral)
    write(14,*) 0.d0,0.d0
    ! Escrivim el vector propi en un arxiu
    do i=1,N
        write(14,*) i*dx,vectphi(i)
    enddo
    call write(14)

    ! Càlcul del segon autovalor
    E1=-14.d0
    E2=-13.d0
    call artiller(E1,E2,nequs,N,vectphi,13)
    call write(13)
    call trapezoids(-L/2.d0,L/2.d0,N,abs(vectphi)**2.d0,integral)
    vectphi=vectphi/dsqrt(integral)
    write(14,*) 0.d0,0.d0
    ! Escrivim el vector propi en un arxiu
    do i=1,N
        write(14,*) i*dx,vectphi(i)
    enddo
    call write(14)

    ! Càlcul del tercer autovalor
    E1=-8.d0
    E2=-7.5d0
    call artiller(E1,E2,nequs,N,vectphi,13)
    call trapezoids(-L/2.d0,L/2.d0,N,abs(vectphi)**2.d0,integral)
    vectphi=vectphi/dsqrt(integral)
    write(14,*) 0.d0,0.d0
    ! Escrivim el vector propi en un arxiu
    do i=1,N
        write(14,*) i*dx,vectphi(i)
    enddo

    ! -------------------------------- Apartat 3 --------------------------------------- !
    E1=-21.d0
    E2=-20.5d0
    open(15,file="aux4.dat")
    open(16,file="P8-1920-res1.dat")

    ! Primer valor per beta
    beta=0.d0
    call artiller(E1,E2,nequs,N,vectphi,15) 
    call trapezoids(-L/2.d0,L/2.d0,N,abs(vectphi)**2.d0,integral)
    vectphi=vectphi/dsqrt(integral)
    write(15,*) 0.d0,0.d0
    ! Escrivim el vector propi en un arxiu
    do i=1,N
        write(15,*) i*dx,vectphi(i)
    enddo
    call write(15)
    ! Càlcul de la probabilitat de trobar l'electró entre -1.3,1.3
    call trapezoids(-1.3d0,1.3d0,N,abs(vectphi)**2.d0,integral)
    write(16,*) "Probabilitat beta=0: ",integral
    call write(16)

    ! Segon valor per beta
    beta=1.d0
    call artiller(E1,E2,nequs,N,vectphi,15) 
    call trapezoids(-L/2.d0,L/2.d0,N,abs(vectphi)**2.d0,integral)
    vectphi=vectphi/dsqrt(integral)
    write(15,*) 0.d0,0.d0
    ! Escrivim el vector propi en un arxiu
    do i=1,N
        write(15,*) i*dx,vectphi(i)
    enddo
    call write(15)
    ! Càlcul de la probabilitat de trobar l'electró entre -1.3,1.3
    call trapezoids(-1.3d0,1.3d0,N,abs(vectphi)**2.d0,integral)
    write(16,*) "Probabilitat beta=1: ",integral
    call write(16)

    ! Tercer valor per beta
    beta=5.d0
    call artiller(E1,E2,nequs,N,vectphi,15) 
    call trapezoids(-L/2.d0,L/2.d0,N,abs(vectphi)**2.d0,integral)
    vectphi=vectphi/dsqrt(integral)
    write(15,*) 0.d0,0.d0
    ! Escrivim el vector propi en un arxiu
    do i=1,N
        write(15,*) i*dx,vectphi(i)
    enddo
    call write(15)
    ! Càlcul de la probabilitat de trobar l'electró entre -1.3,1.3
    call trapezoids(-1.3d0,1.3d0,N,abs(vectphi)**2.d0,integral)
    write(16,*) "Probabilitat beta=5: ",integral
    call write(16)

    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
end program practica8

! Subrutina RLSTN3 --> Calcula un pas del mètode de Ralston de tercer ordre per un sistema
! d'n equacions de primer ordre acoblades
subroutine RLSTN3(x,dx,nequs,yyin,yyout)
    ! x --> Variable independent del problema
    ! dx --> Pas 
    ! nequs --> Nombre d'equacions 
    ! yyin --> Vector amb les dades del punt anterior
    ! yyout --> Vector amb les dades del punt següent
    implicit none
    integer nequs,i
    double precision x,dx,yyin(nequs),yyout(nequs)
    double precision k1(nequs),k2(nequs),k3(nequs)
    double precision dyout(nequs) ! Vector mut necessari per cridar la subrutina derivades

    ! Càlcul dels vectors k1,k2,k3 
    call derivades(nequs,x,yyin,dyout)
    k1=dyout
    call derivades(nequs,x+dx/2.d0,yyin+dx/2.d0*k1,dyout)
    k2=dyout
    call derivades(nequs,x+dx*0.75d0,yyin+(3.d0*dx/4.d0)*k2,dyout)
    k3=dyout

    ! Càlcul del vector yyout 
    yyout=yyin+dx/9.d0*(2.d0*k1+3.d0*k2+4.d0*k3)

    return 
end subroutine RLSTN3

! Subrutina derivades --> Calcula la derivada de la funció (vectorial) a trobar
subroutine derivades(nequ,x,yin,dyout)
    ! nequ --> nombre d'equacions
    ! x --> Variable independent del problema
    ! yin --> Vector amb les dades del punt on es calcula la derivada
    ! dyout --> Vector amb les derivades
    implicit none
    integer nequ
    double precision x,yin(nequ),dyout(nequ)
    double precision V_0,E,alfa,delta,a,beta,V
    common/cts/V_0,E,alfa,delta,a
    common/b/beta

    dyout(1)=yin(2)
    dyout(2)=(2.d0/a)*(V(x)+beta*x**2.d0-E)*yin(1)

    return
end subroutine derivades

! Subrutina artiller --> Resol una equació diferencial amb condicions de contorn pel mètode de tir
subroutine artiller(E1,E2,nequs,npassos,vectphi,arxiu)
    ! nequs,npassos --> nombre d'equacions del sistema a resoldre i nombre de passos pel solver
    ! E1,E2 --> Valors inicials aleatoris per inicialitzar el mètode de la secant
    ! yyin,yyout --> vectors amb les conidions inicials i resultats finals
    ! vectphi --> vector amb la solució a l'equació 
    ! arxiu --> arxiu on s'escriuen els valors d'energia
    implicit none
    integer i,j
    double precision E1,E2,yyin(nequs),yyout(nequs)
    double precision V_0,E,alfa,delta,a
    integer nequs,npassos,N,arxiu
    double precision phiE1,phiE2,phiE3,vectphi(npassos),E3,dx,x
    common/cts/V_0,E,alfa,delta,a

    N=npassos
    dx=1.d0/dble(N)

    do j=1,10000 
        ! Trobem una primera solució per E=E1
        yyin(1)=0.d0
        yyin(2)=2.000006d0
        E=E1
        do i=1,N
            call RLSTN3(x,dx,nequs,yyin,yyout)
            vectphi(i)=yyout(1)
            ! Reescrivim variables
            yyin=yyout
        enddo
        phiE1=vectphi(N)

        ! Trobem una segona solució per E=E2
        yyin(1)=0.d0
        yyin(2)=2.000006d0
        E=E2
        do i=1,N
            call RLSTN3(x,dx,nequs,yyin,yyout)
            vectphi(i)=yyout(1)
            ! Reescrivim variables
            yyin=yyout
        enddo
        phiE2=vectphi(N)

        ! Amb una nova E calculada amb el mètode de la secant, trobem una solució 
        ! aproximada
        yyin(1)=0.d0
        yyin(2)=2.000006d0
        E3=(E1*phiE2-E2*phiE1)/(phiE2-phiE1)
        E=E3
        do i=1,N
            call RLSTN3(x,dx,nequs,yyin,yyout)
            vectphi(i)=yyout(1)
            ! Reescrivim variables
            yyin=yyout
        enddo
        phiE3=vectphi(N)
        write(arxiu,*) E3,phiE3 

        if (abs(phiE3).lt.0.00006d0) then
            exit
        else
            E1=E2
            E2=E3
        endif 
    enddo

    return
end subroutine artiller

double precision function V(x)
    implicit none
    double precision x,V_0,E,alfa,delta,a
    common/cts/V_0,E,alfa,delta,a

    V=V_0*dsinh(alfa/delta)/(dcosh(alfa/delta)+dcosh(x/delta))

    return
end function V

! Subrutina trapezoids --> Calcula una integral 1-D per trapezis
subroutine trapezoids(x1,x2,ndim,funci,integral)
    ! x1,x2 --> Punts inicial i Final 
    ! ndim --> Nombre de dimensions del vector funció
    ! funci --> Vector funció (conté totes les imatges)
    ! integral --> Valor de la integral a calcular
    implicit none
    double precision x1,x2,integral,funci(ndim)
    integer i,ndim

    integral = 0.d0
    do i=1,ndim-1
        integral = integral + funci(i)
    enddo
    integral = (integral + funci(1)/2.d0 + funci(ndim)/2.d0)*((x2-x1)/dble(ndim))

    return
end subroutine trapezoids

! Subrutina write --> Escriu dues línies en blanc en un arxiu
subroutine write(arxiu)
    ! arxiu --> número de l'arxiu
    implicit none
    integer arxiu

    write(arxiu,*) ""
    write(arxiu,*) ""

    return
end subroutine