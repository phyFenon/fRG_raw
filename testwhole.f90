program ktest
    !!!!---------------!!!!
    !    o1(k1,w1,s1)->---->o3(k3,w3,s3)
    !      |     |
    !    o2(k2,w2,s2)->---->o4(k4,w4,s4)
    !------------------!!!!
    implicit none
    integer,parameter::k_pn=3,nspin=2
    integer::k1,k2,k3,k4,o1,o2,o3,o4,k,oo1,oo2,oo3,oo4
    integer::kout,kindex
    complex,dimension(k_pn**3,k_pn**3,k_pn**3,nspin,nspin,nspin,nspin)::vert0,vert
    complex,dimension(k_pn**3,k_pn**3,k_pn**3,nspin,nspin,nspin,nspin)::pp,phc,phd
    real,dimension(k_pn**3,3)::k_point
    real::hubbu=8.0,lambda=30.0,vertmax=40
    real::dlambda,start,finish
    integer::step,stepmax=500
    complex::lpp,lphc,lphd
    

   dlambda=1.0
   call cpu_time(start)

   open(11,file='initialvalue of vertex.txt')

   ! to get the coarse kmesh
   call kmesh(k_pn,k_point)


   ! initial value for the vertex
   vert0=0.0
   vert0(:,:,:,1,2,1,2)=hubbu
   vert0(:,:,:,2,1,2,1)=hubbu
   vert0(:,:,:,1,2,2,1)=-hubbu
   vert0(:,:,:,2,1,1,2)=-hubbu

!to check and instore the initial vertex
! write(11,*)"#k1   k2   k3   o1   o2   o3   o4   V.real  "
! do k1=1,k_pn**3
!      do k2=1,k_pn**3
!         do k3=1,k_pn**3
!             do o1=1,nspin
!                 do o2=1,nspin
!                     do o3=1,nspin
!                          do o4=1,nspin
!                           write(11,*)k1,k2,k3,o1,o2,o3,o4,real(vert(k1,k2,k3,o1,o2,o3,o4))
!                           enddo
!                      enddo
!                   enddo
!               enddo
!           enddo
!       enddo
!     enddo                    
! close(11)
   
  
step=0
print*,"start"
   do k1=1,k_pn**3
    print*,k1
         do k2=1,k_pn**3
           do k3=1,k_pn**3
      
              do o1=1,nspin
                   do o2=1,nspin
                        do o3=1,nspin
                            do o4=1,nspin
                            pp=0.0
                            phd=0.0
                            phc=0.0
   !first do the k in BZ integration 
   !1/N SUM_k
   do k=1,k_pn**3
     kindex=kout(k_pn,k1,k2,k,k_point)
     
   ! Do the summation of sigma_i
    do oo1=1,nspin
        do oo2=1,nspin
            do oo3=1,nspin
                do oo4=1,nspin

     call lpp_cal(k_pn,k_point,k1,k2,k,oo1,oo2,oo3,oo4,lambda,lpp)

          pp(k1,k2,k3,o1,o2,o3,o4)=pp(k1,k2,k3,o1,o2,o3,o4)+0.5*vert0(k1,k2,k,o1,o2,oo1,oo2)&
                                  &*lpp*vert0(k,kindex,k3,oo3,oo4,o3,o4)

     call lph_cal(k_pn,k_point,k1,k3,k,oo1,oo2,oo3,oo4,lambda,lphd)
          kindex=kout(k_pn,k,k3,k1,k_point)

          phd(k1,k2,k3,o1,o2,o3,o4)=phd(k1,k2,k3,o1,o2,o3,o4)+vert0(k1,kindex,k3,o1,oo4,o3,oo1)&
                                   &*lphd*vert0(k,k2,kindex,oo3,o2,oo2,o4)
          
          k4=kout(k_pn,k1,k2,k3,k_point)

     call lph_cal(k_pn,k_point,k1,k4,k,oo1,oo2,oo3,oo4,lambda,lphc)
          
          kindex=kout(k_pn,k,k4,k1,k_point)
          phc(k1,k2,k4,o1,o2,o3,o4)= phc(k1,k2,k4,o1,o2,o3,o4)+vert0(k1,kindex,k4,o1,oo4,o4,oo1)&
                                     &*lphc*vert0(k,k2,kindex,oo3,o2,oo2,o3)
                enddo
            enddo
        enddo
    enddo
enddo

enddo
enddo
enddo

enddo
enddo
enddo
enddo

! the variation of V 
lambda=lambda-dlambda
vert=vert0-(-pp+phd-phc)*dlambda
print*,maxval(real(vert)),minval(real(vert))
call cpu_time(finish)
print*,"total time:",finish-start,"s"

! if(maxval(abs(vert))<vertmax) then
!     vert0=vert
!     step=step+1
!     print*,"step",step,lambda,maxval(real(vert)),minval(real(vert))
! goto 100
! else
! print*,"divergency:",maxval(real(vert)),minval(real(vert))
! endif

end
!=====================the end main===============================!

!-----------------------------------------------------------!
! input: k_pn discresization of BZ
! output: the matrix of k_point for coarse discretization in bz
subroutine kmesh(k_pn,k_point)
    implicit none
    integer::i,j,k,k_pn
    real,parameter::pi=3.1415926
    real,dimension(k_pn**3,3)::k_point
       
    open(10,file="kmesh.txt")
   
    do i=0,k_pn-1
      do j=0,k_pn-1
        do k=0,k_pn-1
    
         k_point(i*k_pn**2+j*k_pn+k+1,3)=(-1.0+2.0*k/(k_pn-1))*pi
         k_point(i*k_pn**2+j*k_pn+k+1,2)=(-1.0+2.0*j/(k_pn-1))*pi
         k_point(i*k_pn**2+j*k_pn+k+1,1)=(-1.0+2.0*i/(k_pn-1))*pi
       
            enddo
        enddo
    enddo

    do i=1,k_pn**3
        write(10,*)i,k_point(i,1:3)
    enddo

    close(10) 
end
!---------------------the end---------------------------------------!

!----------calculate energy and its eigenvectors----------------!
subroutine energy_cal(k_pn,kindex,k_point,h,enerval)
 implicit none
 Integer, Parameter :: LWORK = 5000 !LWORK>=MAX(2*2*N-1)
 Integer :: INFO
 real, Parameter :: pi = atan(1.0)*4
 complex,dimension(2,2)::h
 real, Dimension (2) :: enerval
 Complex, Dimension (LWORK) :: WORK
 real, Dimension (3*2-2) :: RWORK
 Complex, Dimension(2,2) :: sigmax,sigmay,sigmaz,sigma0
 real,dimension(k_pn**3,3)::k_point
 integer::kindex,k_pn
 real::gamma,m0,kx0,tx

 gamma=0.0  ! This parameter can be tunned to change typeI to typeII.
 m0=2.0
 tx=1.0
 kx0=pi/2.0
 sigma0=0.0
 sigmax=0.0
 sigmay=0.0
 sigmaz=0.0 
 sigma0(1,1)=1.0
 sigma0(2,2)=1.0
 sigmax(1,2)=1.0
 sigmax(2,1)=1.0
 sigmaz(1,1)=1.0
 sigmaz(2,2)=-1.0
 sigmay(1,2)=cmplx(0.0,-1.0)
 sigmay(2,1)=cmplx(0.0,1.0)    

  h=gamma*(cos(k_point(kindex,1))-cos(kx0))*sigma0-(m0*(2.0-cos(k_point(kindex,2))&
  &-cos(k_point(kindex,3)))+2.0*tx*(cos(k_point(kindex,1))-cos(kx0)))*sigmax-2.0*sin(k_point(kindex,2))*sigmay&
  &-2.0*sin(k_point(kindex,3))*sigmaz
 
  Call cheev ('V', 'U', 2, h, 2, enerval, WORK, LWORK, RWORK, &
     & INFO)
 end
!---------------------the end---------------------------------------!

!----------calculate the every element of Lpp matrix----------------!
!labled by Lpp(k1,k2,k,o1,o2,o3,o4)
subroutine lpp_cal(k_pn,k_point,k1,k2,k,o1,o2,o3,o4,lambda,lpp)
    implicit none
    !-------------------------------------------------------------------!
    ! input indexs are the lables of k point k1+k2-k
    !
    integer::k1,k2,k,kindex,o1,o2,o3,o4
    integer::k_pn,i,j
    integer::kout
    real::lambda,delta
    real,dimension(k_pn**3,3)::k_point
    complex,dimension(2,2)::hk,hkindex
    real,dimension(2)::valk,valkindex
    complex::lpp,lppfre ! the final result for every lpp(k1,k2,k,o1,o2,o3,o4)

    lpp=0.0
    delta=0.0000001
    call energy_cal(k_pn,k,k_point,hk,valk)
    !print*,valk(:)
    kindex=kout(k_pn,k1,k2,k,k_point)
    ! print*,"k1,k2,k,kindex",k1,k2,k,kindex
    call energy_cal(k_pn,kindex,k_point,hkindex,valkindex)
    !print*,valkindex(:)
    !------------pp bubble in the band basis------------------
    !lpp(k(i),k1+k2-k(kindex(j)),o1,o2,o3,o4)
    
    do i=1,2
        do j=1,2
    !print*,valk(i),valkindex(j)
    if(valk(i)>0.0.and.valkindex(j)<=0.0) then
        !print*,"case 1"
        lppfre=(valk(i)**2*valkindex(j)+lambda**3)*(3.0*valkindex(j)-lambda)+&
               &valk(i)*lambda*(valkindex(j)**2+6*valkindex(j)*lambda-3*lambda**2)
        lppfre=lppfre/((valk(i)+lambda)**3)
        lppfre=lppfre/((valkindex(j)-lambda)**3)

        else if(valk(i)>0.0.and.valkindex(j)>0.0) then
            !print*,"case 2"
        lppfre=valk(i)*(valk(i)*3.0+lambda)/((valk(i)+lambda)**3)+&
               &valkindex(j)*(valkindex(j)*3.0+lambda)/((valkindex(j)+lambda)**3)
        lppfre=lppfre/(valk(i)+valkindex(j)+cmplx(0.0,delta))

        else if(valk(i)<=0.0.and.valkindex(j)>0.0)then
            !print*,"case 3"
        lppfre=(valk(i)**2*valkindex(j)-lambda**3)*(3.0*valkindex(j)-lambda)+&
               &valk(i)*lambda*(-valkindex(j)**2+6*valkindex(j)*lambda+3*lambda**2)
        lppfre=lppfre/((valk(i)-lambda)**3)
        lppfre=lppfre/((valkindex(j)+lambda)**3)
        else
            !print*,"case 4"
        lppfre=valk(i)*(valk(i)*3.0-lambda)/((valk(i)-lambda)**3)+&
               &valkindex(j)*(valkindex(j)*3.0-lambda)/((valkindex(j)-lambda)**3)
        lppfre=lppfre/(valk(i)+valkindex(j)+cmplx(0.0,delta))
        endif

        lppfre=-0.25*lppfre
    
        lpp=lpp+lppfre*hk(o3,i)*conjg(hk(o1,i))*hkindex(o4,j)*conjg(hkindex(o2,j))
       enddo
       enddo
    
end
!---------------------the end---------------------------------------!

!----------calculate the every element of Lph matrix----------------!
!labled by Lph(k1,k3,k,o1,o2,o3,o4)
subroutine lph_cal(k_pn,k_point,k1,k3,k,o1,o2,o3,o4,lambda,lph)
    implicit none
    !-------------------------------------------------------------------!
    ! input indexs are the lables of k point k+k3-k1
    !
    integer::k1,k3,k,kindex,o1,o2,o3,o4
    integer::k_pn,i,j
    integer::kout
    real::lambda,delta
    real,dimension(k_pn**3,3)::k_point
    complex,dimension(2,2)::hk,hkindex
    real,dimension(2)::valk,valkindex
    complex::lph,lphfre ! the final result for every lpp(k1,k2,k,o1,o2,o3,o4)

    lph=0.0
    delta=0.0000001
    call energy_cal(k_pn,k,k_point,hk,valk)
    
    kindex=kout(k_pn,k,k3,k1,k_point)

    call energy_cal(k_pn,kindex,k_point,hkindex,valkindex)

    !print*,valkindex(:)
    !------------ph bubble in the band basis------------------
    !lph(k(i),k+k3-k1(kindex(j)),o1,o2,o3,o4)
    
    do i=1,2
        do j=1,2
    !print*,valk(i),valkindex(j)
    if(valk(i)>0.0.and.valkindex(j)<=0.0) then
       
        lphfre=valk(i)*(valk(i)*3.0+lambda)/((valk(i)+lambda)**3)-&
             &valkindex(j)*(valkindex(j)*3.0-lambda)/((valkindex(j)-lambda)**3)
        lphfre=lphfre/(valk(i)-valkindex(j))
 
        else if(valk(i)>0.0.and.valkindex(j)>0.0) then
            lphfre=(-valk(i)**2*valkindex(j)+lambda**3)*(3.0*valkindex(j)+lambda)+&
                    &valk(i)*lambda*(-valkindex(j)**2+6*valkindex(j)*lambda+3*lambda**2)
            lphfre=lphfre/((valk(i)+lambda)**3)
            lphfre=lphfre/((valkindex(j)+lambda)**3)

        else if(valk(i)<=0.0.and.valkindex(j)>0.0)then
            
            lphfre=valk(i)*(valk(i)*3.0-lambda)/((valk(i)-lambda)**3)-&
                  &valkindex(j)*(valkindex(j)*3.0+lambda)/((valkindex(j)+lambda)**3)
            lphfre=lphfre/(valk(i)-valkindex(j))
        else
            
            lphfre=(-valk(i)**2*valkindex(j)-lambda**3)*(3.0*valkindex(j)-lambda)+&
                  &valk(i)*lambda*(valkindex(j)**2+6*valkindex(j)*lambda-3*lambda**2)
            lphfre=lphfre/((valk(i)-lambda)**3)
            lphfre=lphfre/((valkindex(j)-lambda)**3)
        endif

        lphfre=0.25*lphfre
    
        lph=lph+lphfre*hk(o3,i)*conjg(hk(o1,i))*hkindex(o4,j)*conjg(hkindex(o2,j))
       enddo
       enddo
    
 end
!---------------------the end---------------------------------------!

!-------------------the find the index of k1+k2-k-------------------!
function kout(k_pn,k1,k2,k,k_point)
    implicit none
    !calculate the index of k1+k2-k
    !
    integer::i,k_pn,k1,k2,k,kout
    real, Parameter :: pi = atan(1.0)*4
    real,dimension(k_pn**3,3)::k_point
    real,dimension(3)::kadd
    integer,dimension(3)::ik

    kout=0
    kadd(:)=k_point(k1,:)+k_point(k2,:)-k_point(k,:)
    
    do i=1,3
    if(anint(kadd(i)/pi)<-1) then
        kadd(i)=kadd(i)+2.0*pi
    else if(anint(kadd(i)/pi)>1) then
        kadd(i)=kadd(i)-2.0*pi
    else
        kadd(i)=kadd(i)
    end if
    enddo

    do i=1,3
        ik(i)=anint((kadd(i)/pi+1)*(k_pn-1)/2.0)
    enddo

    do i=1,3
    kout=kout+ik(i)*k_pn**(3-i)
    enddo

    kout=kout+1
    !write(11,*)k1,k2,k,kout
   
end function
