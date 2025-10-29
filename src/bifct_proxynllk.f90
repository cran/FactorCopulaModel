! PROGRAM test   ! main function
! implicit none
! 	integer npar,mgrp,dvar,n,i
! 	double precision,dimension(:,:),allocatable::udata,vlat
! 	double precision,dimension(:),allocatable::th
! 	integer,dimension(:),allocatable::grsize
! 	double precision nllk
!     read *,npar
!     read *,dvar
!     read *,mgrp
!     allocate(grsize(mgrp))
!     read *,grsize
!     read*,n
!     allocate(th(npar),udata(n,dvar),vlat(n,mgrp+1))
!  	read *,th
!  	read *, (udata(i,:), i=1,n) 
!  	! do i=1,n
!     !      print *, udata(i,:)
!     !    end do
!  	read *, (vlat(i,:), i=1,n) 
!  	! do i=1,n
!      !     print *, vlat(i,:)
!      !   end do
!     
!  	call frkproxynllk(npar,th,mgrp,n,dvar,grsize,udata,vlat,nllk)
!  	
!  	print *, "the negative log-likelihood is"
!  	print *,nllk
!  	deallocate(grsize,udata,vlat,th)
! END PROGRAM


subroutine frkproxynllk(npar,th,mgrp,n,dvar,grsize,udata,vlat,nllk)
  implicit none
  integer npar,mgrp,dvar,n
  double precision th(npar),udata(n,dvar),uvec(dvar),vlat(n,mgrp+1)
  integer grsize(mgrp)
  double precision nllk,lk1,lk2,liki,ccdf(dvar),pdf(npar)
  integer i,jg,mj,ind1,ind2,ind,int0
  double precision intj(mgrp)

   nllk=0.d0; 
   do i =1,n 
       uvec = udata(i,:) 
       int0 = 1.d0; ind = 0.d0; intj = 1.d0;  
      do jg =1,mgrp   !jth group                         
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);lk1 = 1.d0; lk2 = 1.d0;  
          do mj = ind1,ind2 ! within group index
            call pcondfrk(uvec(mj),vlat(i,1),th(mj),ccdf(mj)) !C_{ij|V_0}
            call dfrk(uvec(mj),vlat(i,1),th(mj),pdf(mj)) !c_{ij,V_0} 
            call dfrk(ccdf(mj),vlat(i,jg+1),th(dvar+mj),pdf(dvar+mj)) !c_{i,V_j;V_0}  
    
            lk1 = lk1*pdf(mj)
            lk2 = lk2*pdf(dvar+mj)
          end do
          intj(jg) = intj(jg)*lk1*lk2 !intj value of the jth integral
        end do   
      liki = product(intj)   !product of j inner integrals  
      nllk = nllk - log(liki)       !updating loglikelihood
  end do
  return
  end

subroutine gumproxynllk(npar,th,mgrp,n,dvar,grsize,udata,vlat,nllk)
  implicit none
  integer npar,mgrp,dvar,n
  double precision th(npar),udata(n,dvar),uvec(dvar),vlat(n,mgrp+1)
  integer grsize(mgrp)
  double precision nllk,lk1,lk2,liki,ccdf(dvar),pdf(npar)
  integer i,jg,mj,ind1,ind2,ind,int0
  double precision intj(mgrp)

  nllk=0.d0; 
  do i =1,n 
       uvec = udata(i,:) 
       int0 = 1.d0; ind = 0.d0; intj = 1.d0;  
      do jg =1,mgrp   !jth group                         
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);lk1 = 1.d0; lk2 = 1.d0;  
          do mj = ind1,ind2 ! within group index
            call pcondgum(uvec(mj),vlat(i,1),th(mj),ccdf(mj)) !C_{ij|V_0}
            call dgum(uvec(mj),vlat(i,1),th(mj),pdf(mj)) !c_{ij,V_0} 
            call dgum(ccdf(mj),vlat(i,jg+1),th(dvar+mj),pdf(dvar+mj)) !c_{i,V_j;V_0}  
            lk1 = lk1*pdf(mj)
            lk2 = lk2*pdf(dvar+mj)
          end do
          intj(jg) = intj(jg)*lk1*lk2 !intj value of the jth integral
        end do   
      liki = product(intj)   !product of j inner integrals  
      nllk = nllk - log(liki)       !updating loglikelihood
  end do
  return
  end
  



