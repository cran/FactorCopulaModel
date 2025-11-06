
!  main program
! PROGRAM test3   ! main function
!   implicit none
!   integer npar,mgrp,dvar,n,i
!   integer,dimension(:),allocatable::family
!   double precision,dimension(:,:),allocatable::udata,vlat,hess
!   double precision,dimension(:),allocatable::th,grad
!   integer,dimension(:),allocatable::grsize
!   double precision nllk
  
!      read *,npar

!      print*,npar

!      read *,dvar

!      print*,dvar
!      read *,mgrp
!      print*,mgrp

!     allocate(grsize(mgrp))
!     allocate(family(2*dvar))

!     read *,grsize
!     print*,grsize

!     read*,family
!     print*,family

!     read*,n
!     print*,n
!     allocate(th(npar),udata(n,dvar),vlat(n,mgrp+1))
!     allocate(grad(npar),hess(npar,npar))
!     read *,th
!     print*,th

!     read *, (udata(i,:), i=1,n) 

!      do i=1,n
!         print *, udata(i,:)
!       end do

!     read *, (vlat(i,:), i=1,n) 

!     do i=1,n
!            print *, vlat(i,:)
!       end do
    
!     call bifactproxy(npar,th,mgrp,n,dvar,family,grsize,udata,vlat,nllk,grad,hess)
!    ! call strbb1frk2proxy(npar,th,mgrp,n,dvar,grsize,udata,vlat,nllk,grad,hess)
  
!    print *, "the negative log-likelihood is"
!    print *,nllk

!     print*,"the gradient is"
!     print*, grad

!    ! print*,"the hessian matrix is"
!    ! print*,hess
!   deallocate(grsize,family,udata,vlat,th,grad,hess)
! END PROGRAM



! subroutine isinfinite(number,res)
!   USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE      
!   DOUBLE PRECISION :: number
!   integer res

!   IF (IEEE_IS_FINITE(number)) THEN
!      res=0
!   ELSE
!      res=1
!   END IF         
!   end
subroutine strbb1frk2proxy(npar,th,mgrp,n,dvar,grsize,udata,vlat,nllk,grad,hess)
  implicit none
  integer npar,mgrp,dvar,n,ip,jp,ind,m1,m2,flag
  double precision tem
  double precision th(npar),udata(n,dvar),uvec(dvar),vlat(n,mgrp+1)
  double precision lderu(4),lderv(4),lderuv
  double precision lpdf(2*dvar),der1(npar),der2(dvar+npar)
  double precision ccdf(dvar),cder1(2*dvar),cder2(npar),der1bb1(2),der2bb1(3),cder1bb1(2),cder2bb1(3)
  double precision nllk,liki,llk,lk,lk2,grad(npar),hess(npar,npar)
  integer i,jg,mj,mj2,ind1,ind2,grsize(mgrp),ibb1(2)
  double precision intj(mgrp), intjj(npar), grdj(npar), der2j(npar), int0, grd0(npar), grdi(npar)
  double precision der1u(dvar), der2u(dvar),derumix(dvar), dermxj(npar)
  double precision hss0(npar,npar), hssi(npar,npar), hessi(npar,npar), hssaux(npar,npar)

  ! npar = 3*dvar; dvar = sum(grsize)
  nllk=0.d0; grad = 0.d0; hess=0.d0 
  do i =1,n 
    uvec = udata(i,:)    
    liki = 1.d0; grdi = 1.d0; hssi = 1.d0; 
    !i0 ---> iq
    !do iq =1,nq
      int0 = 1.d0; grd0 = 1.d0; hssaux = 0.d0 
      ind = 0; intj = 1.d0;  grdj = 0.d0; der2j = 0.d0; dermxj = 0.d0 
      !j ---> jg
      do jg =1,mgrp   !jth group                         
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
        !ij ---> iq2
        !do iq2 =1,nq   
          lk = 1.d0; lk2 = 1.d0 
          do mj = ind1,ind2
            !call cbb1derivs(uvec(mj),vlat(i,1),th((/mj,mj+dvar/)),ccdf(mj),cder11,cder22)
            call cbb1derivs(uvec(mj),vlat(i,1),th((/mj,mj+dvar/)),ccdf(mj),cder1bb1,cder2bb1)
            !C_{ij|V_0}
            !call cbb1derivs(uvec(mj),vlat(i,1),th((/mj,mj+dvar/)),ccdf(mj),cder1bb1,cder2bb1)                                          
            cder1((/mj,mj+dvar/)) = cder1bb1
            cder2((/mj,mj+2*dvar,mj+dvar/)) = cder2bb1
            !c_{i,V_j;V_0} 
            call lfrk2derivs(ccdf(mj),vlat(i,jg+1),th(2*dvar+mj),lpdf(dvar+mj),der1(2*dvar+mj),der2(3*dvar+mj),&
                 der1u(mj),der2u(mj),derumix(mj))   
            !c_{ij,V_0}
            call lbb1derivs2(uvec(mj),vlat(i,1),th((/mj,mj+dvar/)),lpdf(mj),der1bb1,der2bb1,lderu,lderv,lderuv)
            !call lbb1derivs(uvec(mj),vlat(i,1),th((/mj,mj+dvar/)),lpdf(mj),der1bb1,der2bb1)                                           
            der1((/mj,mj+dvar/)) = der1bb1
            der2((/mj,mj+2*dvar,mj+dvar/)) = der2bb1  
            lk = lk*exp(lpdf(mj))
            lk2 = lk2*exp(lpdf(dvar+mj))
          end do
          llk = lk*lk2
          intj(jg) = intj(jg)*llk !intj value of the jth integral
          intjj(ind1:ind2) = intj(jg)  !intjj: repeat intj(j) grsize[j] times
          intjj((dvar+ind1):(dvar+ind2)) = intj(jg)
          intjj((2*dvar+ind1):(2*dvar+ind2)) = intj(jg)
          do mj = ind1,ind2
            ibb1 = (/mj,mj+dvar/) 
            !1st o.d. wrt th.bb1,dl.bb1
            grdj(ibb1) = grdj(ibb1) + llk*(der1u(mj)*cder1(ibb1)+der1(ibb1))                                        
            !1st o.d. wrt th.frk
            grdj(2*dvar+mj) = grdj(2*dvar+mj)+llk*der1(2*dvar+mj)
            !2nd o.d. wrt th.bb1,dl.bb1
            der2j(ibb1) = der2j(ibb1) + llk*(der2u(mj)+der1u(mj)*der1u(mj))*cder1(ibb1)*cder1(ibb1)                  
            der2j(ibb1) = der2j(ibb1) + llk*(der1u(mj)*cder2(ibb1)+der2(ibb1)+der1(ibb1)*der1(ibb1))       
            der2j(ibb1) = der2j(ibb1) + llk*2*der1(ibb1)*der1u(mj)*cder1(ibb1)                            
            !2nd o.d. wrt th.frk
            der2j(2*dvar+mj) = der2j(2*dvar+mj)+ llk*(der2(3*dvar+mj)+der1(2*dvar+mj)*der1(2*dvar+mj)) 
            !mixed d. wrt th.bb1 and dl.bb1
            dermxj(ibb1) = dermxj(ibb1) + llk &
                   *(der1(ibb1)*der1(2*dvar+mj)+(derumix(mj)+der1u(mj)*der1(2*dvar+mj))*cder1(ibb1))  
            dermxj(2*dvar+mj) = dermxj(2*dvar+mj) + llk &
                   *(der2(mj+2*dvar)+der1(mj)*der1(mj+dvar)+der1u(mj)*der1(dvar+mj)*cder1(mj))
            dermxj(2*dvar+mj) = dermxj(2*dvar+mj) + llk &
                   *(der1u(mj)*der1(mj)*cder1(mj+dvar)+der1u(mj)*cder2(2*dvar+mj))
            dermxj(2*dvar+mj) = dermxj(2*dvar+mj) + llk &
                   *((der2u(mj)+der1u(mj)*der1u(mj))*cder1(mj)*cder1(dvar+mj))
            do mj2 = ind1,(mj-1)  !2nd order der. of the jth int. wrt th(mj), th(mj2)
              hssaux(mj,mj2)  = hssaux(mj,mj2) + llk*(der1u(mj)*cder1(mj)+der1(mj))*(der1u(mj2)*cder1(mj2)+der1(mj2))
              hssaux(mj+dvar,mj2) = hssaux(mj+dvar,mj2) + llk &
                    *(der1u(mj)*cder1(mj+dvar)+der1(mj+dvar))*(der1u(mj2)*cder1(mj2)+der1(mj2))
              hssaux(mj,mj2+dvar) = hssaux(mj,mj2+dvar) + llk &
                    *(der1u(mj)*cder1(mj)+der1(mj))*(der1u(mj2)*cder1(mj2+dvar) &
                    +der1(mj2+dvar))
              hssaux(mj+dvar,mj2+dvar) = hssaux(mj+dvar,mj2+dvar) + llk &
                    *(der1u(mj)*cder1(mj+dvar)+der1(mj+dvar))*(der1u(mj2)*cder1(mj2+dvar)+der1(mj2+dvar))
              hssaux(mj+2*dvar,mj2) = hssaux(mj+2*dvar,mj2) + llk*der1(2*dvar+mj)*(der1u(mj2)*cder1(mj2)+der1(mj2))
              hssaux(mj2+2*dvar,mj) = hssaux(mj2+2*dvar,mj) + llk*der1(2*dvar+mj2)*(der1u(mj)*cder1(mj)+der1(mj))
              hssaux(mj+2*dvar,mj2+dvar) = hssaux(mj+2*dvar,mj2+dvar) + llk &
                    *der1(2*dvar+mj)*(der1u(mj2)*cder1(mj2+dvar)+der1(mj2+dvar))
              hssaux(mj2+2*dvar,mj+dvar) = hssaux(mj2+2*dvar,mj+dvar) + llk &
                    *der1(2*dvar+mj2)*(der1u(mj)*cder1(mj+dvar)+der1(mj+dvar))
              hssaux(mj+2*dvar,mj2+2*dvar) = hssaux(mj+2*dvar,mj2+2*dvar) + llk*der1(2*dvar+mj)*der1(2*dvar+mj2)
              hssaux(mj2,mj)  = hssaux(mj,mj2)
              hssaux(mj2,mj+dvar) = hssaux(mj+dvar,mj2)
              hssaux(mj2+dvar,mj) = hssaux(mj,mj2+dvar) 
              hssaux(mj2+dvar,mj+dvar) = hssaux(mj+dvar,mj2+dvar)
              hssaux(mj2,mj+2*dvar) = hssaux(mj+2*dvar,mj2)
              hssaux(mj,mj2+2*dvar) = hssaux(mj2+2*dvar,mj) 
              hssaux(mj2+dvar,mj+2*dvar) = hssaux(mj+2*dvar,mj2+dvar) 
              hssaux(mj+dvar,mj2+2*dvar) = hssaux(mj2+2*dvar,mj+dvar)  
              hssaux(mj2+2*dvar,mj+2*dvar) = hssaux(mj+2*dvar,mj2+2*dvar)  
            end do                    
          end do
        !end do          
      end do

      int0 = product(intj) !product of j inner integrals  
      grd0 = int0*(grdj/intjj)!gradient of the product
      !m1 ---> ip; m2 ---> jp
      do ip=2,npar !hessian of the product
        do jp=1,ip-1
          hss0(ip,jp) = int0*grdj(ip)*grdj(jp)/intjj(ip)/intjj(jp)
          hss0(jp,ip) = hss0(ip,jp)
        end do
      end do
      ind = 0
      !j ---> jg
      do jg=1,mgrp    !2nd order derivatives within groups 
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg)
        hss0(ind1:ind2,ind1:ind2)=int0*hssaux(ind1:ind2,ind1:ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,dvar+ind1:dvar+ind2)=int0*hssaux(dvar+ind1:dvar+ind2,dvar+ind1:dvar+ind2)/intj(jg)
        hss0(2*dvar+ind1:2*dvar+ind2,2*dvar+ind1:2*dvar+ind2)=int0*hssaux(2*dvar+ind1:2*dvar+ind2,2*dvar+ind1:2*dvar+ind2)/intj(jg)
        hss0(ind1:ind2,dvar+ind1:dvar+ind2)=int0*hssaux(ind1:ind2,dvar+ind1:dvar+ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,ind1:ind2)=int0*hssaux(dvar+ind1:dvar+ind2,ind1:ind2)/intj(jg)
        hss0(2*dvar+ind1:2*dvar+ind2,ind1:ind2)=int0*hssaux(2*dvar+ind1:2*dvar+ind2,ind1:ind2)/intj(jg)
        hss0(ind1:ind2,2*dvar+ind1:2*dvar+ind2)=int0*hssaux(ind1:ind2,2*dvar+ind1:2*dvar+ind2)/intj(jg)
        hss0(2*dvar+ind1:2*dvar+ind2,dvar+ind1:dvar+ind2)=int0*hssaux(2*dvar+ind1:2*dvar+ind2,dvar+ind1:dvar+ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,2*dvar+ind1:2*dvar+ind2)=int0*hssaux(dvar+ind1:dvar+ind2,2*dvar+ind1:2*dvar+ind2)/intj(jg)
      end do

      do ip=1,npar
        hss0(ip,ip)  = int0*der2j(ip)/intjj(ip)
      end do
      do ip=1,dvar
        hss0(ip,ip+dvar)  = int0*dermxj(ip+2*dvar)/intjj(ip)
        hss0(ip+2*dvar,ip)  = int0*dermxj(ip)/intjj(ip)
        hss0(ip+2*dvar,ip+dvar) = int0*dermxj(ip+dvar)/intjj(ip)
        hss0(ip+dvar,ip)  = hss0(ip,ip+dvar)
        hss0(ip,ip+2*dvar)  = hss0(ip+2*dvar,ip)
        hss0(ip+dvar,ip+2*dvar) = hss0(ip+2*dvar,ip+dvar)
      end do
      liki = liki*int0   !integrating over V0
      grdi = grdi*grd0
      hssi = hssi*hss0
    !end do
    ! nllk = nllk - log(liki)       !updating loglikelihood
    ! grad = grad - grdi/liki       !updating gradient
    do ip=1,npar                  !updating hessian
      do jp=1,ip
        hessi(ip,jp)=hssi(ip,jp)/liki-grdi(ip)*grdi(jp)/liki/liki 
        hessi(jp,ip)=hessi(ip,jp)
      end do
    end do 

    !add some code to handle the infinite values and nan values
    flag=0;
    do m1=1,npar
      do m2=1,m1
        tem=hessi(m1,m2)
        ! call isinfinite(tem,res)
        if(isnan(tem).or.(tem.eq.(tem-1.0d0))) then !nan values
          flag=1;
          ! print*,tem
          ! print*, "there is inifnite values"
          exit;
        else
          flag=0;
        end if 
       end do 
       
      if(flag==1) then
        exit;
      end if  
    end do
    ! print*,i
    ! print*,flag
    if(flag==1) then
      nllk=nllk
      grad=grad
      hess=hess
    else
      nllk = nllk - log(liki)       !updating loglikelihood
      grad = grad - grdi/liki       !updating gradient
      hess=hess-hessi
    end if 
    
  end do
  return
  end


subroutine bifactproxy(npar,th,mgrp,n,dvar,family,grsize,udata,vlat,nllk,grad,hess)
  implicit none
  integer npar,mgrp,dvar,n,ip,jp,ind
  integer flag,m1,m2
  double precision tem
  integer family(2*dvar)
  double precision th(npar),udata(n,dvar),uvec(dvar),vlat(n,mgrp+1)
  double precision llpdf(2*dvar),der1(npar),der2(dvar+npar)
  double precision lder2v,lder2uv,lder1v
  double precision ccdf(dvar),cder1(2*dvar),cder2(npar),ldermixuu(2),ldermixvv(2)
  double precision lder1u,lder2u,lder11(2),lder22(3),cder1vec(2),cder2vec(3)
  double precision nllk,liki,llk,lk,lk2,grad(npar),hess(npar,npar)
  integer i,jg,mj,mj2,ind1,ind2,grsize(mgrp),ibb1(2)
  double precision intj(mgrp), intjj(npar), grdj(npar), der2j(npar), int0, grd0(npar), grdi(npar)
  double precision der1u(dvar), der2u(dvar),derumix(dvar), dermxj(npar)
  double precision hss0(npar,npar), hssi(npar,npar), hessi(npar,npar), hssaux(npar,npar)
 

  ! npar = 3*dvar; dvar = sum(grsize)
  nllk=0.d0; grad = 0.d0; hess=0.d0 
  do i =1,n 
    uvec = udata(i,:)    
    liki = 1.d0; grdi = 1.d0; hssi = 1.d0; 
    !i0 ---> iq
    !do iq =1,nq
      int0 = 1.d0; grd0 = 1.d0; hssaux = 0.d0 
      ind = 0; intj = 1.d0;  grdj = 0.d0; der2j = 0.d0; dermxj = 0.d0 
      !j ---> jg
      do jg =1,mgrp   !jth group                         
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
        !ij ---> iq2
        !do iq2 =1,nq   
          lk = 1.d0; lk2 = 1.d0 
          do mj = ind1,ind2
            call ccopderiv(uvec(mj),vlat(i,1),th((/mj,mj+dvar/)),family(mj),&
              ccdf(mj),cder1vec,cder2vec)
            cder1((/mj,mj+dvar/)) = cder1vec
            cder2((/mj,mj+2*dvar,mj+dvar/)) = cder2vec
            
            call lfrk2derivs(ccdf(mj),vlat(i,jg+1),th(2*dvar+mj),&
                llpdf(dvar+mj),der1(2*dvar+mj),der2(3*dvar+mj),&
                der1u(mj),der2u(mj),derumix(mj))   

            call  lcop2derivt(uvec(mj),vlat(i,1),th((/mj,mj+dvar/)),family(mj),&
              llpdf(mj),lder11,lder22,lder1u,lder2u,&
                  ldermixuu,lder1v,lder2v,lder2uv,ldermixvv)

              der1((/mj,mj+dvar/)) = lder11
              der2((/mj,mj+2*dvar,mj+dvar/)) = lder22

        
            lk = lk*exp(llpdf(mj))
            lk2 = lk2*exp(llpdf(dvar+mj))
           end do

          llk = lk*lk2

          intj(jg) = intj(jg)*llk !intj value of the jth integral
          intjj(ind1:ind2) = intj(jg)  !intjj: repeat intj(j) grsize[j] times
          intjj((dvar+ind1):(dvar+ind2)) = intj(jg)
          intjj((2*dvar+ind1):(2*dvar+ind2)) = intj(jg)
          
          do mj = ind1,ind2
            ibb1 = (/mj,mj+dvar/) 
            !1st o.d. wrt th.bb1,dl.bb1
            grdj(ibb1) = grdj(ibb1) + llk*(der1u(mj)*cder1(ibb1)+der1(ibb1))                                        
            !1st o.d. wrt th.frk
            grdj(2*dvar+mj) = grdj(2*dvar+mj)+llk*der1(2*dvar+mj)
            !2nd o.d. wrt th.bb1,dl.bb1
            der2j(ibb1) = der2j(ibb1) + llk*(der2u(mj)+der1u(mj)*der1u(mj))*cder1(ibb1)*cder1(ibb1)                  
            der2j(ibb1) = der2j(ibb1) + llk*(der1u(mj)*cder2(ibb1)+der2(ibb1)+der1(ibb1)*der1(ibb1))       
            der2j(ibb1) = der2j(ibb1) + llk*2*der1(ibb1)*der1u(mj)*cder1(ibb1)                            
            !2nd o.d. wrt th.frk
            der2j(2*dvar+mj) = der2j(2*dvar+mj)+ llk*(der2(3*dvar+mj)+der1(2*dvar+mj)*der1(2*dvar+mj)) 
            !mixed d. wrt th.bb1 and dl.bb1
            dermxj(ibb1) = dermxj(ibb1) + llk &
                   *(der1(ibb1)*der1(2*dvar+mj)+(derumix(mj)+der1u(mj)*der1(2*dvar+mj))*cder1(ibb1))  

            dermxj(2*dvar+mj) = dermxj(2*dvar+mj) + llk &
                   *(der2(mj+2*dvar)+der1(mj)*der1(mj+dvar)+der1u(mj)*der1(dvar+mj)*cder1(mj))


            dermxj(2*dvar+mj) = dermxj(2*dvar+mj) + llk &
                   *(der1u(mj)*der1(mj)*cder1(mj+dvar)+der1u(mj)*cder2(2*dvar+mj))
                   
            dermxj(2*dvar+mj) = dermxj(2*dvar+mj) + llk &
                   *((der2u(mj)+der1u(mj)*der1u(mj))*cder1(mj)*cder1(dvar+mj))
            do mj2 = ind1,(mj-1)  !2nd order der. of the jth int. wrt th(mj), th(mj2)
              hssaux(mj,mj2)  = hssaux(mj,mj2) + llk*(der1u(mj)*cder1(mj)+der1(mj))*(der1u(mj2)*cder1(mj2)+der1(mj2))
              hssaux(mj+dvar,mj2) = hssaux(mj+dvar,mj2) + llk &
                    *(der1u(mj)*cder1(mj+dvar)+der1(mj+dvar))*(der1u(mj2)*cder1(mj2)+der1(mj2))
              hssaux(mj,mj2+dvar) = hssaux(mj,mj2+dvar) + llk &
                    *(der1u(mj)*cder1(mj)+der1(mj))*(der1u(mj2)*cder1(mj2+dvar) &
                    +der1(mj2+dvar))
              hssaux(mj+dvar,mj2+dvar) = hssaux(mj+dvar,mj2+dvar) + llk &
                    *(der1u(mj)*cder1(mj+dvar)+der1(mj+dvar))*(der1u(mj2)*cder1(mj2+dvar)+der1(mj2+dvar))


              hssaux(mj+2*dvar,mj2) = hssaux(mj+2*dvar,mj2) + llk*der1(2*dvar+mj)*(der1u(mj2)*cder1(mj2)+der1(mj2))
              hssaux(mj2+2*dvar,mj) = hssaux(mj2+2*dvar,mj) + llk*der1(2*dvar+mj2)*(der1u(mj)*cder1(mj)+der1(mj))


              hssaux(mj+2*dvar,mj2+dvar) = hssaux(mj+2*dvar,mj2+dvar) + llk &
                    *der1(2*dvar+mj)*(der1u(mj2)*cder1(mj2+dvar)+der1(mj2+dvar))
              hssaux(mj2+2*dvar,mj+dvar) = hssaux(mj2+2*dvar,mj+dvar) + llk &
                    *der1(2*dvar+mj2)*(der1u(mj)*cder1(mj+dvar)+der1(mj+dvar))
              hssaux(mj+2*dvar,mj2+2*dvar) = hssaux(mj+2*dvar,mj2+2*dvar) + llk*der1(2*dvar+mj)*der1(2*dvar+mj2)



              hssaux(mj2,mj)  = hssaux(mj,mj2)
              hssaux(mj2,mj+dvar) = hssaux(mj+dvar,mj2)
              hssaux(mj2+dvar,mj) = hssaux(mj,mj2+dvar) 
              hssaux(mj2+dvar,mj+dvar) = hssaux(mj+dvar,mj2+dvar)


              hssaux(mj2,mj+2*dvar) = hssaux(mj+2*dvar,mj2)
              hssaux(mj,mj2+2*dvar) = hssaux(mj2+2*dvar,mj) 
              hssaux(mj2+dvar,mj+2*dvar) = hssaux(mj+2*dvar,mj2+dvar) 
              hssaux(mj+dvar,mj2+2*dvar) = hssaux(mj2+2*dvar,mj+dvar)  
              hssaux(mj2+2*dvar,mj+2*dvar) = hssaux(mj+2*dvar,mj2+2*dvar)  
            end do                    
          end do
          
        
        !end do          
      end do


      int0 = product(intj) !product of j inner integrals  
      grd0 = int0*(grdj/intjj)!gradient of the product

      do ip=2,npar !hessian of the product
        do jp=1,ip-1
          hss0(ip,jp) = int0*grdj(ip)*grdj(jp)/intjj(ip)/intjj(jp)
          hss0(jp,ip) = hss0(ip,jp)
        end do
      end do
      ind = 0
      !j ---> jg
      do jg=1,mgrp    !2nd order derivatives within groups 
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg)
        hss0(ind1:ind2,ind1:ind2)=int0*hssaux(ind1:ind2,ind1:ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,dvar+ind1:dvar+ind2)=int0*hssaux(dvar+ind1:dvar+ind2,dvar+ind1:dvar+ind2)/intj(jg)
        hss0(2*dvar+ind1:2*dvar+ind2,2*dvar+ind1:2*dvar+ind2)=int0*hssaux(2*dvar+ind1:2*dvar+ind2,2*dvar+ind1:2*dvar+ind2)/intj(jg)
        hss0(ind1:ind2,dvar+ind1:dvar+ind2)=int0*hssaux(ind1:ind2,dvar+ind1:dvar+ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,ind1:ind2)=int0*hssaux(dvar+ind1:dvar+ind2,ind1:ind2)/intj(jg)
        hss0(2*dvar+ind1:2*dvar+ind2,ind1:ind2)=int0*hssaux(2*dvar+ind1:2*dvar+ind2,ind1:ind2)/intj(jg)
        hss0(ind1:ind2,2*dvar+ind1:2*dvar+ind2)=int0*hssaux(ind1:ind2,2*dvar+ind1:2*dvar+ind2)/intj(jg)
        hss0(2*dvar+ind1:2*dvar+ind2,dvar+ind1:dvar+ind2)=int0*hssaux(2*dvar+ind1:2*dvar+ind2,dvar+ind1:dvar+ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,2*dvar+ind1:2*dvar+ind2)=int0*hssaux(dvar+ind1:dvar+ind2,2*dvar+ind1:2*dvar+ind2)/intj(jg)
      end do
      
    
      do ip=1,npar
        hss0(ip,ip)  = int0*der2j(ip)/intjj(ip)
      end do
      do ip=1,dvar
        hss0(ip,ip+dvar)  = int0*dermxj(ip+2*dvar)/intjj(ip)
        hss0(ip+2*dvar,ip)  = int0*dermxj(ip)/intjj(ip)
        hss0(ip+2*dvar,ip+dvar) = int0*dermxj(ip+dvar)/intjj(ip)
        hss0(ip+dvar,ip)  = hss0(ip,ip+dvar)
        hss0(ip,ip+2*dvar)  = hss0(ip+2*dvar,ip)
        hss0(ip+dvar,ip+2*dvar) = hss0(ip+2*dvar,ip+dvar)
      end do
      liki = liki*int0   !integrating over V0
      grdi = grdi*grd0
      hssi = hssi*hss0

     
    !end do
   
    do ip=1,npar                  !updating hessian
      do jp=1,ip
        hessi(ip,jp)=hssi(ip,jp)/liki-grdi(ip)*grdi(jp)/liki/liki 
        hessi(jp,ip)=hessi(ip,jp)
      end do
    end do 
    
    !add some code to handle the infinite values and nan values
    flag=0;
    do m1=1,npar
      do m2=1,m1
        tem=hessi(m1,m2)
        ! call isinfinite(tem,res)
        if(isnan(tem).or.(tem.eq.(tem-1.0d0))) then !nan values
          flag=1;
          ! print*,tem
          ! print*, "there is inifnite values"
          exit;
        else
          flag=0;
        end if 
       end do 
       
      if(flag==1) then
        exit;
      end if  
    end do
   

    ! print*,i
    ! print*,flag
    if(flag==1) then
      nllk=nllk
      grad=grad
      hess=hess
    else
      nllk = nllk - log(liki)       !updating loglikelihood
      grad = grad - grdi/liki       !updating gradient
      hess=hess-hessi
    end if 
    
  end do
  return
  end

 !the below code has problems
          ! if(family(mj)==4) then
          !      call cgumderivs(uvec(mj),vlat(i,1),th(mj),ccdf(mj),ccder1,ccder2)
          !      cder1vec(1)=ccder1;cder1vec(2)=0.0d0;
          !      cder2vec(1)=ccder2;cder2vec(2)=0.0d0;cder2vec(3)=0.0d0
          !      cder1((/mj,mj+dvar/)) = cder1vec
          !      cder2((/mj,mj+2*dvar,mj+dvar/)) = cder2vec
                
          !     call lfrk2derivs(ccdf(mj),vlat(i,jg+1),th(2*dvar+mj),&
          !       llpdf(dvar+mj),der1(2*dvar+mj),der2(3*dvar+mj),&
          !                       der1u(mj),der2u(mj),derumix(mj))   

          !      ! call lcop2derivt(ccdf(mj),vlat(i,jg+1),th(2*dvar+mj),family(dvar+mj), &
          !      !    llpdf(dvar+mj),lder11,lder22,der1u(mj),der2u(mj),derumix(mj),&
          !      !    lder1v,lder2v,lder2uv,ldermixvv) !c_{i,V_j;V_0}  

        
          !       ! der1(2*dvar+mj)=lder11(1)
          !       ! der2(3*dvar+mj)=lder22(1)
               
          !      call lgum2derivs(uvec(mj),vlat(i,1),th(mj),llpdf(mj),lder1,lder2,lder1u,lder2u,ldermix)
          !      der1(mj)=lder1; der1(mj+dvar)=0.0d0
          !      der2(mj)=lder2;der2(mj+2*dvar)=0.0d0;der2(mj+dvar)=0.0d0;
    
          
          !     else if (family(mj)==5) then
          !      call clfrkderivs(uvec(mj),vlat(i,1),th(mj),ccdf(mj),ccder1,ccder2) !C_{ij|V_0}
          !      cder1vec(1)=ccder1;cder1vec(2)=0.0d0;
          !      cder2vec(1)=ccder2;cder2vec(2)=0.0d0;cder2vec(3)=0.0d0
          !      cder1((/mj,mj+dvar/)) = cder1vec
          !      cder2((/mj,mj+2*dvar,mj+dvar/)) = cder2vec

          !      call lfrk2derivs(ccdf(mj),vlat(i,jg+1),th(2*dvar+mj),&
          !       llpdf(dvar+mj),der1(2*dvar+mj),der2(3*dvar+mj),&
          !                       der1u(mj),der2u(mj),derumix(mj))
             
          !       call lfrk2derivs(uvec(mj),vlat(i,1),th(mj),llpdf(mj),lder1,lder2,lder1u,lder2u,ldermix)
          !       der1(mj)=lder1; der1(mj+dvar)=0.0d0
          !       der2(mj)=lder2;der2(mj+2*dvar)=0.0d0;der2(mj+dvar)=0.0d0;
               
                 
          !     else if(family(mj)==1) then
          !        rho=th(mj);nu=5000.0d0;
          !        t1=qt(uvec(mj),nu);t2=qt(vlat(i,1),nu);
          !        call ctderivs(t1,t2,rho,nu,ccdf(mj),ccder1,ccder2)
          !         cder1vec(1)=ccder1;cder1vec(2)=0.0d0
          !         cder2vec(1)=ccder2;cder2vec(2)=0.0d0;cder2vec(3)=0.0d0
          !         cder1((/mj,mj+dvar/)) = cder1vec
          !         cder2((/mj,mj+2*dvar,mj+dvar/)) = cder2vec
                 
          !         call lfrk2derivs(ccdf(mj),vlat(i,jg+1),th(2*dvar+mj),&
          !         llpdf(dvar+mj),der1(2*dvar+mj),der2(3*dvar+mj),&
          !                       der1u(mj),der2u(mj),derumix(mj))
               
          !         call lgau2derivt(uvec(mj),vlat(i,1),th(mj),llpdf(mj),lder1,lder2,lder1u,lder2u,ldermixu,lder1v,&
          !          lder2v,lder2uv,ldermixv) !c_{ij,V_0} 
          !         der1(mj)=lder1; der1(mj+dvar)=0.0d0
          !         der2(mj)=lder2;der2(mj+2*dvar)=0.0d0;der2(mj+dvar)=0.0d0;
               

          !     else if(family(mj)==2) then
          !        call ctderivs(uvec(mj),vlat(i,1),th(mj),th(mj+dvar),ccdf(mj),ccder1,ccder2) !C_{ij|V_0}
          !        cder1vec(1)=ccder1;cder1vec(2)=0.0d0
          !        cder2vec(1)=ccder2;cder2vec(2)=0.0d0;cder2vec(3)=0.0d0
                  
          !        cder1((/mj,mj+dvar/)) = cder1vec
          !        cder2((/mj,mj+2*dvar,mj+dvar/)) = cder2vec

          !       call lfrk2derivs(ccdf(mj),vlat(i,jg+1),th(2*dvar+mj),&
          !       llpdf(dvar+mj),der1(2*dvar+mj),der2(3*dvar+mj),&
          !                       der1u(mj),der2u(mj),derumix(mj))

          !        call lt2derivs(uvec(mj),vlat(i,1),th(mj),th(mj+dvar),llpdf(mj),lder1,lder2,ltder1u,ltder2u,ltdermix) 
          !        der1(mj)=lder1; der1(mj+dvar)=0.0d0
          !        der2(mj)=lder2;der2(mj+2*dvar)=0.0d0;der2(mj+dvar)=0.0d0;

          !     else if(family(mj)==7) then
          !        call cbb1derivs(uvec(mj),vlat(i,1),th((/mj,mj+dvar/)),ccdf(mj),cder1vec,cder2vec) !C_{ij|V_0}
          !        cder1((/mj,mj+dvar/)) = cder1vec
          !        cder2((/mj,mj+2*dvar,mj+dvar/)) = cder2vec

          !        call lfrk2derivs(ccdf(mj),vlat(i,jg+1),th(2*dvar+mj),&
          !         llpdf(dvar+mj),der1(2*dvar+mj),der2(3*dvar+mj),&
          !                       der1u(mj),der2u(mj),derumix(mj))
    
          !         call lbb1derivs2(uvec(mj),vlat(i,1),th((/mj,mj+dvar/)),llpdf(mj),lder11,lder22,lderu,lderv,lderuv) !c_{ij,V_0}  
          !          der1((/mj,mj+dvar/)) = lder11
          !          der2((/mj,mj+2*dvar,mj+dvar/)) = lder22
                

          !     else if(family(mj)==17) then
          !        call cbb1derivs(1.0d0-uvec(mj),1.0d0-vlat(i,1),th((/mj,mj+dvar/)),1.0d0-ccdf(mj),cder1vec,cder2vec) !C_{ij|V_0}
          !        cder1((/mj,mj+dvar/)) = -cder1vec
          !        cder2((/mj,mj+2*dvar,mj+dvar/)) = -cder2vec

          !        call lfrk2derivs(ccdf(mj),vlat(i,jg+1),th(2*dvar+mj),&
          !         llpdf(dvar+mj),der1(2*dvar+mj),der2(3*dvar+mj),&
          !                       der1u(mj),der2u(mj),derumix(mj))
                
          !        call lbb1derivs2(1.0d0-uvec(mj),1.0d0-vlat(i,1),th((/mj,mj+dvar/)),llpdf(mj),lder11,lder22,lderu,lderv,lderuv) !c_{ij,V_0}
          !          der1((/mj,mj+dvar/)) = lder11
          !          der2((/mj,mj+2*dvar,mj+dvar/)) = lder22
                 
             
          !     else if(family(mj)==10) then
          !        call cbb8derivs(uvec(mj),vlat(i,1),th((/mj,mj+dvar/)),ccdf(mj),cder1vec,cder2vec)!C_{ij|V_0}
          !        cder1((/mj,mj+dvar/)) = cder1vec
          !        cder2((/mj,mj+2*dvar,mj+dvar/)) = cder2vec

          
          !        call lfrk2derivs(ccdf(mj),vlat(i,jg+1),th(2*dvar+mj),&
          !         llpdf(dvar+mj),der1(2*dvar+mj),der2(3*dvar+mj),&
          !                       der1u(mj),der2u(mj),derumix(mj))
               
          !        call lbb8derivs(uvec(mj),vlat(i,1),th((/mj,mj+dvar/)),llpdf(mj),lder11,lder22,lderu,lderv,lderuv) !c_{ij,V_0} 
          !          der1((/mj,mj+dvar/)) = lder11
          !          der2((/mj,mj+2*dvar,mj+dvar/)) = lder22

          !     else if (family(mj)==14) then

          !      call cgumderivs(1.0d0-uvec(mj),1.0d0-vlat(i,1),th(mj),1.0d0-ccdf(mj),ccder1,ccder2)
          !      cder1vec(1)=-ccder1;cder1vec(2)=0.0d0;
          !      cder2vec(1)=-ccder2;cder2vec(2)=0.0d0;cder2vec(3)=0.0d0
          !      cder1((/mj,mj+dvar/)) = cder1vec
          !      cder2((/mj,mj+2*dvar,mj+dvar/)) = cder2vec
                
          !     call lfrk2derivs(ccdf(mj),vlat(i,jg+1),th(2*dvar+mj),&
          !       llpdf(dvar+mj),der1(2*dvar+mj),der2(3*dvar+mj),&
          !                       der1u(mj),der2u(mj),derumix(mj))   

               
          !      call lgum2derivs(1.0d0-uvec(mj),1.0d0-vlat(i,1),th(mj),llpdf(mj),lder1,lder2,lder1u,lder2u,ldermix)
          !      der1(mj)=lder1; der1(mj+dvar)=0.0d0
          !      der2(mj)=lder2;der2(mj+2*dvar)=0.0d0;der2(mj+dvar)=0.0d0;
                
          !     else if(family(mj)==20) then
          !        call cbb8derivs(1.0d0-uvec(mj),1.0d0-vlat(i,1),th((/mj,mj+dvar/)),1.0d0-ccdf(mj),cder1vec,cder2vec)!C_{ij|V_0}!C_{ij|V_0}
          !        cder1((/mj,mj+dvar/)) = -cder1vec
          !        cder2((/mj,mj+2*dvar,mj+dvar/)) = -cder2vec
                  
    
          !        call lfrk2derivs(ccdf(mj),vlat(i,jg+1),th(2*dvar+mj),&
          !         llpdf(dvar+mj),der1(2*dvar+mj),der2(3*dvar+mj),&
          !                       der1u(mj),der2u(mj),derumix(mj))
               
                  
          !        call lbb8derivs(1.0d0-uvec(mj),1.0d0-vlat(i,1),th((/mj,mj+dvar/)),llpdf(mj),lder11,lder22,lderu,lderv,lderuv) !c_{ij,V_0}
          !          der1((/mj,mj+dvar/)) = lder11
          !          der2((/mj,mj+2*dvar,mj+dvar/)) = lder22  
              
             
          !     end if 
                       

                       
! subroutine isinf(A)
!     double precision, intent(in) :: A
!     !print*, "Test (A-1 == A)"
!     if (A-1 .eq. A) then
!         !print*, "    INFINITY!!!"
!         return TRUE
!     else
!      !   print*, "    NOT infinite"
!      return FALSE
!     endif
! end subroutine
! ! ! ! outer product of two vectors
! ! ! input
! ! !   na = length of avec
! ! !   nb = length of bvec
! ! !   avec = vector 1
! ! !   bvec = vector 2
! ! ! output
! ! !   abmat = na x nb matrix with outer product of avec and bvec 
! ! subroutine outer(na,nb,avec,bvec,abmat)
! !   implicit none
! !   integer na,nb,i
! !   double precision avec(na),bvec(nb),abmat(na,nb)
! !   ! na=size(avec); nb=size(bvec)
! !   do i =1,na  
! !     abmat(i,:)=avec(i)*bvec
! !   end do
! !   return
! !   end


! subroutine lfrk2derivs(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermix)
!   implicit none
!   double precision u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermix
!   double precision den,den1,den2,den1u,den2u,denmix,t0,t1,t2

!   t0 = exp(-cpar);
!   t1 = exp(-cpar*u1);
!   t2 = exp(-cpar*u2);
!   den = t1+t2-t0-t1*t2;
!   den1 = -u1*t1-u2*t2+t0+(u1+u2)*t1*t2;
!   den2 = u1*u1*t1+u2*u2*t2-t0-(u1+u2)*(u1+u2)*t1*t2;
!   den1u = -cpar*t1*(1.d0-t2);  
!   den2u = cpar*cpar*t1*(1.d0-t2);
!   denmix = t1*(-1.d0+cpar*u1+t2-(u1+u2)*cpar*t2); 
   
!   !lpdf = log(abs(cpar))+log(abs(1-t0))-cpar*(u1+u2)-2*log(abs(den));
!   ! 121023 maybe later add the limits as cpar->0
!   ! pdf = cpar*(1-t0)/den^2 where
!   !    1-t0 has same sign as cpar,  den has same sign as cpar
!   lpdf = log(cpar*(1.d0-t0))-cpar*(u1+u2)-2.d0*log(abs(den));
!   lder1 = 1.d0/cpar+t0/(1.d0-t0)-(u1+u2)-2.d0*den1/den;
!   lder2 = -1.d0/(cpar*cpar)-t0/((1.d0-t0)*(1.d0-t0))-2.d0*den2/den+2.d0*den1*den1/(den*den); 
!   lder1u = -cpar-2.d0*den1u/den;
!   lder2u = -2.d0*den2u/den+2.d0*den1u*den1u/(den*den);
!   ldermix = -1.d0-2.d0*denmix/den+2.d0*den1u*den1/(den*den);
!   return
!   end


! subroutine lgum2derivs(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermix)
!   implicit none
!   double precision u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermix
!   double precision x,y,tx,ty,xd,yd,s,m,logs,logm,msq,dlsq,dlcu
!   double precision sder1,sder2,mder1,mder2,den,den2
!   double precision mu,m2u,u1sq,muder1
!   x = -log(u1); y = -log(u2);
!   tx = log(x); ty = log(y); 
!   xd = x**cpar; yd = y**cpar; s = xd+yd; m = s**(1.d0/cpar);
!   logs=log(s);
!   dlsq=cpar*cpar; dlcu=dlsq*cpar
!   ! for 1-factor and 2-factor models
!   sder1 = xd*tx+yd*ty;
!   sder2 = xd*tx*tx+yd*ty*ty;
!   mder1 = m*sder1/(s*cpar)-m*logs/dlsq;
!   mder2 = -mder1*logs/dlsq-2.d0*m*sder1/(s*dlsq)+2.d0*m*logs/dlcu;
!   mder2 = mder2+sder2*m/(s*cpar)+(mder1/s-m*sder1/(s*s))*sder1/cpar;
!   den = m+cpar-1.d0; den2=den*den
!   logm=log(m); msq=m*m
!   lpdf = -m+log(den)+(1.d0-2.d0*cpar)*logm+(cpar-1.d0)*(tx+ty)+x+y;
!   lder1 = -mder1+(mder1+1.d0)/den-2.d0*logm+(1.d0-2.d0*cpar)*mder1/m+tx+ty;
!   lder2 = -mder2+mder2/den-(mder1+1.d0)**2/den2-4.d0*mder1/m+(1.d0-2.d0*cpar)*(mder2/m-mder1*mder1/msq);
!   ! for 2-factor model
!   u1sq=u1*u1;
!   mu = -m*xd/(u1*s*x); 
!   m2u = (1.d0-cpar)*m*xd*xd/(u1*s*x)**2+(cpar-1.d0)*m*xd/(s*x*x*u1sq)+m*xd/(s*x*u1sq);
!   muder1 = -(mder1/s-m*sder1/(s*s))*xd/(x*u1)-m*xd*tx/(s*u1*x);
!   lder1u = -mu+mu/den+(1.d0-2.d0*cpar)*mu/m-(cpar-1.d0)/(u1*x)-1.d0/u1;
!   lder2u = -m2u+m2u/den-mu*mu/den2+(1.d0-2.d0*cpar)*(m2u/m-mu*mu/msq)+(cpar-1.d0)/(u1sq*x); 
!   lder2u = lder2u-(cpar-1.d0)/(x*x*u1sq)+1.d0/u1sq;
!   ldermix = -muder1+muder1/den-mu*(mder1+1.d0)/den2-2.d0*mu/m+(1.d0-2.d0*cpar)*(muder1/m-mu*mder1/msq)-1.d0/(x*u1);
!   return
!   end


! subroutine cgumderivs(u1,u2,cpar,ccdf,cder1,cder2)
!   implicit none
!   double precision u1,u2,cpar,ccdf,cder1,cder2
!   double precision x,y,tx,ty,xd,yd,s,m,logs,logm,msq,dlsq,dlcu
!   double precision sder1,sder2,mder1,mder2
!   double precision lccdf,lcder1,lcder2
!   x = -log(u1); y = -log(u2);
!   tx = log(x); ty = log(y); 
!   xd = x**cpar; yd = y**cpar; s = xd+yd; m = s**(1.d0/cpar);
!   logs=log(s);
!   dlsq=cpar*cpar; dlcu=dlsq*cpar
!   sder1 = xd*tx+yd*ty;
!   sder2 = xd*tx*tx+yd*ty*ty;
!   mder1 = m*sder1/(s*cpar)-m*logs/dlsq;
!   mder2 = -mder1*logs/dlsq-2.d0*m*sder1/(s*dlsq)+2.d0*m*logs/dlcu;
!   mder2 = mder2+sder2*m/(s*cpar)+(mder1/s-m*sder1/(s*s))*sder1/cpar;
!   logm=log(m); msq=m*m
!   lccdf = y-m+(1.d0-cpar)*(logm-ty)
!   lcder1 = -mder1+(1.d0-cpar)*mder1/m-logm+tx;
!   lcder1 = -mder1+(1.d0-cpar)*mder1/m-logm+ty;  ! given y
!   lcder2 = -mder2-2.d0*mder1/m+(1.d0-cpar)*(mder2/m-mder1*mder1/msq);
!   ccdf = exp(lccdf)
!   cder1 = ccdf*lcder1
!   cder2 = ccdf*(lcder1*lcder1+lcder2)
!   return
!   end



! subroutine clfrkderivs(u1,u2,cpar,ccdf,cder1,cder2)
!   implicit none
!   double precision u1,u2,cpar,ccdf,cder1,cder2
!   double precision t0,t1,t2,den,den1,den2,lccdf,lcder1,lcder2
  
!   t0 = exp(-cpar);
!   t1 = exp(-cpar*u1);
!   t2 = exp(-cpar*u2);
!   den = t1+t2-t0-t1*t2;
!   den1 = -u1*t1-u2*t2+t0+(u1+u2)*t1*t2;
!   den2 = u1*u1*t1+u2*u2*t2-t0-(u1+u2)*(u1+u2)*t1*t2;
!   lccdf = -cpar*u2+log(1-t1)-log(den);
!   lccdf = -cpar*u2+log((1.d0-t1)/den);
!   lcder1 = -u2+u1*t1/(1.d0-t1)-den1/den;
!   lcder2 = -u1*u1*t1/((1.d0-t1)*(1.d0-t1))-den2/den+den1*den1/(den*den);
!   ccdf = exp(lccdf)
!   cder1 = ccdf*lcder1
!   cder2 = ccdf*(lcder1*lcder1+lcder2)
!   return
!   end


! subroutine cbb1derivs(u1,u2,cparv,ccdf,cder11,cder22)
!   implicit none
!   double precision u1,u2,cparv(2),ccdf,cder11(2),cder22(3)
!   double precision theta,delta
!   double precision m,mder1th,mder1dl,mder2th,mder2dl,mderthd
!   double precision t2,tu2,cf1,logm
!   double precision mp1,msq,thsq,dl1n,t2a,lt2a,lmp1
!   double precision lcdf,lder1th,lder1dl,lder2th,lder2dl,lderthd
!   double precision cder1th,cder1dl,cder2th,cder2dl,cderthd

!   theta=cparv(1); delta=cparv(2)

!   call mderivs(u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd)
!   mp1=1.d0+m; msq=m*m; thsq=theta*theta
!   lmp1=log(mp1); logm=log(m);

!   t2 = u2**(-theta); tu2 = -log(u2);
!   t2a=t2-1.d0; lt2a=log(t2a)
!   cf1 = 1.d0+1.d0/theta;
!   dl1n=1.d0-delta

!   lcdf = -cf1*lmp1+(dl1n)*(logm-lt2a)+(theta+1.d0)*tu2;
!   lder1th = lmp1/thsq-cf1*mder1th/mp1;
!   lder1th = lder1th+(dl1n)*(mder1th/m-t2*tu2/t2a)+tu2;
!   lder1dl = -cf1*mder1dl/mp1+(dl1n)*mder1dl/m-logm+lt2a;
!   lder2th = -2.d0*lmp1/(thsq*theta)+(2.d0/thsq)*mder1th/mp1-cf1*(mder2th/mp1&
!   -mder1th**2/(mp1*mp1));
!   lder2th = lder2th+(dl1n)*(mder2th/m-mder1th**2/msq+tu2*tu2*t2/(t2a*t2a));
!   lder2dl = -cf1*(mder2dl/mp1-mder1dl**2/(mp1*mp1))-2.d0*mder1dl/m;
!   lder2dl = lder2dl+(dl1n)*(mder2dl/m-mder1dl**2/msq);
!   lderthd = (1.d0/thsq)*mder1dl/mp1-cf1*(mderthd/mp1-mder1th*mder1dl/(mp1*mp1))-mder1th/m;
!   lderthd = lderthd+(dl1n)*(mderthd/m-mder1th*mder1dl/msq)+t2*tu2/t2a;

!   ccdf = exp(lcdf);
!   cder1th=ccdf*lder1th;
!   cder1dl=ccdf*lder1dl;
!   cder2dl=ccdf*(lder2dl+lder1dl**2);
!   cder2th=ccdf*(lder2th+lder1th**2);
!   cderthd=ccdf*(lderthd+lder1dl*lder1th);
!   cder11(1)=cder1th; cder11(2)=cder1dl  ! gradient of log pdf
!   cder22(1)=cder2th; cder22(2)=cderthd; cder22(3)=cder2dl; ! Hessian terms
!   return
! end




! subroutine cbb8derivs(u1,u2,cparv,ccdf,cder11,cder22)
!  implicit none
!   double precision u1,u2,cparv(2),ccdf,cder11(2),cder22(3)
!   double precision nu,delta,x,y,tem,eta,eta1nu,eta1del
!   double precision x1nu,x1del,y1nu,y1del,tem1nu,tem1del,mm
!   double precision tem1,tem2,tem3,tem4,eta2nu,eta2del,eta2nudel
!   double precision x2nu,y2nu,x2del,y2del,x2nudel,y2nudel
!   double precision tem2nu,tem2del, tem2nudel
!   double precision lpdf,lpdf1nu,lpdf1del,lpdf2del,lpdf2nu,lpdf2nudel
!   double precision cder1nu,cder1del,cder2nu,cder2del,xy1mixdel,cder2nudel
  
!   nu=cparv(1)
!   delta=cparv(2)
!   eta=1.0d0-(1.0d0-delta)**nu
!   x=1.0d0-(1.0d0-delta*u1)**nu
!   y=1.0d0-(1.0d0-delta*u2)**nu
!   tem=x*y/eta
!   mm=1.0d0/nu-1.0d0

!   eta1nu=-log(1.0d0-delta)*((1.0d0-delta)**nu) !\p eta/\p nu
!   eta1del=nu*((1.0d0-delta)**(nu-1.0d0))!\p eta/\p delta
!   x1nu=-log(1.0d0-delta*u1)*((1-delta*u1)**nu) !\p x/\p nu
!   x1del=(nu*(1.0d0-delta*u1)**(nu-1))*u1!\p x/\p delta
!   y1nu=-log(1.0d0-delta*u2)*((1-delta*u2)**nu)!\p y/\p nu
!   y1del=nu*((1.0d0-delta*u2)**(nu-1))*u2!\p y/\p delta
!   tem1nu=-(x*y/eta**2)*eta1nu+(x1nu*y+y1nu*x)/eta  !\p tem/\p nu
!   tem1del=(x1del*y+y1del*x)/eta-eta1del*x*y/eta**2!\p tem/\p delta
 
!   eta2nu=-((log(1.0d0-delta))**2)*((1.0d0-delta)**nu) !\p^2 eta/\p nu^2
!   eta2del=-nu*(nu-1.0d0)*(1.0d0-delta)**(nu-2.0d0) !\p^2 eta/\p delta^2
!  !\p^2 eta/\p delta*nu
!   eta2nudel=nu*log(1-delta)*((1.0d0-delta)**(nu-1.0d0))+(1.0d0-delta)**(nu-1.0d0)
!   x2nu=-(log(1.0d0-delta*u1)**2)*((1.0d0-delta*u1)**nu)!\p^2 x/\p nu^2
!   y2nu=-(log(1.0d0-delta*u2)**2)*((1.0d0-delta*u2)**nu)!\p^2 y/\p nu^2
!   x2del=-u1*u1*nu*(nu-1.0d0)*((1-delta*u1)**(nu-2.0d0))!\p^2 x/\p delta^2
!   y2del=-u2*u2*nu*(nu-1.0d0)*((1-delta*u2)**(nu-2.0d0))!\p^2 y/\p delta^2
!  !\p^2 x/\p delta*nu
!   x2nudel=u1*nu*(1-delta*u1)**(nu-1.0d0)*log(1-delta*u1)+u1*((1-delta*u1)**nu)/(1-delta*u1)
! !\p^2 y/\p delta*nu
!   y2nudel=u2*nu*(1-delta*u2)**(nu-1.0d0)*log(1-delta*u2)+u2*((1-delta*u2)**nu)/(1-delta*u2)
 
!   tem1=(2.0d0/eta**3)*x*y*((eta1nu)**2)-(eta2nu*x*y+(x1nu*y+y1nu*x)*eta1nu)/eta**2
!   tem2=-eta1nu*x1nu*y/eta**2+(x2nu*y+x1nu*y1nu)/eta
!   tem3=-eta1nu*y1nu*x/eta**2+(y2nu*x+x1nu*y1nu)/eta
!   tem2nu=tem1+tem2+tem3 !\p^2 tem/\p nu^2
 
!   xy1mixdel=x1del*y1del
!   tem1=(2.0d0*x*y*((eta1del)**2))/eta**3-(eta2del*x*y+(x1del*y+y1del*x)*eta1del)/eta**2
!   tem2=-eta1del*y*x1del/eta**2+(xy1mixdel+x2del*y)/eta
!   tem3=-eta1del*x*y1del/eta**2+(xy1mixdel+y2del*x)/eta
!   tem2del=tem1+tem2+tem3!\p^2 tem/\p delta^2
!   tem1=2.0d0*eta1nu*eta1del*x*y/eta**3.0d0
!   tem2=-(eta2nudel*x*y+(x1del*y+y1del*x)*eta1nu)/eta**2.0d0
!   tem3=-eta1del*x1nu*y/eta**2+(x2nudel*y+y1del*x1nu)/eta
!   tem4=-eta1del*y1nu*x/eta**2+(y2nudel*x+x1del*y1nu)/eta
!   tem2nudel=tem1+tem2+tem3+tem4 !\p^2 tem/\p nu*delta
 
!   lpdf=-log(eta)+log(x)+(1.0d0/nu-1.0d0)*log(1.0d0-tem)+(1.0d0-1.0d0/nu)*log(1.0d0-y)
!   ccdf=exp(lpdf)
!  ccdf=x*(1-tem)**(1.0d0/nu-1.0d0)*(1.0d0-y)**(1.0d0-1.0d0/nu)/eta
 
!   tem1=-eta1nu/eta+x1nu/x-log(1.0d0-tem)/nu**2-tem1nu*mm/(1.0d0-tem)
!   tem2=log(1.0d0-y)/nu**2+y1nu*mm/(1.0d0-y)
!   lpdf1nu=tem1+tem2
!   cder1nu=ccdf*lpdf1nu !\p ccdf/\p nu
 
!   lpdf1del=-eta1del/eta+x1del/x+(1.0d0/nu-1.0d0)*(-tem1del)/(1.0d0-tem)&
!   +(1.0d0-1.0d0/nu)*(-y1del)/(1.0d0-y)
!   cder1del=ccdf*lpdf1del!\p ccdf/\p delta
 
!   tem1=(eta1nu)**2/eta**2-eta2nu/eta-(x1nu)**2/x**2+x2nu/x
!   tem2=2*log(1.0d0-tem)/nu**3+tem1nu/(1.0d0-tem)/nu**2
!   tem3=-mm*tem1nu**2/(1.0d0-tem)**2-(tem2nu*mm-tem1nu/nu**2)/(1.0d0-tem)
!   tem4=-2.0d0*log(1.0d0-y)/nu**3-y1nu/(1.0d0-y)/nu**2-(-y2nu*mm/(1.0d0-y)&
!   +(1.0d0/(1.0d0-y)/nu**2-y1nu*mm/(1.0d0-y)**2)*y1nu)
!   lpdf2nu=tem1+tem2+tem3+tem4 !\p^2 lpdf/\p nu^2
!   cder2nu=ccdf*(lpdf1nu)**2+lpdf2nu*ccdf!\p^2 ccdf/\p nu^2
 
!   tem1=eta1del**2/eta**2-eta2del/eta-x1del**2/x**2+x2del/x
!   tem2=(-tem1del**2/(1.0d0-tem)**2-tem2del/(1.0d0-tem))*(1.0d0/nu-1.0d0)
!   tem3=-mm*(-y2del/(1.0d0-y)-y1del**2/(1.0d0-y)**2)
!   lpdf2del=tem1+tem2+tem3!\p^2 lpdf/\p delta^2
!   cder2del=ccdf*(lpdf1del)**2+lpdf2del*ccdf!\p^2 ccdf/\p delta^2
 
!   tem1=eta1nu*eta1del/eta**2-eta2nudel/eta
!   tem2=-x1del*x1nu/x**2+x2nudel/x+tem1del/(1.0d0-tem)/nu**2
!   tem3=-(mm*tem1nu*tem1del/(1.0d0-tem)**2+mm*tem2nudel/(1.0d0-tem))
!   tem4=-y1del/(1.0d0-y)/nu**2+mm*(y2nudel/(1.0d0-y)+y1del*y1nu/(1.0d0-y)**2)
!   lpdf2nudel=tem1+tem2+tem3+tem4!\p^2 lpdf/\p delta*nu
!   cder2nudel=ccdf*lpdf1nu*lpdf1del+lpdf2nudel*ccdf!\p^2 ccdf/\p nu*delta
!   cder11(1)=cder1nu;cder11(2)=cder1del
!   cder22(1)=cder2nu;cder22(2)=cder2nudel;cder22(3)=cder2del
! return
! end



! subroutine lbb1derivs2(u1,u2,cparv,lpdf,lder11,lder22,lderu,lderv,lderuv)
!   implicit none
!   double precision u1,u2,cparv(2),lpdf,lder11(2),lder22(3)
!   double precision theta,delta,der1th,der1dl,der2th,der2dl,derthd
!   double precision t10,t1der1th,t1der1dl,t1der2th,t1der2dl,t1derthd
!   double precision t20,t2der1th,t2der1dl,t2der2th,t2der2dl,t2derthd
!   double precision t30,t3der1th,t3der1dl,t3der2th,t3der2dl,t3derthd
!   double precision t40,t4der1th,t4der1dl,t4der2th,t4der2dl,t4derthd
!   double precision m,mder1th,mder1dl,mder2th,mder2dl,mderthd,mderu(4),mderv(4)
!   double precision den,coef,t1,t2,tu1,tu2
!   double precision mp1,msq,thsq,den2,thtem,dltem,dl1,t1a,t2a,smlog
!   double precision mder1u,mder2u,mder2uth,mder2udl,mder1v,mder2v,mder2vth,mder2vdl,mder2uv
!   double precision lderu(4),lderv(4),lderuv,tem,tem1u,tem1v,lpdf1u,lpdf2u,lpdf2uth,lpdf2udl
!   double precision lpdf2uv,lpdf1v,lpdf2v,lpdf2vth,lpdf2vdl,dl2,th1
!   double precision tem1,tem2,tem3,th2,tem1th,tem1dl,tem2uth,tem2vth
!   theta=cparv(1); delta=cparv(2);

  
!   call  mderivs2(u1,u2,cparv,m,mder1th,mder1dl,mder2th,mder2dl,mderthd,mderu,mderv,mder2uv)
  
!   mder1u=mderu(1);mder2u=mderu(2);mder2uth=mderu(3);mder2udl=mderu(4)
!   mder1v=mderv(1);mder2v=mderv(2);mder2vth=mderv(3);mder2vdl=mderv(4);
!   mp1=1.d0+m; msq=m*m; thsq=theta*theta
!   thtem=2.d0+1.d0/theta
!   t10 = -(thtem)*log(mp1)
!   t1der1th = log(mp1)/thsq - (thtem)*mder1th/mp1
!   t1der1dl = -(thtem)*mder1dl/mp1
!   t1der2th = -2.d0*log(mp1)/theta**3+(2.d0/thsq)*mder1th/mp1
!   t1der2th = t1der2th-(thtem)*(mder2th/mp1-mder1th*mder1th/(mp1*mp1))
!   t1der2dl = -(thtem)*(mder2dl/mp1-mder1dl*mder1dl/(mp1*mp1))
!   t1derthd = mder1dl/(thsq*(mp1))-(thtem)*(mderthd/(mp1)-mder1th*mder1dl/(mp1*mp1))
  
!   tem=theta*(delta-1.0d0)+(theta*delta+1.d0)*m
!   tem1u=(theta*delta+1.0d0)*mder1u
!   tem1v=(theta*delta+1.0d0)*mder1v
!   tem1th=delta-1.0d0+mder1th+delta*(m+mder1th*theta)
!   tem1dl=theta+mder1dl+theta*(m+mder1dl*delta)
!   tem2uth=delta*mder1u+mder2uth*(theta*delta+1.0d0)
!   tem2vth=delta*mder1v+mder2vth*(theta*delta+1.0d0)
  
!   dltem=1.d0-2.d0*delta; dl1=delta-1.d0
!   t20 = (dltem)*log(m)
!   t2der1th = (dltem)*mder1th/m
!   t2der1dl = -2.d0*log(m)+(dltem)*mder1dl/m
!   t2der2th = (dltem)*(mder2th/m-mder1th*mder1th/msq)
!   t2der2dl = -4.d0*mder1dl/m+(dltem)*(mder2dl/m-mder1dl*mder1dl/msq)
!   t2derthd = -2.d0*mder1th/m+(dltem)*(mderthd/m-mder1th*mder1dl/msq)

!   coef = theta*delta+1.d0
!   den = theta*(dl1)+(coef)*m
!   den2=den*den
!   t30 = log(theta*(dl1)+coef*m)
!   t3der1th = (dl1+delta*m+coef*mder1th)/den
!   t3der1dl = (theta+theta*m+coef*mder1dl)/den
!   t3der2th = (2.d0*delta*mder1th+coef*mder2th)/den
!   t3der2th = t3der2th-(dl1+delta*m+coef*mder1th)**2/den2
!   t3der2dl = (2.d0*theta*mder1dl+coef*mder2dl)/den
!   t3der2dl = t3der2dl-(theta+theta*m+coef*mder1dl)**2/den2;
!   t3derthd = (1.d0+m+theta*mder1th+delta*mder1dl+coef*mderthd)/den;
!   t3derthd = t3derthd - (theta+theta*m+coef*mder1dl)*(dl1+delta*m+coef*mder1th)/den2

!   t1 = u1**(-theta); t2 = u2**(-theta)
!   t1a=t1-1.d0; t2a=t2-1.d0; smlog=log(t1a)+log(t2a)
!   tu1 = -log(u1); tu2 = -log(u2)
!   t40 = (dl1)*(smlog)+(theta+1.d0)*(tu1+tu2)
!   t4der1th = (dl1)*(t1*tu1/t1a+t2*tu2/t2a)+tu1+tu2  
!   t4der1dl = smlog
!   t4der2th = -(dl1)*(t1*tu1*tu1/t1a**2+t2*tu2*tu2/t2a**2)
!   t4der2dl = 0.d0
!   t4derthd = t1*tu1/t1a+t2*tu2/t2a

!   lpdf = t10+t20+t30+t40
!   der1th = t1der1th+t2der1th+t3der1th+t4der1th
!   der1dl = t1der1dl+t2der1dl+t3der1dl+t4der1dl
!   der2th = t1der2th+t2der2th+t3der2th+t4der2th
!   der2dl = t1der2dl+t2der2dl+t3der2dl+t4der2dl
!   derthd = t1derthd+t2derthd+t3derthd+t4derthd
!   lder11(1)=der1th; lder11(2)=der1dl  ! gradient of log pdf
!   lder22(1)=der2th; lder22(2)=derthd; lder22(3)=der2dl; ! Hessian terms
  
!   th2=(-1.0d0/theta)-2.0d0
!   dl2=1.0d0-2.0d0*delta
!   th1=-theta-1.0d0
  
!   tem1=th2*mder1u/(m+1.0d0)+(dl2)*mder1u/m
!   tem2=coef*mder1u/tem+dl1*(-theta*u1**(-theta-1.0d0)/(u1**(-theta)-1.0d0))&
!   -(theta+1.0d0)/u1
!   lpdf1u=tem1+tem2 !correct
  
!   tem1=th2*mder1v/(m+1.0d0)+(dl2)*mder1v/m
!   tem2=coef*mder1v/tem+dl1*(-theta*u2**(-theta-1.0d0)/(u2**(-theta)-1.0d0))&
!   -(theta+1.0d0)/u2
!   lpdf1v=tem1+tem2 !correct
  
!   tem1=th2*(-mder1u**2.0d0/(m+1.0d0)**2.0d0+mder2u/(m+1.0d0))&
!   +dl2*(-mder1u**2.0d0/m**2+mder2u/m)
!   tem2=coef*(-tem1u*mder1u/tem**2.0d0+mder2u/tem)+(theta+1.0d0)/u1**2.0d0
!   tem3=(dl1/(u1**(-theta)-1)**2.d0)*(-theta*th1*u1**(th1-1.0d0)*(u1**(-theta)-1.0d0)&
!   -(theta*u1**th1)**2.0d0)
!   lpdf2u=tem1+tem2+tem3 !correct
  
  
!   tem1=th2*(-mder1v**2.0d0/mp1**2.0d0+mder2v/mp1)+dl2*(-mder1v**2.0d0/m**2+mder2v/m)
!   tem2=coef*(-tem1v*mder1v/tem**2.0d0+mder2v/tem)+(theta+1.0d0)/u2**2.0d0
!   tem3=(dl1/(u2**(-theta)-1)**2.d0)*(-theta*th1*u2**(th1-1.0d0)*(u2**(-theta)-1.0d0)&
!   -(theta*u2**th1)**2.0d0)
!   lpdf2v=tem1+tem2+tem3 !correct
  
!   tem1=th2*(-mder1u*mder1v/mp1**2+mder2uv/mp1)+dl2*(-mder1u*mder1v/m**2+mder2uv/m)
!   tem2=coef*(-tem1v*mder1u/tem**2+mder2uv/tem)
!   lpdf2uv=tem1+tem2!correct
  
!   tem1=mder1u/(mp1*thsq)+(-1.0d0/theta-2.0d0)*(-mder1th*mder1u/mp1**2+mder2uth/mp1)&
!   +dltem*(-mder1th*mder1u/m**2+mder2uth/m)
!   tem2=(-tem1th*tem1u/tem/tem)+tem2uth/tem
!   tem3=(-t1/u1+t1*log(u1)*theta/u1)*t1a-t1*log(u1)*theta*t1/u1
!   lpdf2uth=tem1+tem2+dl1*tem3/t1a**2-1.0d0/u1
  
!   tem1=mder1v/(mp1*thsq)+(-1.0d0/theta-2.0d0)*(-mder1th*mder1v/mp1**2+mder2vth/mp1)&
!   +dltem*(-mder1th*mder1v/m**2+mder2vth/m)
!   tem2=(-tem1th*tem1v/tem/tem)+tem2vth/tem
!   tem3=(-t2/u2+t2*log(u2)*theta/u2)*t2a-t2*log(u2)*theta*t2/u2
!   lpdf2vth=tem1+tem2+dl1*tem3/t2a**2-1.0d0/u2
  
!   tem1=th2*(-mder1dl*mder1u/mp1**2+mder2udl/mp1)-2.0d0*mder1u/m
!   tem2=dl2*(-mder1dl*mder1u/m**2+mder2udl/m)-tem1dl*coef*mder1u/tem**2.0d0
!   tem3=(theta*mder1u+coef*mder2udl)/tem
!   lpdf2udl=tem1+tem2+tem3-theta*(t1/u1)/(t1-1.0d0)!correct
  
!   tem1=th2*(-mder1dl*mder1v/mp1**2+mder2vdl/mp1)-2.0d0*mder1v/m
!   tem2=dl2*(-mder1dl*mder1v/m**2+mder2vdl/m)-tem1dl*coef*mder1v/tem**2.0d0
!   tem3=(theta*mder1v+coef*mder2vdl)/tem
!   lpdf2vdl=tem1+tem2+tem3-theta*(t2/u2)/(t2-1.0d0)!correct
  
!   lderu(1)=lpdf1u;lderu(2)=lpdf2u;lderu(3)=lpdf2uth;lderu(4)=lpdf2udl;
!   lderv(1)=lpdf1v;lderv(2)=lpdf2v;lderv(3)=lpdf2vth;lderv(4)=lpdf2vdl;
  
!   lderuv=lpdf2uv
!   return
! end





! subroutine  lgum2derivt(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,&
! ldermixu,lder1v,lder2v,lder2uv,ldermixv)
!   implicit none
!   double precision u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,lder2uv
!   double precision x,y,tx,ty,xd,yd,s,m,logs,logm,msq,dlsq,dlcu
!   double precision sder1,sder2,mder1,mder2,den,den2
!   double precision mu,m2u,u1sq,muder1,u2sq,lder1v,lder2v,ldermixv,ldermixu
!   double precision mv,m2v,mvder1,m2uv
!   x = -log(u1); y = -log(u2);
!   tx = log(x); ty = log(y); 
!   xd = x**cpar; yd = y**cpar; s = xd+yd; m = s**(1.d0/cpar);
!   logs=log(s);
!   dlsq=cpar*cpar; dlcu=dlsq*cpar
!   !for 1-factor and 2-factor models
!   sder1 = xd*tx+yd*ty;
!   sder2 = xd*tx*tx+yd*ty*ty;
!   mder1 = m*sder1/(s*cpar)-m*logs/dlsq;
!   mder2 = -mder1*logs/dlsq-2.d0*m*sder1/(s*dlsq)+2.d0*m*logs/dlcu;
!   mder2 = mder2+sder2*m/(s*cpar)+(mder1/s-m*sder1/(s*s))*sder1/cpar;
!   den = m+cpar-1.d0; den2=den*den
!   logm=log(m); msq=m*m
!   lpdf = -m+log(den)+(1.d0-2.d0*cpar)*logm+(cpar-1.d0)*(tx+ty)+x+y;
!   lder1 = -mder1+(mder1+1.d0)/den-2.d0*logm+(1.d0-2.d0*cpar)*mder1/m+tx+ty;
!   lder2 = -mder2+mder2/den-(mder1+1.d0)**2/den2-4.d0*mder1/m&
!   +(1.d0-2.d0*cpar)*(mder2/m-mder1*mder1/msq);
!   !for 2-factor model
!   u1sq=u1*u1;
!   mu = -m*xd/(u1*s*x); 
!   m2u = (1.d0-cpar)*m*xd*xd/(u1*s*x)**2+(cpar-1.d0)*m*xd/(s*x*x*u1sq)+m*xd/(s*x*u1sq);
!   muder1 = -(mder1/s-m*sder1/(s*s))*xd/(x*u1)-m*xd*tx/(s*u1*x);
!   lder1u = -mu+mu/den+(1.d0-2.d0*cpar)*mu/m-(cpar-1.d0)/(u1*x)-1.d0/u1;
!   lder2u = -m2u+m2u/den-mu*mu/den2+(1.d0-2.d0*cpar)*(m2u/m-mu*mu/msq)+(cpar-1.d0)/(u1sq*x); 
!   lder2u = lder2u-(cpar-1.d0)/(x*x*u1sq)+1.d0/u1sq;
!   ldermixu = -muder1+muder1/den-mu*(mder1+1.d0)/den2-2.d0*mu/m&
!   +(1.d0-2.d0*cpar)*(muder1/m-mu*mder1/msq)-1.d0/(x*u1);
  
!   u2sq=u2*u2;
!   mv = -m*yd/(u2*s*y); 
!   m2v = (1.d0-cpar)*m*yd*yd/(u2*s*y)**2+(cpar-1.d0)*m*yd/(s*y*y*u2sq)+m*yd/(s*y*u2sq);
!   mvder1 = -(mder1/s-m*sder1/(s*s))*yd/(y*u2)-m*yd*ty/(s*u2*y);
!   lder1v = -mu+mu/den+(1.d0-2.d0*cpar)*mu/m-(cpar-1.d0)/(u2*y)-1.d0/u2;
!   lder2v = -m2u+m2u/den-mu*mu/den2+(1.d0-2.d0*cpar)*(m2u/m-mu*mu/msq)+(cpar-1.d0)/(u2sq*y); 
!   lder2v = lder2v-(cpar-1.d0)/(y*y*u2sq)+1.d0/u2sq;
!   ldermixv = -muder1+muder1/den-mu*(mder1+1.d0)/den2-2.d0*mu/m&
!   +(1.d0-2.d0*cpar)*(muder1/m-mu*mder1/msq)-1.d0/(y*u2);
  
!   m2uv=(x**(cpar-1.0d0)/u1/u2)*(cpar*(1.0d0/cpar-1.0d0)*m*y**(cpar-1.0d0)/s/s)
!   lder2uv=-m2uv+(-mu*mv/den**2+m2uv/den)+(1.0d0-2*cpar)*(m2uv/m-mu*mv/m**2)
!   return
!   end
  


! subroutine lfrk2derivt(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu,&
!  lder1v,lder2v,lder2uv,ldermixv)
!   implicit none
!   double precision u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu
!   double precision lder1v,lder2v,lder2uv,ldermixv
!   double precision den,den1,den2,den1u,den2u,denmixu,t0,t1,t2
!   double precision den1v,den2v,denmixv

!   t0 = exp(-cpar);
!   t1 = exp(-cpar*u1);
!   t2 = exp(-cpar*u2);
!   den = t1+t2-t0-t1*t2;
!   den1 = -u1*t1-u2*t2+t0+(u1+u2)*t1*t2;
!   den2 = u1*u1*t1+u2*u2*t2-t0-(u1+u2)*(u1+u2)*t1*t2;
!   den1u = -cpar*t1*(1.d0-t2);  
!   den1v = -cpar*t2*(1.d0-t1);   ! added
!   den2u = cpar*cpar*t1*(1.d0-t2);
!   den2v = cpar*cpar*t2*(1.d0-t1); ! added
!   denmixu = t1*(-1.d0+cpar*u1+t2-(u1+u2)*cpar*t2); 
!   denmixv = t2*(-1.d0+cpar*u2+t1-(u1+u2)*cpar*t1); ! added 
   
!   lpdf = log(abs(cpar))+log(abs(1-t0))-cpar*(u1+u2)-2*log(abs(den));
!   !121023 maybe later add the limits as cpar->0
!   !pdf = cpar*(1-t0)/den^2 where
!   !   1-t0 has same sign as cpar,  den has same sign as cpar
!   lpdf = log(cpar*(1.d0-t0))-cpar*(u1+u2)-2.d0*log(abs(den));
!   lder1 = 1.d0/cpar+t0/(1.d0-t0)-(u1+u2)-2.d0*den1/den;
!   lder2 = -1.d0/(cpar*cpar)-t0/((1.d0-t0)*(1.d0-t0))-2.d0*den2/den&
!   +2.d0*den1*den1/(den*den); 
!   lder1u = -cpar-2.d0*den1u/den;
!   lder1v = -cpar-2.d0*den1v/den;  ! added
!   lder2u = -2.d0*den2u/den+2.d0*den1u*den1u/(den*den);
!   lder2v = -2.d0*den2v/den+2.d0*den1v*den1v/(den*den); ! added
!   lder2uv = 2.d0*cpar*cpar*t1*t2/den+2.d0*den1u*den1v/(den*den); ! added
!   ldermixu = -1.d0-2.d0*denmixu/den+2.d0*den1u*den1/(den*den);
!   ldermixv = -1.d0-2.d0*denmixv/den+2.d0*den1v*den1/(den*den); ! added
!   return
!   end



! subroutine lbb8derivs(u1,u2,cparv,lpdf,lder11,lder22,lderu,lderv,lderuv)
!  implicit none
!  double precision u1,u2,cparv(2),lpdf,lder11(2),lder22(3)
!  double precision eta,x,y,delta,dl1,dlu1,dlu2,nu,nu1,nu2
!  double precision eta1nu,eta1del,x1nu,x1del,y1nu,y1del,eta2nu,eta2del,eta2nudel
!  double precision tem,temm,tem1,tem2,tem3,tem4,lpdf1nu,lpdf1del,tem1nu,tem1dl,eta1dl
!  double precision x2nu,y2nu,y2nudel,y2del,x2nudel,x2del,lpdf2nu,lpdf2del,lpdf2nudel
!  double precision xy1mixdel,tem2nu,tem2del,tem2nudel
!  double precision lderu(4),lderv(4),lderuv,lpdf1u,lpdf2u,lpdf1v,lpdf2v
!  double precision lpdf2udl,lpdf2unu,lpdf2vdl,lpdf2vnu
!  double precision tem1u,tem1v,x1u,y1v,tem2u,tem2v,x2u,y2v,x2unu,y2vnu
!  double precision x2udl,y2vdl,tem2unu,tem2udl,tem2vnu,tem2vdl
!  double precision tem2uv
!  nu=cparv(1)
!  delta=cparv(2)
!  dl1=1.0d0-delta
!  nu1=nu-1.0d0
!  nu2=1.0d0/nu-2.0d0
!  dlu1=delta*u1
!  dlu2=delta*u2
!  eta=1.0d0-dl1**nu
!  x=1.0d0-(1.0d0-delta*u1)**nu
!  y=1.0d0-(1.0d0-delta*u2)**nu
!  tem=x*y/eta !x*y/eta
!  temm=1.0d0-tem
!  lpdf=-log(eta)+log(delta)+(nu2)*log(temm)+log(nu-tem)+(nu1)*log(1.0d0-dlu1)+&
!  (nu1)*log(1-dlu2)
 
!  eta1nu=-log(dl1)*((dl1)**nu)!\p eta/\p nu
!  eta1del=nu*((dl1)**(nu1)) !\p eta/\p delta
!  x1nu=-log(1.0d0-dlu1)*((1-dlu1)**nu) !\p x/\p nu
!  x1del=(nu*(1.0d0-dlu1)**(nu1))*u1 !\p x/\p delta
!  y1nu=-log(1.0d0-dlu2)*((1-delta*u2)**nu)  !\p y/\p nu
!  y1del=nu*((1.0d0-dlu2)**(nu1))*u2 !\p y/\p delta
!  tem1nu=-(x*y/eta**2)*eta1nu+(x1nu*y+y1nu*x)/eta !\p tem/\p nu
!  tem1dl=(x1del*y+y1del*x)/eta-eta1del*x*y/eta**2!\p tem/\p delta
 
!  tem1=-eta1nu/eta-log(temm)/nu**2+nu2*(1.0d0/(temm))*(-tem1nu)
!  tem2=(1.0d0/(nu-tem))*(1.0d0-tem1nu)+log(1.0d0-dlu1)+log(1.0d0-dlu2)
!  lpdf1nu=tem1+tem2!\p lpdf/\p nu
 
!  tem1=-eta1del/eta+1.0d0/delta+nu2*(1.0d0/temm)*(-tem1dl)
!  tem2=-(1.0d0/(nu-tem))*(tem1dl)-(nu1)*u1/(1.0d0-delta*u1)-nu1*u2/(1.0d0-delta*u2)
!  lpdf1del=tem1+tem2 !\p lpdf/\p delta
 
!  lder11(1)=lpdf1nu;lder11(2)=lpdf1del;
 
!  eta2nu=-((log(dl1))**2)*(dl1**nu) !\p^2 eta/\p nu^2
!  eta2del=-nu*(nu1)*(dl1)**(nu-2.0d0) !\p^2 eta/\p delta^2
!  eta2nudel=nu*log(dl1)*(dl1**nu1)+dl1**nu1!\p^2 eta/\p nu delta
!  x2nu=-(log(1.0d0-dlu1)**2)*((1.0d0-dlu1)**nu) !\x^2 eta/\p nu^2
!  y2nu=-(log(1.0d0-delta*u2)**2)*((1.0d0-dlu2)**nu)!\y^2 eta/\p nu^2
!  x2del=-u1*u1*nu*nu1*((1-dlu1)**(nu-2.0d0))!\x^2 eta/\p delta^2
!  y2del=-u2*u2*nu*nu1*((1-dlu2)**(nu-2.0d0))!\y^2 eta/\p delta^2
!  x2nudel=u1*nu*(1-dlu1)**(nu1)*log(1-dlu1)+u1*((1-dlu1)**nu)/(1-dlu1) !\x^2 eta/\p nu delta
!  y2nudel=u2*nu*(1-dlu2)**(nu1)*log(1-dlu2)+u2*((1-dlu2)**nu)/(1-dlu2)!\y^2 eta/\p nu delta
 
!  tem1=(2.0d0/eta**3)*x*y*((eta1nu)**2)-(eta2nu*x*y+(x1nu*y+y1nu*x)*eta1nu)/eta**2
!  tem2=-eta1nu*x1nu*y/eta**2+(x2nu*y+x1nu*y1nu)/eta
!  tem3=-eta1nu*y1nu*x/eta**2+(y2nu*x+x1nu*y1nu)/eta
!  tem2nu=tem1+tem2+tem3 !\p^2 tem/\p nu^2
 
!  xy1mixdel=x1del*y1del
!  tem1=(2.0d0*x*y*((eta1del)**2))/eta**3-(eta2del*x*y+(x1del*y+y1del*x)*eta1del)/eta**2
!  tem2=-eta1del*y*x1del/eta**2+(xy1mixdel+x2del*y)/eta
!  tem3=-eta1del*x*y1del/eta**2+(xy1mixdel+y2del*x)/eta
!  tem2del=tem1+tem2+tem3!\p^2 tem/\p delta^2
 
 
!  tem1=2.0d0*eta1nu*eta1del*x*y/eta**3.0d0
!  tem2=-(eta2nudel*x*y+(x1del*y+y1del*x)*eta1nu)/eta**2.0d0
!  tem3=-eta1del*x1nu*y/eta**2+(x2nudel*y+y1del*x1nu)/eta
!  tem4=-eta1del*y1nu*x/eta**2+(y2nudel*x+x1del*y1nu)/eta
!  tem2nudel=tem1+tem2+tem3+tem4 !\p^2 tem/\p nu delta
 
!  tem1=eta1nu**2/eta**2-eta2nu/eta+2.0d0*log(temm)/nu**3.0d0
!  tem2=tem1nu/(temm)/nu**2.0d0+tem1nu/(temm)/nu**2.0d0
!  tem3=(1.0d0/nu-2.0d0)*(-(tem1nu**2.0d0/(temm)**2.0d0)-tem2nu/(temm))
!  tem4=-(1.0d0-tem1nu)**2.0d0/(nu-tem)**2.0d0-tem2nu/(nu-tem)
!  lpdf2nu=tem1+tem2+tem3+tem4!\p^2 lpdf/\p nu^2 
 
!  tem1=eta1del**2/eta**2-eta2del/eta-1.0d0/delta**2
!  !tem2=(-(1.0d0/(1.0d0-tem)**2)*(tem1del)**2+(-tem2del)/(1.0d0-tem)))*(1.0d0/nu-2.0d0)
!  tem2=nu2*(-tem1dl**2/(temm)**2-tem2del/temm)
!  tem3=-(tem1dl**2/(nu-tem)**2+tem2del/(nu-tem))
!  tem4=-(nu1*u1**2.0d0/(1-dlu1)**2.0d0+nu1*u2**2.0d0/(1-dlu2)**2.0d0)
!  lpdf2del=tem1+tem2+tem3+tem4!\p^2 lpdf/\p delta^2 
 
!  tem1=eta1nu*eta1del/eta**2-eta2nudel/eta+tem1dl/nu**2/(1.0d0-tem)
!  tem2=nu2*(-(tem1nu*tem1dl/(1-tem)**2)-tem2nudel/(1.0d0-tem))
!  tem3=(1.0d0-tem1nu)*tem1dl/(nu-tem)**2-tem2nudel/(nu-tem)
!  tem4=-u1/(1.0d0-dlu1)-u2/(1.0d0-dlu2)
!  lpdf2nudel=tem1+tem2+tem3+tem4
 
!  lder22(1)=lpdf2nu
!  lder22(2)=lpdf2nudel
!  lder22(3)=lpdf2del
 
!  x1u=delta*nu*(1-dlu1)**nu1
!  y1v=delta*nu*(1-dlu2)**nu1
!  tem1u=y*x1u/eta
!  tem1v=x*y1v/eta
 
!  lpdf1u=-nu2*tem1u/temm-tem1u/(nu-tem)-nu1*delta/(1-dlu1)
!  lpdf1v=-nu2*tem1v/temm-tem1v/(nu-tem)-nu1*delta/(1-dlu2)
 
!  x2u=-(delta**2)*nu*nu1*(1-dlu1)**(nu-2.0d0)
!  y2v=-(delta**2)*nu*nu1*(1-dlu2)**(nu-2.0d0)
!  tem2u=y*x2u/eta
!  tem2v=x*y2v/eta
 
!  lpdf2u=-nu2*tem1u**2/temm**2-nu2*(tem2u/temm)-tem1u**2/(nu-tem)**2&
!  -tem2u/(nu-tem)-nu1*delta**2/(1-dlu1)**2
!  lpdf2v=-nu2*tem1v**2/temm**2-nu2*(tem2v/temm)-tem1v**2/(nu-tem)**2&
!  -tem2v/(nu-tem)-nu1*delta**2/(1-dlu2)**2
 
!  x2unu=delta*(1-dlu1)**nu1+delta*nu*log(1-dlu1)*(1-dlu1)**nu1
!  tem2unu=-eta1nu*y*x1u/eta**2+(y1nu*x1u+x2unu*y)/eta 
!  tem1=tem1u/nu**2/temm+(-tem1nu*tem1u/temm**2-tem2unu/temm)*nu2
!  tem2=tem1u*(1-tem1nu)/(nu-tem)**2-tem2unu/(nu-tem)-delta/(1-dlu1)
!  lpdf2unu=tem1+tem2
 
!  eta1dl=nu*dl1**nu1
!  x2udl=nu*(1.0d0-dlu1)**nu1-delta*nu*nu1*u1*(1.0d0-dlu1)**(nu1-1.0d0)
!  tem2udl=-eta1dl*y*x1u/eta**2+(y1del*x1u+x2udl*y)/eta
!  tem1=-nu2*tem1dl*tem1u/temm**2-nu2*tem2udl/temm
!  tem2=-tem1dl*tem1u/(nu-tem)**2-tem2udl/(nu-tem)-nu1/(1-dlu1)**2
!  lpdf2udl=tem1+tem2
 
!  y2vnu=delta*(1-dlu2)**nu1+delta*nu*log(1-dlu2)*(1-dlu2)**nu1
!  tem2vnu=-eta1nu*x*y1v/eta**2+(x1nu*y1v+y2vnu*x)/eta 
!  tem1=tem1v/nu**2/temm+(-tem1nu*tem1v/temm**2-tem2vnu/temm)*nu2
!  tem2=tem1v*(1-tem1nu)/(nu-tem)**2-tem2vnu/(nu-tem)-delta/(1-dlu2)
!  lpdf2vnu=tem1+tem2
 
!  y2vdl=nu*(1.0d0-dlu2)**nu1-delta*nu*nu1*u2*(1.0d0-dlu2)**(nu1-1.0d0)
!  tem2vdl=-eta1dl*x*y1v/eta**2+(x1del*y1v+y2vdl*x)/eta
!  tem1=-nu2*tem1dl*tem1v/temm**2-nu2*tem2vdl/temm
!  tem2=-tem1dl*tem1v/(nu-tem)**2-tem2vdl/(nu-tem)-nu1/(1-dlu2)**2
!  lpdf2vdl=tem1+tem2
 
!  tem2uv=y1v*x1u/eta
!  tem1=-nu2*tem1u*tem1v/temm**2-tem2uv*nu2/temm
!  tem2=-tem1u*tem1v/(nu-tem)**2-tem2uv/(nu-tem)
!  lderuv=tem1+tem2
 
 
!  lderu(1)=lpdf1u;lderu(2)=lpdf2u;lderu(3)=lpdf2unu;lderu(4)=lpdf2udl;
!  lderv(1)=lpdf1v;lderv(2)=lpdf2v;lderv(3)=lpdf2vnu;lderv(4)=lpdf2vdl;
! return
! end
 





! subroutine mderivs2(u1,u2,cparv,m,mder1th,mder1dl,mder2th,mder2dl,mderthd,&
! mderu,mderv,mder2uv)
!   implicit none
!   double precision u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd
!   double precision t1,t2,tu1,tu2,ttu1,ttu2,td01,td02,td11,td12,td21,td22
!   double precision s,sder1th,sder1dl,sder2th,sder2dl,sderthd,m1,ts,m1der1dl 
!   double precision t1a,t2a,dlsq,dlcu
!   double precision xder1u,yder1v,xder2u,mder1u,mder1v,mder2u,mder2uv,mder2v,del1,th1
!   double precision yder2v,xder2uth,yder2vth,xder1th,xder2udl,yder2vdl,yder1th
!   double precision mder2vdl,mder2vth,mder2udl,mder2uth
!   double precision mderu(4),mderv(4),cparv(2)
!   theta=cparv(1);delta=cparv(2)
!   del1=1.0d0/delta-1.0d0
!   th1=-theta-1.0d0
  
!   t1 = u1**(-theta); t2 = u2**(-theta)
!   t1a=t1-1.d0; t2a=t2-1.d0
!   tu1 = -log(u1); tu2 = -log(u2)
!   ttu1 = log(t1a); ttu2 = log(t2a)
!   td01 = (t1a)**delta; td02 = (t2a)**delta
!   td11=td01/t1a; td12=td02/t2a;
!   td21=td11/t1a; td22=td12/t2a;
 
!   s = td01+td02
!   sder1th = delta*(td11*t1*tu1+td12*t2*tu2)
!   sder1dl = td01*ttu1+td02*ttu2
!   sder2th = delta*(delta-1.d0)*(td21*t1*t1*tu1*tu1+td22*t2*t2*tu2*tu2)
!   sder2th = sder2th+delta*(td11*t1*tu1*tu1+td12*t2*tu2*tu2)
!   sder2dl = td01*ttu1*ttu1+td02*ttu2*ttu2
!   sderthd = sder1th/delta+delta*(td11*ttu1*tu1*t1+td12*ttu2*tu2*t2)

!   m = s**(1.d0/delta); m1 = m/s
!   ts = log(s)
!   dlsq=delta*delta; dlcu=delta*dlsq
!   mder1th = m1*sder1th/delta
!   mder1dl = m1*sder1dl/delta - m*ts/dlsq
!   m1der1dl = mder1dl/s - m*sder1dl/s**2
!   mder2th = (1.d0-delta)*m1*sder1th**2/(dlsq*s)+m1*sder2th/delta
!   mder2dl = 2.d0*m*ts/dlcu-mder1dl*ts/dlsq-2.d0*m1*sder1dl/dlsq
!   mder2dl = mder2dl+sder2dl*m1/delta+sder1dl*m1der1dl/delta
!   mderthd = -m1*sder1th/dlsq+sder1th*m1der1dl/delta+m1*sderthd/delta
  
!   xder1u=-delta*theta*td11*u1**(th1) !correct
!   yder1v=-delta*theta*td12*u2**(th1) !correct
!   xder2u=-delta*theta*(-theta*(delta-1.0d0)*td21*(u1**(2*th1))+&
!   th1*td11*u1**(th1-1.0d0))!correct
!   yder2v=-delta*theta*(-theta*(delta-1.0d0)*td22*(u2**(2*th1))+&
!   th1*td12*u2**(th1-1.0d0))!correct
!   xder2uth=-delta*td11*u1**th1&
!   -delta*theta*(-(delta-1.0d0)*td21*u1**(-theta)*log(u1)*u1**th1&
!   -u1**th1*log(u1)*td11)!correct
!   yder2vth=-delta*td12*u2**th1-&
!   delta*theta*(-(delta-1.0d0)*td22*u2**(-theta)*log(u2)*u2**th1&
!   -u2**th1*log(u2)*td12)!correct
!   xder1th=delta*(td11*t1*tu1)
!   yder1th=delta*(td12*t2*tu2)
!   xder2udl=-theta*u1**th1*(td11+delta*td11*ttu1)
!   yder2vdl=-theta*u2**th1*(td12+delta*td12*ttu2)
  
!   mder1u=xder1u*s**del1/delta !correct
!   mder1v=yder1v*s**del1/delta !correct
!   mder2uv=del1*(s**(del1-1))*yder1v*xder1u/delta!correct
!   mder2u=del1*s**(del1-1.0d0)*xder1u**2/delta+xder2u*s**del1/delta!correct
!   mder2v=del1*(s**(del1-1.0d0))*(yder1v)**2/delta+yder2v*s**(del1)/delta!correct
!   mder2uth=(xder2uth*s**del1+del1*s**(del1-1.0d0)*sder1th*xder1u)/delta !correct
!   mder2udl=xder2udl*s**del1/delta+(-s**del1/delta**2+&
!   (1.0d0/delta)*s**(del1)*(-log(s)/delta**2+sder1dl*del1/s))*xder1u
!   mder2vth=(yder2vth*s**del1+del1*s**(del1-1.0d0)*sder1th*yder1v)/delta!correct
!   mder2vdl=yder2vdl*s**del1/delta+(-s**del1/delta**2+&
!   (1.0d0/delta)*s**(del1)*(-log(s)/delta**2+sder1dl*del1/s))*yder1v
  
!   mderu(1)=mder1u;mderu(2)=mder2u;mderu(3)=mder2uth;mderu(4)=mder2udl;
!   mderv(1)=mder1v;mderv(2)=mder2v;mderv(3)=mder2vth;mderv(4)=mder2vdl;

!   return
  
!   end



! subroutine mderivs(u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd)
!   implicit none
!   double precision u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd
!   double precision t1,t2,tu1,tu2,ttu1,ttu2,td01,td02,td11,td12,td21,td22
!   double precision s,sder1th,sder1dl,sder2th,sder2dl,sderthd,m1,ts,m1der1dl 
!   double precision t1a,t2a,dlsq,dlcu

!   t1 = u1**(-theta); t2 = u2**(-theta)
!   t1a=t1-1.d0; t2a=t2-1.d0
!   tu1 = -log(u1); tu2 = -log(u2)
!   ttu1 = log(t1a); ttu2 = log(t2a)
!   td01 = (t1a)**delta; td02 = (t2a)**delta
!   td11=td01/t1a; td12=td02/t2a;
!   td21=td11/t1a; td22=td12/t2a;
 
!   s = td01+td02
!   sder1th = delta*(td11*t1*tu1+td12*t2*tu2)
!   sder1dl = td01*ttu1+td02*ttu2
!   sder2th = delta*(delta-1.d0)*(td21*t1*t1*tu1*tu1+td22*t2*t2*tu2*tu2)
!   sder2th = sder2th+delta*(td11*t1*tu1*tu1+td12*t2*tu2*tu2)
!   sder2dl = td01*ttu1*ttu1+td02*ttu2*ttu2
!   sderthd = sder1th/delta+delta*(td11*ttu1*tu1*t1+td12*ttu2*tu2*t2)

!   m = s**(1.d0/delta); m1 = m/s
!   ts = log(s)
!   dlsq=delta*delta; dlcu=delta*dlsq
!   mder1th = m1*sder1th/delta
!   mder1dl = m1*sder1dl/delta - m*ts/dlsq
!   m1der1dl = mder1dl/s - m*sder1dl/s**2
!   mder2th = (1.d0-delta)*m1*sder1th**2/(dlsq*s)+m1*sder2th/delta
!   mder2dl = 2.d0*m*ts/dlcu-mder1dl*ts/dlsq-2.d0*m1*sder1dl/dlsq
!   mder2dl = mder2dl+sder2dl*m1/delta+sder1dl*m1der1dl/delta
!   mderthd = -m1*sder1th/dlsq+sder1th*m1der1dl/delta+m1*sderthd/delta
!   return
! end
