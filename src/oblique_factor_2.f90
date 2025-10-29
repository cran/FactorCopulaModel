! 
! 
!  !main program
! PROGRAM test3   ! main function
!   implicit none
!   integer npar,mgrp,dvar,n,i,nq
!   integer,dimension(:),allocatable::fam
!   double precision,dimension(:,:),allocatable::udata,lat,hess
!   double precision,dimension(:),allocatable::th,grad
!   double precision,dimension(:),allocatable::xl,wl
!   integer,dimension(:),allocatable::grsize
!   double precision nllk
!   
!      read *,npar
! 
!      print*,npar
! 
!      read *,dvar
! 
!      print*,dvar
!      read *,mgrp
!      print*,mgrp
! 
!     read*,nq
!     print*,nq
! 
!     allocate(grsize(mgrp))
!     allocate(fam(dvar))
!     allocate(xl(nq),wl(nq))
! 
!     read*,xl
! 
!     print*,xl
! 
! 
!     read*,wl
!     print*,wl
! 
!     read *,grsize
!     print*,grsize
! 
!     read*,fam
!     print*,fam
! 
!     read*,n
!     print*,n
! 
!     allocate(th(npar),udata(n,dvar),lat(n,mgrp))
!     allocate(grad(npar),hess(npar,npar))
!     read *,th
!     print*,th
! 
!     read *, (udata(i,:), i=1,n) 
! 
!     do i=1,n
!       print *, udata(i,:)
!     end do
! 
!     read *, (lat(i,:), i=1,n) 
! 
!     do i=1,n
!        print *, lat(i,:)
!     end do
!     
!   !   !call bifactproxy(npar,th,mgrp,n,dvar,family,grsize,udata,vlat,nllk,grad,hess)
!    !call strbb1frk2proxy(npar,th,mgrp,n,dvar,grsize,udata,vlat,nllk,grad,hess)
!    call strfrkdif(npar,th,mgrp,n,dvar,fam,grsize,udata,lat,nq,wl,xl,nllk,grad,hess)
!   
!     print *, "the negative log-likelihood is"
!     print *,nllk
! 
!     print*,"the gradient is"
!     print*, grad
! 
!     print*,"the hessian matrix is"
!     print*,hess
!   deallocate(grsize,udata,lat,th,grad,hess)
! END PROGRAM


subroutine strfrkdif(npar,th,mgrp,n,dvar,fam,grsize,udata,lat,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,mgrp,dvar,n,nq,ip,jp, ind
  integer fam(dvar)
  double precision th(npar),udata(n,dvar),wl(nq),xl(nq),uvec(dvar),lat(n,mgrp)
  !double precision lderu(4), lderv(4), lderuv
  double precision lpdf(dvar+mgrp),der1(npar),der2(npar+dvar),llk, lder11(2),lder22(3) !der1bb1(2), der2bb1(3)
  double precision nllk,liki,lk,grad(npar),hess(npar,npar)
  double precision lder1u,lder2u,ldermixuu(2),lder1v,lder2v,lder2uv,ldermixvv(2)
  integer i,iq, jg,jg2,mj,mj2,nj,nj2,ind1,ind2,grsize(mgrp),ibb1(2)
  double precision intj(mgrp), intjj(npar), grdj(npar), der2j(npar), int0, grd0(npar), grdi(npar)
  double precision dermxj(dvar), hss0(npar,npar), hssi(npar,npar), hessi(npar,npar), hssaux(npar,npar)

  ! npar = 2*dvar+mgrp; dvar = sum(grsize)
  nllk=0.d0; grad = 0.d0; hess=0.d0 
  do i =1,n 
    uvec = udata(i,:)    
    liki = 0.d0; grdi = 0.d0; hssi = 0.d0; 
    do iq =1,nq
      int0 = 1.d0; grd0 = 1.d0; dermxj = 0.d0; hssaux = 0.d0 
      ind = 0; intj = 0.d0;  grdj = 0.d0; der2j = 0.d0  
      do jg =1,mgrp   !jth group
        !jg2=jg+dvar                         
        jg2=jg
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
      
          lk = 1.d0 
          do mj = ind1,ind2
            nj=mj+mgrp
            call lcop2derivt(uvec(mj),lat(i,jg),th((/nj,nj+dvar/)),fam(mj),lpdf(nj),lder11,lder22,lder1u,lder2u,&
              ldermixuu,lder1v,lder2v,lder2uv,ldermixvv)
            !call lbb1derivs(uvec(mj),lat(i,jg),th((/nj,nj+dvar/)),lpdf(nj),der1bb1,der2bb1,lderu,lderv,lderuv)   !c_{ij,V_j}   
            !der1((/nj,nj+dvar/)) = der1bb1
            !der2((/nj,nj+2*dvar,nj+dvar/)) = der2bb1
            der1((/nj,nj+dvar/))=lder11
            der2((/nj,nj+2*dvar,nj+dvar/))=lder22
            lk = lk*exp(lpdf(nj)) 
            end do
          call lfrk1derivs(lat(i,jg),xl(iq),th(jg2),lpdf(jg2),der1(jg2),der2(jg2))    !c_{V_j,V_0}  
          llk = exp(lpdf(jg2))*lk !llk: values of the jth integrand
          !intj(jg) = intj(jg) + wl(iq2)*llk !intj value of the jth integral
          intj(jg)=llk
          intjj(ind1+mgrp:ind2+mgrp) = intj(jg) !intjj: repeat intj(j) grsize[j] times
          do mj = ind1,ind2
            nj=mj+mgrp
            ibb1 = (/nj,nj+dvar/)
            grdj(ibb1) = llk*der1(ibb1)  !grdj: der. of the jth int. wrt th(mj)
            der2j(ibb1)= llk*(der2(ibb1)+der1(ibb1)*der1(ibb1)) !der2j: second order der. wrt th(mj)^2 (th(m+mj)^2)
            dermxj(mj) = llk*(der2(nj+2*dvar)+der1(nj)*der1(nj+dvar)) !dermxj: mixed der. wrt th(mj)=th and th(m+mj)=dl
            hssaux(jg2,ibb1) = hssaux(jg2,ibb1) + llk*der1(ibb1)*der1(jg2)
            hssaux(ibb1,jg2) = hssaux(jg2,ibb1)
            if (mj < ind2) then
              do mj2 = (mj+1),ind2  !2nd order der. of the jth int. wrt th(mj), th(mj2)
                nj2=mj2+mgrp                                  
                hssaux(nj,nj2) = hssaux(nj,nj2) + llk*der1(nj)*der1(nj2) !th,th
                hssaux(nj+dvar,nj2+dvar) = hssaux(nj+dvar,nj2+dvar) + llk*der1(nj+dvar)*der1(nj2+dvar)  !dl,dl
                hssaux(nj+dvar,nj2) = hssaux(nj+dvar,nj2) + llk*der1(nj+dvar)*der1(nj2) !dl,th
                hssaux(nj,nj2+dvar) = hssaux(nj,nj2+dvar) + llk*der1(nj)*der1(nj2+dvar) !th,dl
                hssaux(nj2,nj) = hssaux(nj,nj2)
                hssaux(nj2+dvar,nj+dvar) = hssaux(nj+dvar,nj2+dvar)
                hssaux(nj2,nj+dvar) = hssaux(nj+dvar,nj2)
                hssaux(nj2+dvar,nj) = hssaux(nj,nj2+dvar)
              end do 
            endif         
          end do
          grdj(jg2) = llk*der1(jg2) 
          der2j(jg2)= llk*(der2(jg2)+der1(jg2)*der1(jg2))
      
      end do
      intjj(1+dvar+mgrp:dvar+dvar+mgrp)=intjj(1+mgrp:mgrp+dvar)
      intjj(1:mgrp)=intj
      int0 = product(intj) !product of j inner integrals  
      grd0 = int0*(grdj/intjj) !gradient of the product
      do ip=2,npar !hessian of the product
        do jp=1,ip-1
          hss0(ip,jp) = int0*grdj(ip)*grdj(jp)/intjj(ip)/intjj(jp)
          hss0(jp,ip) = hss0(ip,jp)
        end do
      end do
      !ind = 0
      ind = mgrp
      do jg=1,mgrp  !2nd order derivatives within groups  
        !jg2=jg+dvar
        jg2=jg                           
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg)
        hss0(ind1:ind2,ind1:ind2)=int0*hssaux(ind1:ind2,ind1:ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,dvar+ind1:dvar+ind2)=int0*hssaux(dvar+ind1:dvar+ind2,dvar+ind1:dvar+ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,ind1:ind2)=int0*hssaux(dvar+ind1:dvar+ind2,ind1:ind2)/intj(jg)
        hss0(ind1:ind2,dvar+ind1:dvar+ind2)=int0*hssaux(ind1:ind2,dvar+ind1:dvar+ind2)/intj(jg)
        hss0(ind1:ind2,jg2)=int0*hssaux(ind1:ind2,jg2)/intj(jg)
        hss0(jg2,ind1:ind2)=int0*hssaux(jg2,ind1:ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,jg2)=int0*hssaux(dvar+ind1:dvar+ind2,jg2)/intj(jg)
        hss0(jg2,dvar+ind1:dvar+ind2)=int0*hssaux(jg2,dvar+ind1:dvar+ind2)/intj(jg)
      end do

      do ip=1,npar
        hss0(ip,ip) = int0*der2j(ip)/intjj(ip)
      end do
      do ip=1+mgrp,dvar+mgrp
        hss0(ip,ip+dvar) = int0*dermxj(ip-mgrp)/intjj(ip)
        hss0(ip+dvar,ip) = hss0(ip,ip+dvar)
      end do
      liki = liki + wl(iq)*int0   !integrating over V0
      grdi = grdi + wl(iq)*grd0
      hssi = hssi + wl(iq)*hss0
    end do
    nllk = nllk - log(liki)       !updating loglikelihood
    grad = grad - grdi/liki       !updating gradient
    do ip=1,npar                  !updating hessian
      do jp=1,ip
        hessi(ip,jp)=hssi(ip,jp)/liki-grdi(ip)*grdi(jp)/liki/liki 
        hessi(jp,ip)=hessi(ip,jp)
      end do
    end do 
    hess = hess - hessi 
  end do
  return
  end

! subroutine strfrkdif(npar,th,mgrp,n,dvar,fam,grsize,udata,lat,nq,wl,xl,nllk,grad,hess)
!   implicit none
!   integer npar,mgrp,dvar,n,nq,ip,jp,ind,fam(dvar),jj0(2),dd,jj2(2),jj3(2)
!   double precision lderu(4),lderv(4),lderuv
!   double precision lder11(2),lder22(3),tem(2),th_tem(2,1),temmat1(2,2)
!   double precision th(npar),udata(n,dvar),wl(nq),xl(nq),uvec(dvar),lat(n,mgrp)
!   double precision lpdf(npar),der1(npar),der2(dvar+npar+mgrp),llk
!   double precision nllk,liki,lk,grad(npar),hess(npar,npar)
!   !double precision lder1u,lder2u,ldermixuu,lder1v,lder2v,lder2uv,ldermixvv
!   integer i,iq,jg,jg2,mj,mj2,nj,nj2,ind1,ind2,grsize(mgrp)
!   double precision intj(mgrp), intjj(npar), grdj(npar), der2j(npar), int0, grd0(npar), grdi(npar)
!   double precision hss0(npar,npar), hssi(npar,npar), hessi(npar,npar), hssaux(npar,npar)

!   ! npar = dvar+mgrp; dvar = sum(grsize)
!   dd=dvar+mgrp;
!   nllk=0.d0; grad = 0.d0; hess=0.d0 
!   do i =1,n 
!     uvec = udata(i,:)  
!     liki = 0.d0; grdi = 0.d0; hssi = 0.d0; 
!     do iq =1,nq
!       int0 = 1.d0; grd0 = 1.d0; hssaux = 0.d0 
!       ind = 0; intj = 0.d0;  grdj = 0.d0; der2j = 0.d0 
!       do jg =1,mgrp   !jth group       
!       print*,jg                  
!         !jg2=jg+dvar
!         jg2=jg
!         jj2=(/2*jg2-1,2*jg2/)
!         ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
!         !do iq2 =1,nq   
!           lk = 1.d0      
!           do mj = ind1,ind2 ! index for within subgroup
!             nj=mj+mgrp
!             tem=th((/nj,nj+dvar+mgrp/))
!             th_tem=RESHAPE(tem,(/2,1/))
!             call lbb1derivs(uvec(mj),lat(i,jg),th_tem,lpdf(nj),lder11,lder22,lderu,lderv,lderuv)
!             ! call lcop2derivt(uvec(mj),lat(i,jg),th_tem,fam(mj),lpdf(nj),lder11,lder22,&
!             ! lder1u,lder2u,ldermixuu,lder1v,lder2v,lder2uv,ldermixvv) 
!              der1((/nj,nj+dd/))=lder11
!              der2((/nj,nj+2*dd,nj+dd/))=lder22
!             !call lgum1derivs(uvec(mj),lat(i,jg),th(nj),lpdf(nj),der1(nj),der2(nj))   !c_{ij,V_j}   
!             lk = lk*exp(lpdf(nj))
!           end do
!           call lfrk1derivs(lat(i,jg),xl(iq),th(jg2),lpdf(jg2),der1(jg2),der2(jg2))   !c_{V_j,V_0}  
!           der1(jg2+dd)=0.0d0; der2(jg2+2*dd)=0.0d0; der2(jg2+dd)=0.0d0;
!           llk = exp(lpdf(jg2))*lk  !llk: values of the jth integrand
!           intj(jg) = llk !intj value of the jth integral
!           !intjj(ind1:ind2) = intj(jg) !intjj: repeat intj(j) grsize[j] times
!           intjj(ind1+mgrp:ind2+mgrp) = intj(jg) !intjj: repeat intj(j) grsize[j] times
!           intjj(ind1+mgrp+dd:ind2+mgrp+dd) = intj(jg)
!           do mj = ind1,ind2
!             nj=mj+mgrp
!             jj0=(/nj,nj+dd/)
!             grdj(jj0) = llk*der1(jj0) !grdj: der. of the jth int. wrt th(mj)
!             der2j(jj0)= llk*(der2(jj0)+der1(jj0)*der1(jj0)) !der2j: second order der. wrt th(mj)
!             call outer(2,2,der1(jj0),der1(jj2),temmat1)
!             hssaux(jj2,jj0) = llk*temmat1
!             hssaux(jj0,jj2) = hssaux(jj2,jj0)
!             if (mj < ind2) then
!               do mj2 = (mj+1),ind2  !2nd order der. of the jth int. wrt th(mj), th(mj2)
!                 nj2=mj2+mgrp
!                 jj2=(/nj2,nj2+dd/)
!                 call outer(2,2,der1(jj0),der1(jj2),temmat1)
!                 hssaux(jj0,jj2) = llk*temmat1
!                 hssaux(jj2,jj0) = hssaux(jj0,jj2)
!               end do 
!             endif         
!           end do
!           jj3=(/jg2,jg2+dd/)
!           grdj(jj3) = llk*der1(jj3) 
!           der2j(jg2)= llk*(der2(jg2)+der1(jg2)*der1(jg2))
!         end do
!       !end do
!       !intjj(dvar+1:dvar+mgrp)=intj
!       intjj(1:mgrp)=intj
!       int0 = product(intj)   !product of j inner integrals  
!       grd0 = int0*(grdj/intjj) !gradient of the product
!       do ip=2,npar !hessian of the product
!         do jp=1,ip-1
!           hss0(ip,jp) = int0*grdj(ip)*grdj(jp)/intjj(ip)/intjj(jp)
!           hss0(jp,ip) = hss0(ip,jp)
!         end do
!       end do
!       !ind = 0
!       ind = mgrp
!       do jg=1,mgrp   !2nd order derivatives within groups 
!         !jg2=jg+dvar
!         jg2=jg
!         jj2=(/2*jg2-1,2*jg2/)
!         ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg)
!         hss0(ind1:ind2,ind1:ind2)=int0*hssaux(ind1:ind2,ind1:ind2)/intj(jg)
!         hss0(ind1:ind2,jj2)=int0*hssaux(ind1:ind2,jj2)/intj(jg)
!         hss0(jj2,ind1:ind2)=int0*hssaux(jj2,ind1:ind2)/intj(jg)
!       end do
!       do ip=1,npar
!         hss0(ip,ip) = int0*der2j(ip)/intjj(ip)
!       end do
!       liki = liki + wl(iq)*int0   !integrating over V0
!       grdi = grdi + wl(iq)*grd0
!      !  hssi = hssi + wl(iq)*hss0
! !     end do
! !     nllk = nllk - log(liki)       !updating loglikelihood
! !     grad = grad - grdi/liki       !updating gradient
! !     do ip=1,npar                  !updating hessian
! !       do jp=1,ip
! !         hessi(ip,jp)=hssi(ip,jp)/liki-grdi(ip)*grdi(jp)/liki/liki 
! !         hessi(jp,ip)=hessi(ip,jp)
! !       end do
! !     end do 
! !     hess = hess - hessi 
! !   end do
! !   return
! !   end
! 
! 
!   
! ! 1st and 2nd order derivatives of m = s^(1/delta) wrt theta and delta
! ! inputs
! !   u1,u2 = values in (0,1)
! !   theta >0
! !   delta >1
! ! outputs
! !   m = s^(1/delta); s =  [u1^(-theta)-1]^delta + [u2^(-theta)-1]^delta
! !   mder1th = \p m/\p theta
! !   mder1dl = \p m/\p delta
! !   mder2th = \p^2 m/\p theta^2
! !   mder2dl = \p^2 m/\p delta^2
! !   mderthd = \p^2 m/\p theta \p delta
! subroutine mderivs(u1,u2,cparv,m,mder1th,mder1dl,mder2th,mder2dl,mderthd,&
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
!   
!   t1 = u1**(-theta); t2 = u2**(-theta)
!   t1a=t1-1.d0; t2a=t2-1.d0
!   tu1 = -log(u1); tu2 = -log(u2)
!   ttu1 = log(t1a); ttu2 = log(t2a)
!   td01 = (t1a)**delta; td02 = (t2a)**delta
!   td11=td01/t1a; td12=td02/t2a;
!   td21=td11/t1a; td22=td12/t2a;
!  
!   s = td01+td02
!   sder1th = delta*(td11*t1*tu1+td12*t2*tu2)
!   sder1dl = td01*ttu1+td02*ttu2
!   sder2th = delta*(delta-1.d0)*(td21*t1*t1*tu1*tu1+td22*t2*t2*tu2*tu2)
!   sder2th = sder2th+delta*(td11*t1*tu1*tu1+td12*t2*tu2*tu2)
!   sder2dl = td01*ttu1*ttu1+td02*ttu2*ttu2
!   sderthd = sder1th/delta+delta*(td11*ttu1*tu1*t1+td12*ttu2*tu2*t2)
! 
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
!   
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
!   
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
!   
!   mderu(1)=mder1u;mderu(2)=mder2u;mderu(3)=mder2uth;mderu(4)=mder2udl;
!   mderv(1)=mder1v;mderv(2)=mder2v;mderv(3)=mder2vth;mderv(4)=mder2vdl;
! 
!   return
!   
!   end
! 
! ! this routine is in bb1facts.f90
! !Function with bb8 logpdf =log f(u1,u2,th,delta) 
! !lder1, lder2 (partial wrt cparv, 1st and 2nd order) lder1 is 2-dimenisonal vector
! !lder2 is 
!  !also lderu (partial wrt u1, 1st and 2nd order and also wrt parmaters)
! ! ***  also lderv (partial wrt u2, 1st order and 2nd order and also wrt parameters) and others
! ! and lderuv  (partial wrt u and v) 
! ! outputs 
! !lpdf = log pdf, 
! !lder11 = \p lpdf/\p nu  \p lpdf/\p delta !2-dimensional vector
! !lder22 = \p^2 lpdf/\p^2 th^2, \p^2 lpdf/\p^2 th delta, lpdf/\p^2 delta^2!  3d-vector
! !lderu = \p lpdf/\p u, \p^2 lpdf/\p u^2, \p^2 lpdf/\p u th, \p^2 lpdf/\p u delta  4d-vector
! !lderv = \p lpdf/\p v, \p^2 lpdf/\p v^2, \p^2 lpdf/\p v th, \p^2 lpdf/\p v delta  4d-vector
! !lderuv= \p^2 lpdf/\p uv
! subroutine lbb1derivs(u1,u2,cparv,lpdf,lder11,lder22,lderu,lderv,lderuv)
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
! 
!   
!   call  mderivs(u1,u2,cparv,m,mder1th,mder1dl,mder2th,mder2dl,mderthd,mderu,mderv,mder2uv)
!   !print*,mderu
!   !print*,mderv
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
!   
!   tem=theta*(delta-1.0d0)+(theta*delta+1.d0)*m
!   tem1u=(theta*delta+1.0d0)*mder1u
!   tem1v=(theta*delta+1.0d0)*mder1v
!   tem1th=delta-1.0d0+mder1th+delta*(m+mder1th*theta)
!   tem1dl=theta+mder1dl+theta*(m+mder1dl*delta)
!   tem2uth=delta*mder1u+mder2uth*(theta*delta+1.0d0)
!   tem2vth=delta*mder1v+mder2vth*(theta*delta+1.0d0)
!  !  print*,tem
! !   print*,tem1th
! !   print*,tem2uth
!   dltem=1.d0-2.d0*delta; dl1=delta-1.d0
!   t20 = (dltem)*log(m)
!   t2der1th = (dltem)*mder1th/m
!   t2der1dl = -2.d0*log(m)+(dltem)*mder1dl/m
!   t2der2th = (dltem)*(mder2th/m-mder1th*mder1th/msq)
!   t2der2dl = -4.d0*mder1dl/m+(dltem)*(mder2dl/m-mder1dl*mder1dl/msq)
!   t2derthd = -2.d0*mder1th/m+(dltem)*(mderthd/m-mder1th*mder1dl/msq)
! 
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
! 
!   t1 = u1**(-theta); t2 = u2**(-theta)
!   t1a=t1-1.d0; t2a=t2-1.d0; smlog=log(t1a)+log(t2a)
!   tu1 = -log(u1); tu2 = -log(u2)
!   t40 = (dl1)*(smlog)+(theta+1.d0)*(tu1+tu2)
!   t4der1th = (dl1)*(t1*tu1/t1a+t2*tu2/t2a)+tu1+tu2  
!   t4der1dl = smlog
!   t4der2th = -(dl1)*(t1*tu1*tu1/t1a**2+t2*tu2*tu2/t2a**2)
!   t4der2dl = 0.d0
!   t4derthd = t1*tu1/t1a+t2*tu2/t2a
! 
!   lpdf = t10+t20+t30+t40
!   der1th = t1der1th+t2der1th+t3der1th+t4der1th
!   der1dl = t1der1dl+t2der1dl+t3der1dl+t4der1dl
!   der2th = t1der2th+t2der2th+t3der2th+t4der2th
!   der2dl = t1der2dl+t2der2dl+t3der2dl+t4der2dl
!   derthd = t1derthd+t2derthd+t3derthd+t4derthd
!   lder11(1)=der1th; lder11(2)=der1dl  ! gradient of log pdf
!   lder22(1)=der2th; lder22(2)=derthd; lder22(3)=der2dl; ! Hessian terms
!   
!   th2=(-1.0d0/theta)-2.0d0
!   dl2=1.0d0-2.0d0*delta
!   th1=-theta-1.0d0
!   
!   tem1=th2*mder1u/(m+1.0d0)+(dl2)*mder1u/m
!   tem2=coef*mder1u/tem+dl1*(-theta*u1**(-theta-1.0d0)/(u1**(-theta)-1.0d0))&
!   -(theta+1.0d0)/u1
!   lpdf1u=tem1+tem2 !correct
!   
!   tem1=th2*mder1v/(m+1.0d0)+(dl2)*mder1v/m
!   tem2=coef*mder1v/tem+dl1*(-theta*u2**(-theta-1.0d0)/(u2**(-theta)-1.0d0))&
!   -(theta+1.0d0)/u2
!   lpdf1v=tem1+tem2 !correct
!   
!   tem1=th2*(-mder1u**2.0d0/(m+1.0d0)**2.0d0+mder2u/(m+1.0d0))&
!   +dl2*(-mder1u**2.0d0/m**2+mder2u/m)
!   tem2=coef*(-tem1u*mder1u/tem**2.0d0+mder2u/tem)+(theta+1.0d0)/u1**2.0d0
!   tem3=(dl1/(u1**(-theta)-1)**2.d0)*(-theta*th1*u1**(th1-1.0d0)*(u1**(-theta)-1.0d0)&
!   -(theta*u1**th1)**2.0d0)
!   lpdf2u=tem1+tem2+tem3 !correct
!   
!   
!   tem1=th2*(-mder1v**2.0d0/mp1**2.0d0+mder2v/mp1)+dl2*(-mder1v**2.0d0/m**2+mder2v/m)
!   tem2=coef*(-tem1v*mder1v/tem**2.0d0+mder2v/tem)+(theta+1.0d0)/u2**2.0d0
!   tem3=(dl1/(u2**(-theta)-1)**2.d0)*(-theta*th1*u2**(th1-1.0d0)*(u2**(-theta)-1.0d0)&
!   -(theta*u2**th1)**2.0d0)
!   lpdf2v=tem1+tem2+tem3 !correct
!   
!   tem1=th2*(-mder1u*mder1v/mp1**2+mder2uv/mp1)+dl2*(-mder1u*mder1v/m**2+mder2uv/m)
!   tem2=coef*(-tem1v*mder1u/tem**2+mder2uv/tem)
!   lpdf2uv=tem1+tem2!correct
!   
!   tem1=mder1u/(mp1*thsq)+(-1.0d0/theta-2.0d0)*(-mder1th*mder1u/mp1**2+mder2uth/mp1)&
!   +dltem*(-mder1th*mder1u/m**2+mder2uth/m)
!   tem2=(-tem1th*tem1u/tem/tem)+tem2uth/tem
!   tem3=(-t1/u1+t1*log(u1)*theta/u1)*t1a-t1*log(u1)*theta*t1/u1
!   lpdf2uth=tem1+tem2+dl1*tem3/t1a**2-1.0d0/u1
!   
!   tem1=mder1v/(mp1*thsq)+(-1.0d0/theta-2.0d0)*(-mder1th*mder1v/mp1**2+mder2vth/mp1)&
!   +dltem*(-mder1th*mder1v/m**2+mder2vth/m)
!   tem2=(-tem1th*tem1v/tem/tem)+tem2vth/tem
!   tem3=(-t2/u2+t2*log(u2)*theta/u2)*t2a-t2*log(u2)*theta*t2/u2
!   lpdf2vth=tem1+tem2+dl1*tem3/t2a**2-1.0d0/u2
!   
!   tem1=th2*(-mder1dl*mder1u/mp1**2+mder2udl/mp1)-2.0d0*mder1u/m
!   tem2=dl2*(-mder1dl*mder1u/m**2+mder2udl/m)-tem1dl*coef*mder1u/tem**2.0d0
!   tem3=(theta*mder1u+coef*mder2udl)/tem
!   lpdf2udl=tem1+tem2+tem3-theta*(t1/u1)/(t1-1.0d0)!correct
!   
!   tem1=th2*(-mder1dl*mder1v/mp1**2+mder2vdl/mp1)-2.0d0*mder1v/m
!   tem2=dl2*(-mder1dl*mder1v/m**2+mder2vdl/m)-tem1dl*coef*mder1v/tem**2.0d0
!   tem3=(theta*mder1v+coef*mder2vdl)/tem
!   lpdf2vdl=tem1+tem2+tem3-theta*(t2/u2)/(t2-1.0d0)!correct
!   
!   lderu(1)=lpdf1u;lderu(2)=lpdf2u;lderu(3)=lpdf2uth;lderu(4)=lpdf2udl;
!   lderv(1)=lpdf1v;lderv(2)=lpdf2v;lderv(3)=lpdf2vth;lderv(4)=lpdf2vdl;
!   
!   lderuv=lpdf2uv
!   return
! end
! 
! 
! ! Function with Frank lpdf, lder1, lder2 (partial wrt cpar and 2nd order)
! !  for factor 1.
! ! inputs
! !   u1,u2 = values in (0,1)
! !   cpar = scalar in (-oo,oo) but could be underflow problems for cpar>35
! ! outputs 
! !   lpdf = log pdf, 
! !   lder1 = \p lpdf/\p cpar
! !   lder2 = \p^2 lpdf/\p cpar^2
! subroutine lfrk1derivs(u1,u2,cpar,lpdf,lder1,lder2)
!   implicit none
!   double precision u1,u2,cpar,lpdf,lder1,lder2
!   double precision den,den1,den2,t0,t1,t2
!   t0 = exp(-cpar);
!   t1 = exp(-cpar*u1);
!   t2 = exp(-cpar*u2);
!   den = t1+t2-t0-t1*t2;
!   den1 = -u1*t1-u2*t2+t0+(u1+u2)*t1*t2;
!   den2 = u1*u1*t1+u2*u2*t2-t0-(u1+u2)*(u1+u2)*t1*t2;
!   lpdf = log(abs(cpar))+log(abs(1.d0-t0))-cpar*(u1+u2)-2.d0*log(abs(den));
!   lder1 = 1.d0/cpar+t0/(1.d0-t0)-(u1+u2)-2.d0*den1/den;
!   lder2 = -1.d0/(cpar*cpar)-t0/((1.d0-t0)*(1.d0-t0))-2.d0*den2/den+2.d0*den1*den1/(den*den); 
!   return
!   end
! 
! 
! ! outer product of two vectors
! ! input
! !   na = length of avec
! !   nb = length of bvec
! !   avec = vector 1
! !   bvec = vector 2
! ! output
! !   abmat = na x nb matrix with outer product of avec and bvec 
! subroutine outer(na,nb,avec,bvec,abmat)
!  implicit none
!  integer na,nb,i
!  double precision avec(na),bvec(nb),abmat(na,nb)
!   na=size(avec); nb=size(bvec)
!  do i =1,na  
!    abmat(i,:)=avec(i)*bvec
!  end do
!  return
!  end

  
