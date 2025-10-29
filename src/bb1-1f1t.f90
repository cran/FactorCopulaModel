
! This file can serve as a template for 1F1T
!  when the linking copula for the first factor has 2 parameters and
!  the linking copula for the second tree has 1 parameter.
!1F1T with BB1 for tree 1 rooted a latent
!   and Frank or Gaussian for tree 2
! return nllk and grad, hess

! test main program
!program lgauder
!  implicit none
!  double precision u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu, lder1v,lder2v,lder2uv,ldermixv;
!  read *, u1,u2,cpar
!  call lgau2derivt(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu, lder1v,lder2v,lder2uv,ldermixv)
!  print *, lpdf,lder1,lder2,lder1u,lder2u,ldermixu, lder1v,lder2v,lder2uv,ldermixv
!  stop
!  end

! test main program
!program lfrkder
!  implicit none
!  double precision u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu, lder1v,lder2v,lder2uv,ldermixv;
!  read *, u1,u2,cpar
!  call lfrk2derivt(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu, lder1v,lder2v,lder2uv,ldermixv)
!  print *, lpdf,lder1,lder2,lder1u,lder2u,ldermixu, lder1v,lder2v,lder2uv,ldermixv
!  stop
!  end

! test main program
!program bb1frkgau
!  implicit none
!  integer n,d,nq,i,iq,j,npar,ip, fam
!  integer, dimension(:), allocatable :: edg1,edg2
!  double precision, dimension(:,:), allocatable :: udata,hess
!  double precision, dimension(:), allocatable :: xl,wl,grad,param0
!  double precision nllk
!  read *,nq
!  allocate ( xl(nq),wl(nq) )
!  read *, xl(:)
!  read *, wl(:)
!  read *,n,d
!  npar=d*3-1
!  allocate ( udata(n,d), param0(npar), grad(npar), hess(npar,npar) )
!  allocate ( edg1(d-1), edg2(d-1) )
!  do i=1,n
!    read *, udata(i,:)
!  end do
!  read *, (param0(ip),ip=1,npar)
!  read *, edg1(:)
!  read *, edg2(:)
!  fam=2 ! for gaussian, fam=1 for frank
!  !fam=1 ! for frank, 
!  !call bb1frk1f1t(npar,param0,d,n,udata,nq,wl,xl, edg1,edg2, nllk,grad,hess)
!  call bb1m1f1t(fam,npar,param0,d,n,udata,nq,wl,xl, edg1,edg2, nllk,grad,hess)
!  print *, nllk
!  print "(8f10.5)", grad
!  print *, "grad above, hess below"
!  print *," "
!  print "(8f10.5)", hess
!  deallocate (xl, wl, param0, udata, grad, hess, edg1, edg2)
!  stop
!  end



! 1-factor 1-truncated model with BB1/Frank or BB1/Gaussian
! inputs 
!   fam =1 for Frank copula for residual dependence conditional on latent,
!   fam =2 for Gaussian copula for residual dependence conditional on latent,
!   npar = #parameters = 3*d-1
!   param0 = parameter vector (dimension d=npar), 
!       param0 has th1,dl1,th2,dl2, ... as vector theta>0, delta>1
!   d = #variables
!   n = sample size
!   udata = nxd matrix of uniform scores
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes
!   edg1 = node 1 vector of edges on residual tree  (length d-1)
!   edg2 = node 1 vector of edges on residual tree 
! outputs 
!   nllk = negative log-likelihood, 
!   grad = gradient of nllk, 
!   hess = hessian of nllk
subroutine bb1m1f1t(fam,npar,param0,d,n,udata,nq,wl,xl, edg1,edg2, nllk,grad,hess)
  implicit none
  integer fam,npar,d,n,nq, edg1(d-1),edg2(d-1)
  double precision param0(npar),udata(n,d),wl(nq),xl(nq)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,jj(2),jj2(2)
  integer d2,d1,ie,e1,e2,jj1(2)
  double precision, dimension(:), allocatable :: uvec,lpdf,fval1,integl1,grd
  double precision, dimension(:,:), allocatable :: cder1,cder2
  double precision, dimension(:), allocatable :: ccdf
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss,der1,der2
  double precision, dimension(:,:), allocatable :: param
  double precision, dimension(:), allocatable :: lgrad,partr
  double precision, dimension(:,:), allocatable :: lhess,gmat
  double precision fval,integl,ww,ljpdf,tmat(2,2)
  double precision lpdfr, ldercpar,ldercpar2,lder1u,lder1v,lder2u,lder2v,lder2uv, ldermixu,ldermixv

  ! npar=d*3-1
  allocate ( uvec(d), lpdf(d), der1(2,d), der2(3,d), fval1(npar), integl1(npar), grd(npar) )
  allocate ( fval2(npar,npar), integl2(npar,npar), hss(npar,npar), param(2,d) )
  allocate ( ccdf(d), cder1(2,d), cder2(3,d) )
  allocate (lgrad(npar), lhess(npar,npar), gmat(npar,npar), partr(d-1) )
  d2=2*d; d1=d-1
  param=reshape(param0(1:d2),(/2,d/))
  partr=param0((d2+1):npar)
  nllk=0.d0; grad=0.d0; hess=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    integl=0.d0; integl1=0.d0; integl2=0.d0; ! initialize for integrals
    do iq=1,nq ! loop over quadrature points
      lgrad=0.d0; lhess=0.d0; 
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        call lbb1derivs(xl(iq),uvec(j),param(:,j),lpdf(j),der1(:,j),der2(:,j))
        ! der2 with 3 cols : pth^2 pth.pdl, pdl^2
        call cbb1derivs(uvec(j),xl(iq),param(:,j),ccdf(j),cder1(:,j),cder2(:,j))
        !if(iq==1) print "(6f11.6)", lpdf(j),der1(:,j),der2(:,j)
      end do
      ! update integrand, and derivs wrt copula parameters 
      ljpdf=sum(lpdf)
      ! fval1: vector of partial deriv wrt param[j], j=1,...,np
      ! fval2: matrix of 2nd partial deriv wrt param[j], param[j2], 
      do j=1,d
        jj=(/2*j-1,2*j/)
        lgrad(jj)=lgrad(jj)+der1(:,j)
        tmat(1,1)=der2(1,j); tmat(2,2)=der2(3,j);
        tmat(1,2)=der2(2,j); tmat(2,1)=der2(2,j);
        lhess(jj,jj)=lhess(jj,jj)+tmat
      end do
      do ie=1,d1
        e1=edg1(ie); e2=edg2(ie)
        jj1=(/2*e1-1,2*e1/); jj2=(/2*e2-1,2*e2/);
        if (fam==1) then
          call lfrk2derivt(ccdf(e1),ccdf(e2),partr(ie), &
             lpdfr,ldercpar,ldercpar2,lder1u,lder2u,ldermixu, &
             lder1v,lder2v,lder2uv,ldermixv)
        else ! fam==2
          call lgau2derivt(ccdf(e1),ccdf(e2),partr(ie), &
             lpdfr,ldercpar,ldercpar2,lder1u,lder2u,ldermixu, &
             lder1v,lder2v,lder2uv,ldermixv)
        endif
        ljpdf=ljpdf+lpdfr
        lgrad(jj1)=lgrad(jj1)+lder1u*cder1(:,e1)
        lgrad(jj2)=lgrad(jj2)+lder1v*cder1(:,e2)
        call outer(2,2,cder1(:,e1),cder1(:,e1),tmat)
        lhess(jj1,jj1)=lhess(jj1,jj1)+lder2u*tmat
        call outer(2,2,cder1(:,e2),cder1(:,e2),tmat)
        lhess(jj2,jj2)=lhess(jj2,jj2)+lder2v*tmat
        call outer(2,2,cder1(:,e1),cder1(:,e2),tmat)
        lhess(jj1,jj2)=lhess(jj1,jj2)+lder2uv*tmat
        call outer(2,2,cder1(:,e2),cder1(:,e1),tmat)
        lhess(jj2,jj1)=lhess(jj2,jj1)+lder2uv*tmat
        tmat(1,1)=cder2(1,e1); tmat(2,2)=cder2(3,e1);
        tmat(1,2)=cder2(2,e1); tmat(2,1)=cder2(2,e1);
        lhess(jj1,jj1)=lhess(jj1,jj1)+lder1u*tmat
        tmat(1,1)=cder2(1,e2); tmat(2,2)=cder2(3,e2);
        tmat(1,2)=cder2(2,e2); tmat(2,1)=cder2(2,e2);
        lhess(jj2,jj2)=lhess(jj2,jj2)+lder1v*tmat
        lgrad(d2+ie)=lgrad(d2+ie)+ldercpar
        lhess(d2+ie,d2+ie)=lhess(d2+ie,d2+ie)+ldercpar2
        lhess(jj1,d2+ie)=lhess(jj1,d2+ie)+ldermixu*cder1(:,e1)
        lhess(d2+ie,jj1)=lhess(d2+ie,jj1)+ldermixu*cder1(:,e1)
        lhess(jj2,d2+ie)=lhess(jj2,d2+ie)+ldermixv*cder1(:,e2)
        lhess(d2+ie,jj2)=lhess(d2+ie,jj2)+ldermixv*cder1(:,e2)
      end do
      ! update quadrature loops
      fval=exp(ljpdf)
      fval1=fval*lgrad
      call outer(npar,npar,lgrad,lgrad,gmat)
      fval2=fval*(lhess+gmat)
      ww=wl(iq)
      integl=integl+fval*ww
      integl1=integl1+fval1*ww
      integl2=integl2+fval2*ww
    end do
    ! update contribution to negative log-likelihood
    if(integl<=0.d0) integl=1.d-100
    nllk=nllk-log(integl)
    grd=integl1/integl; grad=grad-grd;
    do j=1,npar
      do j2=1,npar
        hss(j,j2)=integl2(j,j2)/integl - grd(j)*grd(j2); 
      end do
    end do
    hess=hess-hss
  end do
  deallocate (uvec, lpdf, der1, der2, fval1, fval2, integl1, integl2, grd, hss)
  deallocate (ccdf, cder1, cder2)
  deallocate ( param, partr )
  deallocate ( lgrad,lhess,gmat )
  return
  end


!--------------------------------------------------------------------

! this routine is in bb1facts.f90
! Function with BB1 lpdf, lder1, lder2 (partial wrt param and 2nd order)
!  for factor 1;
! lder1 is 2-vector, lder2 is 2x2 matrix
! inputs
!   u1,u2 = values in (0,1)
!   param = (theta, delta), theta>0, delta>1
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p param
!   lder2 = \p^2 lpdf/\p param \p param^T
!subroutine lbb1derivs(u1,u2,param,lpdf,lder1,lder2)
!  implicit none
!  double precision u1,u2,param(2),lpdf,lder1(2),lder2(3)
!  double precision theta,delta,der1th,der1dl,der2th,der2dl,derthd
!  double precision t10,t1der1th,t1der1dl,t1der2th,t1der2dl,t1derthd
!  double precision t20,t2der1th,t2der1dl,t2der2th,t2der2dl,t2derthd
!  double precision t30,t3der1th,t3der1dl,t3der2th,t3der2dl,t3derthd
!  double precision t40,t4der1th,t4der1dl,t4der2th,t4der2dl,t4derthd
!  double precision m,mder1th,mder1dl,mder2th,mder2dl,mderthd
!  double precision den,coef,t1,t2,tu1,tu2
!  double precision mp1,msq,thsq,den2,thtem,dltem,dl1,t1a,t2a,smlog

!  theta=param(1); delta=param(2);

!  call mderivs(u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd)
!  mp1=1.d0+m; msq=m*m; thsq=theta*theta
!  thtem=2.d0+1.d0/theta
!  t10 = -(thtem)*log(mp1)
!  t1der1th = log(mp1)/thsq - (thtem)*mder1th/mp1
!  t1der1dl = -(thtem)*mder1dl/mp1
!  t1der2th = -2.d0*log(mp1)/theta**3+(2.d0/thsq)*mder1th/mp1
!  t1der2th = t1der2th-(thtem)*(mder2th/mp1-mder1th*mder1th/(mp1*mp1))
!  t1der2dl = -(thtem)*(mder2dl/mp1-mder1dl*mder1dl/(mp1*mp1))
!  t1derthd = mder1dl/(thsq*(mp1))-(thtem)*(mderthd/(mp1)-mder1th*mder1dl/(mp1*mp1))

!  dltem=1.d0-2.d0*delta; dl1=delta-1.d0
!  t20 = (dltem)*log(m)
!  t2der1th = (dltem)*mder1th/m
!  t2der1dl = -2.d0*log(m)+(dltem)*mder1dl/m
!  t2der2th = (dltem)*(mder2th/m-mder1th*mder1th/msq)
!  t2der2dl = -4.d0*mder1dl/m+(dltem)*(mder2dl/m-mder1dl*mder1dl/msq)
!  t2derthd = -2.d0*mder1th/m+(dltem)*(mderthd/m-mder1th*mder1dl/msq)

!  coef = theta*delta+1.d0
!  den = theta*(dl1)+(coef)*m
!  den2=den*den
!  t30 = log(theta*(dl1)+coef*m)
!  t3der1th = (dl1+delta*m+coef*mder1th)/den
!  t3der1dl = (theta+theta*m+coef*mder1dl)/den
!  t3der2th = (2.d0*delta*mder1th+coef*mder2th)/den
!  t3der2th = t3der2th-(dl1+delta*m+coef*mder1th)**2/den2
!  t3der2dl = (2.d0*theta*mder1dl+coef*mder2dl)/den
!  t3der2dl = t3der2dl-(theta+theta*m+coef*mder1dl)**2/den2;
!  t3derthd = (1.d0+m+theta*mder1th+delta*mder1dl+coef*mderthd)/den;
!  t3derthd = t3derthd - (theta+theta*m+coef*mder1dl)*(dl1+delta*m+coef*mder1th)/den2

!  t1 = u1**(-theta); t2 = u2**(-theta)
!  t1a=t1-1.d0; t2a=t2-1.d0; smlog=log(t1a)+log(t2a)
!  tu1 = -log(u1); tu2 = -log(u2)
!  t40 = (dl1)*(smlog)+(theta+1.d0)*(tu1+tu2)
!  t4der1th = (dl1)*(t1*tu1/t1a+t2*tu2/t2a)+tu1+tu2  
!  t4der1dl = smlog
!  t4der2th = -(dl1)*(t1*tu1*tu1/t1a**2+t2*tu2*tu2/t2a**2)
!  t4der2dl = 0.d0
!  t4derthd = t1*tu1/t1a+t2*tu2/t2a

!  lpdf = t10+t20+t30+t40
!  der1th = t1der1th+t2der1th+t3der1th+t4der1th
!  der1dl = t1der1dl+t2der1dl+t3der1dl+t4der1dl
!  der2th = t1der2th+t2der2th+t3der2th+t4der2th
!  der2dl = t1der2dl+t2der2dl+t3der2dl+t4der2dl
!  derthd = t1derthd+t2derthd+t3derthd+t4derthd
!  lder1(1)=der1th; lder1(2)=der1dl  ! gradient of log pdf
!  lder2(1)=der2th; lder2(2)=derthd; lder2(3)=der2dl; ! Hessian terms
!  return
!  end

! this routine is in bb1facts.f90
! Function with BB1 condcdf C_{1|2}(u1|u2;theta,delta) for factor 1, and
! ccdf cder1(2), cder2(3) (partial wrt param, 1st and 2nd order),
! second deriv has theta^2, mix, delta^2
! inputs
!   u1,u2 = values in (0,1)
!   param = (theta, delta), theta>0, delta>1
! outputs
!   ccdf = conditional cdf, 
!   cder1 = \p ccdf/\p param
!   cder2 = \p^2 ccdf/\p param \p param^T
!subroutine cbb1derivs(u1,u2,param,ccdf,cder1,cder2)
!  implicit none
!  double precision u1,u2,param(2),ccdf,cder1(2),cder2(3)
!  double precision theta,delta
!  double precision m,mder1th,mder1dl,mder2th,mder2dl,mderthd
!  double precision t2,tu2,cf1,logm
!  double precision mp1,msq,thsq,dl1n,t2a,lt2a,lmp1
!  double precision lcdf,lder1th,lder1dl,lder2th,lder2dl,lderthd
!  double precision cder1th,cder1dl,cder2th,cder2dl,cderthd

!  theta=param(1); delta=param(2)
!
!  call mderivs(u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd)
!  mp1=1.d0+m; msq=m*m; thsq=theta*theta
!  lmp1=log(mp1); logm=log(m);

!  t2 = u2**(-theta); tu2 = -log(u2);
!  t2a=t2-1.d0; lt2a=log(t2a)
!  cf1 = 1.d0+1.d0/theta;
!  dl1n=1.d0-delta

!  lcdf = -cf1*lmp1+(dl1n)*(logm-lt2a)+(theta+1.d0)*tu2;
!  lder1th = lmp1/thsq-cf1*mder1th/mp1;
!  lder1th = lder1th+(dl1n)*(mder1th/m-t2*tu2/t2a)+tu2;
!  lder1dl = -cf1*mder1dl/mp1+(dl1n)*mder1dl/m-logm+lt2a;
!  lder2th = -2.d0*lmp1/(thsq*theta)+(2.d0/thsq)*mder1th/mp1-cf1*(mder2th/mp1-mder1th**2/(mp1*mp1));
!  lder2th = lder2th+(dl1n)*(mder2th/m-mder1th**2/msq+tu2*tu2*t2/(t2a*t2a));
!  lder2dl = -cf1*(mder2dl/mp1-mder1dl**2/(mp1*mp1))-2.d0*mder1dl/m;
!  lder2dl = lder2dl+(dl1n)*(mder2dl/m-mder1dl**2/msq);
!  lderthd = (1.d0/thsq)*mder1dl/mp1-cf1*(mderthd/mp1-mder1th*mder1dl/(mp1*mp1))-mder1th/m;
!  lderthd = lderthd+(dl1n)*(mderthd/m-mder1th*mder1dl/msq)+t2*tu2/t2a;

!  ccdf = exp(lcdf);
!  cder1th=ccdf*lder1th;
!  cder1dl=ccdf*lder1dl;
!  cder2dl=ccdf*(lder2dl+lder1dl**2);
!  cder2th=ccdf*(lder2th+lder1th**2);
!  cderthd=ccdf*(lderthd+lder1dl*lder1th);
!  cder1(1)=cder1th; cder1(2)=cder1dl  ! gradient of log pdf
!  cder2(1)=cder2th; cder2(2)=cderthd; cder2(3)=cder2dl; ! Hessian terms
!  return
!  end

! this routine is in bb1facts.f90
! 1st and 2nd order derivatives of m = s^(1/delta) wrt theta and delta
! inputs
!   u1,u2 = values in (0,1)
!   theta >0
!   delta >1
! outputs
!   m = s^(1/delta); s =  [u1^(-theta)-1]^delta + [u2^(-theta)-1]^delta
!   mder1th = \p m/\p theta
!   mder1dl = \p m/\p delta
!   mder2th = \p^2 m/\p theta^2
!   mder2dl = \p^2 m/\p delta^2
!   mderthd = \p^2 m/\p theta \p delta
!subroutine mderivs(u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd)
!  implicit none
!  double precision u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd
!  double precision t1,t2,tu1,tu2,ttu1,ttu2,td01,td02,td11,td12,td21,td22
!  double precision s,sder1th,sder1dl,sder2th,sder2dl,sderthd,m1,ts,m1der1dl 
!  double precision t1a,t2a,dlsq,dlcu
!
!  t1 = u1**(-theta); t2 = u2**(-theta)
!  t1a=t1-1.d0; t2a=t2-1.d0
!  tu1 = -log(u1); tu2 = -log(u2)
!  ttu1 = log(t1a); ttu2 = log(t2a)
!  td01 = (t1a)**delta; td02 = (t2a)**delta
!  td11=td01/t1a; td12=td02/t2a;
!  td21=td11/t1a; td22=td12/t2a;
 
!  s = td01+td02
!  sder1th = delta*(td11*t1*tu1+td12*t2*tu2)
!  sder1dl = td01*ttu1+td02*ttu2
!  sder2th = delta*(delta-1.d0)*(td21*t1*t1*tu1*tu1+td22*t2*t2*tu2*tu2)
!  sder2th = sder2th+delta*(td11*t1*tu1*tu1+td12*t2*tu2*tu2)
!  sder2dl = td01*ttu1*ttu1+td02*ttu2*ttu2
!  sderthd = sder1th/delta+delta*(td11*ttu1*tu1*t1+td12*ttu2*tu2*t2)

!  m = s**(1.d0/delta); m1 = m/s
!  ts = log(s)
!  dlsq=delta*delta; dlcu=delta*dlsq
!  mder1th = m1*sder1th/delta
!  mder1dl = m1*sder1dl/delta - m*ts/dlsq
!  m1der1dl = mder1dl/s - m*sder1dl/s**2
!  mder2th = (1.d0-delta)*m1*sder1th**2/(dlsq*s)+m1*sder2th/delta
!  mder2dl = 2.d0*m*ts/dlcu-mder1dl*ts/dlsq-2.d0*m1*sder1dl/dlsq
!  mder2dl = mder2dl+sder2dl*m1/delta+sder1dl*m1der1dl/delta
!  mderthd = -m1*sder1th/dlsq+sder1th*m1der1dl/delta+m1*sderthd/delta
!  return
!  end

! derivatives of Frank log copula density
! March 2017: addition of deriv wrt v as well as deriv wrt u

! in 2derivs.f90
! Function with Frank log pdf c(u1,u2,cpar) for tree 2, and 
!   lder1, lder2 (partial wrt cpar, 1st and 2nd order)
!   also lder1u lder2u (partial wrt u1, 1st and 2nd order)
! ***  also lder1v (partial wrt u2, 1st order) and others
!   and ldermixu  (partial wrt cpar and u1) 
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (-oo,oo)
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p cpar
!   lder2 = \p^2 lpdf/\p cpar^2
!   lder1u = \p lpdf/\p u1
!   lder2u = \p^2 lpdf/\p u1^2
!   ldermixu = \p^2 lpdf/\p u1 \p cpar
!   lder1v = \p lpdf/\p u2
!   lder2v = \p^2 lpdf/\p u2^2
!   lder2uv = \p^2 lpdf/\p u1 \p u2
!   ldermixv = \p^2 lpdf/\p u2 \p cpar
!  some derivs added; function name with ending of 't'
!subroutine lfrk2derivt(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu, lder1v,lder2v,lder2uv,ldermixv)
!  implicit none
!  double precision u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu
!  double precision lder1v,lder2v,lder2uv,ldermixv
!  double precision den,den1,den2,den1u,den2u,denmixu,t0,t1,t2
!  double precision den1v,den2v,denmixv
!
!  t0 = exp(-cpar);
!  t1 = exp(-cpar*u1);
!  t2 = exp(-cpar*u2);
!  den = t1+t2-t0-t1*t2;
!  den1 = -u1*t1-u2*t2+t0+(u1+u2)*t1*t2;
!  den2 = u1*u1*t1+u2*u2*t2-t0-(u1+u2)*(u1+u2)*t1*t2;
!  den1u = -cpar*t1*(1.d0-t2);  
!  den1v = -cpar*t2*(1.d0-t1);   ! added
!  den2u = cpar*cpar*t1*(1.d0-t2);
!  den2v = cpar*cpar*t2*(1.d0-t1); ! added
!  denmixu = t1*(-1.d0+cpar*u1+t2-(u1+u2)*cpar*t2); 
!  denmixv = t2*(-1.d0+cpar*u2+t1-(u1+u2)*cpar*t1); ! added 
!   
!  !lpdf = log(abs(cpar))+log(abs(1-t0))-cpar*(u1+u2)-2*log(abs(den));
!  ! 121023 maybe later add the limits as cpar->0
!  ! pdf = cpar*(1-t0)/den^2 where
!  !    1-t0 has same sign as cpar,  den has same sign as cpar
!  lpdf = log(cpar*(1.d0-t0))-cpar*(u1+u2)-2.d0*log(abs(den));
!  lder1 = 1.d0/cpar+t0/(1.d0-t0)-(u1+u2)-2.d0*den1/den;
!  lder2 = -1.d0/(cpar*cpar)-t0/((1.d0-t0)*(1.d0-t0))-2.d0*den2/den+2.d0*den1*den1/(den*den); 
!  lder1u = -cpar-2.d0*den1u/den;
!  lder1v = -cpar-2.d0*den1v/den;  ! added
!  lder2u = -2.d0*den2u/den+2.d0*den1u*den1u/(den*den);
!  lder2v = -2.d0*den2v/den+2.d0*den1v*den1v/(den*den); ! added
!  lder2uv = 2.d0*cpar*cpar*t1*t2/den+2.d0*den1u*den1v/(den*den); ! added
!  ldermixu = -1.d0-2.d0*denmixu/den+2.d0*den1u*den1/(den*den);
!  ldermixv = -1.d0-2.d0*denmixv/den+2.d0*den1v*den1/(den*den); ! added
!  return
!  end


! Function with bivariate Gaussian copula log pdf c(u1,u2,cpar) for tree 2, and 
!   lder1, lder2 (partial wrt cpar, 1st and 2nd order)
!   also lder1u lder2u (partial wrt u1, 1st and 2nd order)
! ***  also lder1v (partial wrt u2, 1st order) and others
!   and ldermixu  (partial wrt cpar and u1) 
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (-oo,oo)
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p cpar
!   lder2 = \p^2 lpdf/\p cpar^2
!   lder1u = \p lpdf/\p u1
!   lder2u = \p^2 lpdf/\p u1^2
!   ldermixu = \p^2 lpdf/\p u1 \p cpar
!   lder1v = \p lpdf/\p u2
!   lder2v = \p^2 lpdf/\p u2^2
!   lder2uv = \p^2 lpdf/\p u1 \p u2
!   ldermixv = \p^2 lpdf/\p u2 \p cpar
!  some derivs added; function name with ending of 't'
!subroutine lgau2derivt(u1,u2,rho,lpdf,lder1,lder2,lder1u,lder2u,ldermixu, lder1v,lder2v,lder2uv,ldermixv)
!  implicit none
!  double precision u1,u2,rho,lpdf,lder1,lder2,lder1u,lder2u,ldermixu
!  double precision lder1v,lder2v,lder2uv,ldermixv
!  double precision qnorms, dnorms
!  double precision x,y,x2,y2,rh2,den,den2,con,qf1,qf2,dx,dy,dx2,dy2
!
!  x=qnorms(u1); x2=x*x
!  y=qnorms(u2); y2=y*y
!  rh2=rho**2; den=1.d0-rh2; den2=den*den
!  con= -0.5d0*log(den)
!  qf1= -(x2+y2-2.d0*rho*x*y)/(2*den)
!  qf2=(x2+y2)/2.d0
!  lpdf=con+qf1+qf2
!  lder1 = rho/den+2.d0*qf1*rho/den + x*y/den
!  dx=dnorms(x); dy=dnorms(y); dx2=dx*dx; dy2=dy*dy
!  lder1u = (x-(x-rho*y)/den)/dx
!  lder1v = (y-(y-rho*x)/den)/dy
!  lder2u = (rho*(y-rho*x)/den)*x/dx2 - rh2/den/dx2
!  lder2v = (rho*(x-rho*y)/den)*y/dy2 - rh2/den/dy2
!  lder2uv = rho/den/dx/dy
!  lder2=(1.d0+rh2)/den2 + (qf2*2.d0*(-1.d0-3.d0*rh2)+2.d0*rho*(3.d0+rh2)*x*y)/den2/den
!  ldermixu=((1.d0+rh2)*y-2.d0*rho*x)/den2/dx
!  ldermixv=((1.d0+rh2)*x-2.d0*rho*y)/den2/dy
!  return
!  end

