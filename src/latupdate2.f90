
! !main function 
! PROGRAM test
!   implicit none
!   integer nq,n,d,i
!   double precision,dimension(:,:),allocatable::udata
!   integer,dimension(:),allocatable::family
!   double precision,dimension(:),allocatable::xl,wl,lat,th

!   read*,nq
!   print*,nq
!   allocate(xl(nq),wl(nq))
!   read*,xl
!   read*,wl
  
!   print*,xl
!   print*,wl
  
!   read*,n
!   read*,d
!   allocate(udata(n,d),th(2*d),lat(n),family(d))
!   read *, (udata(i,:), i=1,n) 
  
!   read*,th
!   print*,th
!   read*,family
!   print*,family
  
!   call latupdate3(th,n,d,udata,nq,xl,wl,family,lat)
  
!   print*,lat
  
!   deallocate(udata,xl,wl,family)
! end PROGRAM



!latent update in one-factor copula model; different copula families are allowed

subroutine latupdate3 (th,n,d,udata,nq,xl,wl,family,lat)
  implicit none
  integer n,nq,i,iq,j,d,jj(2),jj2
  double precision num,dem,lpdf,rho,nu,t1,t2
  integer family(d)
  double precision lder11(2),lder22(3),lderu(4),lderv(4),lderuv,lder1,lder2,ltder1u,ltder2u,ltdermix
  double precision udata(n,d),th(2*d),wl(nq),xl(nq),uvec(d),val(nq),tem(nq),lat(n)
  double precision ut,vt,qt
  do i=1,n
    uvec=udata(i,:)
    val=0.0d0;
    do j=1,d
      jj=(/2*j-1,2*j/)
      jj2=2*j-1
  ! HJ Feb 2025, added family=14 for survival Gumbel
     do iq=1,nq
      if (family(j)==4) then
        call lgum(uvec(j),xl(iq),th(jj2),lpdf)
      elseif (family(j)==14) then
        ut=1.0d0-uvec(j)
        vt=1.0d0-xl(iq)
        call lgum(ut,vt,th(jj2),lpdf)
      elseif(family(j)==5) then
        call lfrk(uvec(j),xl(iq),th(jj2),lpdf)
      elseif(family(j)==1)then
        call lgau(uvec(j),xl(iq),th(jj2),lpdf)
      elseif(family(j)==2) then
        rho=th(2*j-1);nu=th(2*j);
        t1=qt(uvec(j),nu);t2=qt(xl(iq),nu);
        call lt2derivs(t1,t2,rho,nu,lpdf,lder1,lder2,ltder1u,ltder2u,ltdermix)
      elseif(family(j)==7) then
        call lbb1derivs(uvec(j),xl(iq),th(jj),lpdf,lder11,lder22,lderu,lderv,lderuv)
      elseif(family(j)==17) then
        ut=1.0d0-uvec(j)
        vt=1.0d0-xl(iq)
        call lbb1derivs(ut,vt,th(jj),lpdf,lder11,lder22,lderu,lderv,lderuv)
      elseif(family(j)==10) then
        call lbb8derivs(uvec(j),xl(iq),th(jj),lpdf,lder11,lder22,lderu,lderv,lderuv)
      elseif(family(j)==20) then
        ut=1.0d0-uvec(j)
        vt=1.0d0-xl(iq)
        call lbb8derivs(ut,vt,th(jj),lpdf,lder11,lder22,lderu,lderv,lderuv)
     end if
      val(iq)=val(iq)+lpdf
      end do
    end do
    tem=exp(val)
    num=dot_product(wl,tem*xl)
    dem=dot_product(wl,tem)
    lat(i)=num/dem  
  end do   
return
end subroutine





! !this functions returns the log-density of Gumbel copula
! !inputs: u1,u2,cpar:paramter in the copula 
! !output: lpdf: log density of Gumbel copula 
! subroutine lgum(u1,u2,cpar,lpdf)
!   implicit none
!   double precision u1,u2,cpar,lpdf
!   double precision x,y,tx,ty,logm,den,m,s,xd,yd
  
!   x = -log(u1); y = -log(u2);
!   tx = log(x); ty = log(y);
!   xd = x**cpar; yd = y**cpar; s = xd+yd; m = s**(1.d0/cpar); 
!   den = m+cpar-1.d0
  
!   logm=log(m);
!   lpdf = -m+log(den)+(1.d0-2.d0*cpar)*logm+(cpar-1.d0)*(tx+ty)+x+y;
  
! end subroutine



! this routine is in bb1facts.f90
!Function with bb8 logpdf =log f(u1,u2,th,delta) 
!lder1, lder2 (partial wrt cparv, 1st and 2nd order) lder1 is 2-dimenisonal vector
!lder2 is 
 !also lderu (partial wrt u1, 1st and 2nd order and also wrt parmaters)
! ***  also lderv (partial wrt u2, 1st order and 2nd order and also wrt parameters) and others
! and lderuv  (partial wrt u and v) 
! outputs 
!lpdf = log pdf, 
!lder11 = \p lpdf/\p nu  \p lpdf/\p delta !2-dimensional vector
!lder22 = \p^2 lpdf/\p^2 th^2, \p^2 lpdf/\p^2 th delta, lpdf/\p^2 delta^2!  3d-vector
!lderu = \p lpdf/\p u, \p^2 lpdf/\p u^2, \p^2 lpdf/\p u th, \p^2 lpdf/\p u delta  4d-vector
!lderv = \p lpdf/\p v, \p^2 lpdf/\p v^2, \p^2 lpdf/\p v th, \p^2 lpdf/\p v delta  4d-vector
!lderuv= \p^2 lpdf/\p uv
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

! !This functions returns the log-density of Frank copula
! !inputs: u1,u2,cpar:paramter in the copula 
! !output: lpdf: log density of Frank copula 
! ! subroutine lfrk(u1,u2,cpar,lpdf)
! !   implicit none
! !   double precision u1,u2,cpar,lpdf
! !   double precision den,t0,t1,t2
! ! 
! !   t0 = exp(-cpar);
! !   t1 = exp(-cpar*u1);
! !   t2 = exp(-cpar*u2);
! !   den = t1+t2-t0-t1*t2;
! !   lpdf = log(cpar*(1.d0-t0))-cpar*(u1+u2)-2.d0*log(abs(den));
! !   
! ! end subroutine
! ! 
! ! subroutine lgau(u1,u2,rho,lpdf)
! !   implicit none
! !   double precision u1,u2,rho,lpdf
! !   double precision qnorms
! !   double precision x,y,x2,y2,rh2,con,den,den2,qf1,qf2
! ! 
! !   x=qnorms(u1); x2=x*x
! !   y=qnorms(u2); y2=y*y
! !   rh2=rho**2; den=1.d0-rh2; den2=den*den
! !   con= -0.5d0*log(den)
! !   qf1= -(x2+y2-2.d0*rho*x*y)/(2*den)
! !   qf2=(x2+y2)/2.d0
! !   lpdf=con+qf1+qf2
! !   return
! !   end
