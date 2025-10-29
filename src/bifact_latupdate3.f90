! PROGRAM test2   ! main function
! implicit none
!   integer npar,mgrp,dvar,n,i,nq,iq
!   double precision v0,int11,int22
!   double precision,dimension(:,:),allocatable::udata
!   double precision,dimension(:),allocatable::th
!   double precision,dimension(:),allocatable::vg,uvec
!   integer,dimension(:),allocatable::grsize
!   integer,dimension(:),allocatable::family
!   double precision,dimension(:),allocatable::xl,wl,v0mat
!   double precision,dimension(:,:),allocatable::vgmat

!     read *,nq
!     allocate(xl(nq),wl(nq))
!     read*,xl
!     read*,wl
!     ! do iq =1,nq   ! equidistant and equally weighted for testing
!     !       xl(iq)=iq/(nq+1.d0)
!     !       wl(iq)=1.d0/nq
!     !   end do
!     read *,npar
    
!     read *,dvar
    
!     read *,mgrp
   
!     allocate(grsize(mgrp))
!     read *,grsize
    
!     read*,n
    
!     allocate(th(npar))
!     allocate(family(dvar*2))
!     read*,family

   

!     read *,th
    
!     allocate(udata(n,dvar))
!     read *, (udata(i,:), i=1,n) 
!     !print*,(udata(i,:), i=1,n)
!     allocate(vg(mgrp))
!     !allocate(v0mat(n),vgmat(n,mgrp))
!     v0=0.9977785
!     uvec=udata(1,:)
!     call innerint3(npar,th,mgrp,dvar,family,grsize,uvec,nq,wl,xl,v0,vg,int11,int22)
!     ! print*,int11
!     ! print*,int22
!     ! call latupdatebifact3 (npar,th,mgrp,family,dvar,n,grsize,udata,nq,xl,wl,v0mat,vgmat)
!   !   !call latupdatebifact3 (npar,th,mgrp,dvar,n,grsize,udata,nq,xl,wl,v0mat,vgmat)
!    !   print *, "the estimats for global latent is "
!    !   print *,v0mat
!    !   print *,"the estimates for local latent is"
!    ! do i=1,n
!    !       print *, vgmat(i,:)
!    !   end do
!     deallocate(grsize,th,udata,xl,wl,vg)
! END PROGRAM


 !compute the inner integral for all the combinations of copulas
 subroutine innerint3(npar,th,mgrp,dvar,family,grsize,uvec,nq,wl,xl,v0,vg,int11,int22)
  implicit none
  integer npar,mgrp,dvar,nq,i,ind
  double precision th(npar),wl(nq),xl(nq)
  double precision lder11(2), lder22(3), lder1u, lder2u,ldermixuu(2), lder1v, lder2v
  double precision lder2uv, ldermixvv(2), cder11(2),cder22(3)
  !double precision lderu(4), lderv(4), lderuv, t1, t2, nu, rho
  double precision llpdf(2*dvar),ccdf(dvar),uvec(dvar)
  double precision llk,llk1,llk2,v0,vg(mgrp),int11,int22
  integer iq2,jg,mj,ind1,ind2,grsize(mgrp),family(2*dvar)
  double precision intj(mgrp), intj2(mgrp),int1(mgrp),int2(mgrp)

 
   ind = 0; intj = 0.d0; intj2=0.d0; int1=0.d0; int2=0.d0;
    do jg =1,mgrp   !jth group                         
         ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
         do iq2 =1,nq   
           llk1 = 0.d0; llk2 = 0.d0   
           do mj = ind1,ind2 ! within group index
               call ccopderiv(uvec(mj),v0,th((/mj,mj+dvar/)),family(mj),&
              ccdf(mj),cder11,cder22)
               if(ccdf(mj) < 0.00001) ccdf(mj) = 0.00001
               if(ccdf(mj) > 0.99999) ccdf(mj) = 0.99999   
               
               call lcop2derivt(uvec(mj),v0,th((/mj,mj+dvar/)),family(mj),llpdf(mj),lder11,&
                  lder22,lder1u,lder2u,ldermixuu,lder1v,lder2v,lder2uv,ldermixvv)
              
             
               call lcop2derivt(ccdf(mj),xl(iq2),th((/2*dvar+mj,mj+3*dvar/)),family(dvar+mj),llpdf(dvar+mj),lder11,&
                  lder22,lder1u,lder2u,ldermixuu,lder1v,lder2v,lder2uv,ldermixvv) !c_{i,V_j;V_0} 


            llk1 = llk1+llpdf(mj+dvar)
            llk2 = llk2+llpdf(mj)

           end do

           

           int1(jg)=int1(jg)+wl(iq2)*exp(llk1)
           int2(jg)=exp(llk2)
      

           llk=llk1+llk2
           intj(jg)=intj(jg)+wl(iq2)*exp(llk)*xl(iq2)
           intj2(jg)=intj2(jg)+wl(iq2)*exp(llk)

           
         end do   
       end do


       do i=1,mgrp
        vg(i)=intj(i)/intj2(i)
       end do


       int11=product(int1)
       int22=product(int2)
    return
   end  
! 

subroutine latupdatebifact3 (npar,th,mgrp,family,dvar,n,grsize,udata,nq,xl,wl,v0mat,vgmat)
  implicit none
  integer npar,mgrp,n,nq,iq,i,dvar
  integer grsize(mgrp),family(2*dvar)
  double precision th(npar),udata(n,dvar),wl(nq),xl(nq),uvec(dvar)
  double precision v0,vg(mgrp),int11,int22,int1vec(nq),int2vec(nq)
  double precision num,dem,tem(nq),tem1(nq),v0mat(n),vgmat(n,mgrp)
  do i=1,n
     uvec = udata(i,:)   
     int1vec=0.d0;int2vec=0.d0;
     do iq =1,nq
      v0=xl(iq)
      call innerint3(npar,th,mgrp,dvar,family,grsize,uvec,nq,wl,xl,v0,vg,int11,int22)
      int1vec(iq)=int11;int2vec(iq)=int22
      tem1(iq)=int11*int22
      tem(iq)=xl(iq)*wl(iq)
      end do
     num=dot_product(tem,tem1)
     dem=dot_product(tem1,wl)
     v0mat(i)=num/dem
     ! print*,num
     ! print*,dem
     
     
     v0=v0mat(i)
    
     call innerint3(npar,th,mgrp,dvar,family,grsize,uvec,nq,wl,xl,v0,vg,int11,int22)
     vgmat(i,:)=vg
   end do
 return
end subroutine



! ! Function with Gumbel condcdf C_{1|2}(u1|u2;cpar) for factor 1, and 
! ! ccdf cder1, cder2 (partial wrt cpar, 1st and 2nd order)
! ! inputs
! !   u1,u2 = values in (0,1)
! !   cpar = scalar in (1,oo)
! ! outputs
! !   ccdf = conditional cdf, 
! !   cder1 = \p ccdf/\p cpar
! !   cder2 = \p^2 ccdf/\p cpar^2
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



! !   Function with Frank condcdf C_{1|2}(u1|u2;cpar) for factor 1, and
! ! ccdf cder1, cder2 (partial wrt cpar, 1st and 2nd order)
! ! inputs
! !   u1,u2 = values in (0,1)
! !   cpar = scalar in (-oo,oo)
! ! outputs
! !   ccdf = conditional cdf, 
! !   cder1 = \p ccdf/\p cpar
! !   cder2 = \p^2 ccdf/\p cpar^2
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

! !this routine is in bb1facts.f90
! ! Function with BB1 condcdf C_{1|2}(u1|u2;theta,delta) for factor 1, and
! ! ccdf cder11(2), cder22(3) (partial wrt param, 1st and 2nd order),
! ! second deriv has theta^2, mix, delta^2
! ! inputs
! !   u1,u2 = values in (0,1)
! !   param = (theta, delta), theta>0, delta>1
! ! outputs
! !   ccdf = conditional cdf, 
! !   cder11 = \p ccdf/\p param
! !   cder22= \p^2 ccdf/\p param \p param^T

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

!   call mderivs2(u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd)
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
 
!   eta2nudel=nu*log(1-delta)*((1.0d0-delta)**(nu-1.0d0))+(1.0d0-delta)**(nu-1.0d0)
!   x2nu=-(log(1.0d0-delta*u1)**2)*((1.0d0-delta*u1)**nu)!\p^2 x/\p nu^2
!   y2nu=-(log(1.0d0-delta*u2)**2)*((1.0d0-delta*u2)**nu)!\p^2 y/\p nu^2
!   x2del=-u1*u1*nu*(nu-1.0d0)*((1-delta*u1)**(nu-2.0d0))!\p^2 x/\p delta^2
!   y2del=-u2*u2*nu*(nu-1.0d0)*((1-delta*u2)**(nu-2.0d0))!\p^2 y/\p delta^2

!   x2nudel=u1*nu*(1-delta*u1)**(nu-1.0d0)*log(1-delta*u1)+u1*((1-delta*u1)**nu)/(1-delta*u1)

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
! ! 
! ! 
! ! this routine is in bb1facts.f90
! ! Function with bb8 logpdf =log f(u1,u2,th,delta) 
! ! lder1, lder2 (partial wrt cparv, 1st and 2nd order) lder1 is 2-dimenisonal vector
! ! lder2 is 
! !  also lderu (partial wrt u1, 1st and 2nd order and also wrt parmaters)
! ! ***  also lderv (partial wrt u2, 1st order and 2nd order and also wrt parameters) and others
! ! and lderuv  (partial wrt u and v) 
! ! outputs 
! ! lpdf = log pdf, 
! ! lder11 = \p lpdf/\p nu  \p lpdf/\p delta !2-dimensional vector
! ! lder22 = \p^2 lpdf/\p^2 th^2, \p^2 lpdf/\p^2 th delta, lpdf/\p^2 delta^2!  3d-vector
! ! lderu = \p lpdf/\p u, \p^2 lpdf/\p u^2, \p^2 lpdf/\p u th, \p^2 lpdf/\p u delta  4d-vector
! ! lderv = \p lpdf/\p v, \p^2 lpdf/\p v^2, \p^2 lpdf/\p v th, \p^2 lpdf/\p v delta  4d-vector
! ! lderuv= \p^2 lpdf/\p uv
! ! subroutine lbb1derivs(u1,u2,cparv,lpdf,lder11,lder22,lderu,lderv,lderuv)
! !   implicit none
! !   double precision u1,u2,cparv(2),lpdf,lder11(2),lder22(3)
! !   double precision theta,delta,der1th,der1dl,der2th,der2dl,derthd
! !   double precision t10,t1der1th,t1der1dl,t1der2th,t1der2dl,t1derthd
! !   double precision t20,t2der1th,t2der1dl,t2der2th,t2der2dl,t2derthd
! !   double precision t30,t3der1th,t3der1dl,t3der2th,t3der2dl,t3derthd
! !   double precision t40,t4der1th,t4der1dl,t4der2th,t4der2dl,t4derthd
! !   double precision m,mder1th,mder1dl,mder2th,mder2dl,mderthd,mderu(4),mderv(4)
! !   double precision den,coef,t1,t2,tu1,tu2
! !   double precision mp1,msq,thsq,den2,thtem,dltem,dl1,t1a,t2a,smlog
! !   double precision mder1u,mder2u,mder2uth,mder2udl,mder1v,mder2v,mder2vth,mder2vdl,mder2uv
! !   double precision lderu(4),lderv(4),lderuv,tem,tem1u,tem1v,lpdf1u,lpdf2u,lpdf2uth,lpdf2udl
! !   double precision lpdf2uv,lpdf1v,lpdf2v,lpdf2vth,lpdf2vdl,dl2,th1
! !   double precision tem1,tem2,tem3,th2,tem1th,tem1dl,tem2uth,tem2vth
! !   theta=cparv(1); delta=cparv(2);
! ! 
! !   
! !   call  mderivs(u1,u2,cparv,m,mder1th,mder1dl,mder2th,mder2dl,mderthd,mderu,mderv,mder2uv)
! !   print*,mderu
! !   print*,mderv
! !   mder1u=mderu(1);mder2u=mderu(2);mder2uth=mderu(3);mder2udl=mderu(4)
! !   mder1v=mderv(1);mder2v=mderv(2);mder2vth=mderv(3);mder2vdl=mderv(4);
! !   mp1=1.d0+m; msq=m*m; thsq=theta*theta
! !   thtem=2.d0+1.d0/theta
! !   t10 = -(thtem)*log(mp1)
! !   t1der1th = log(mp1)/thsq - (thtem)*mder1th/mp1
! !   t1der1dl = -(thtem)*mder1dl/mp1
! !   t1der2th = -2.d0*log(mp1)/theta**3+(2.d0/thsq)*mder1th/mp1
! !   t1der2th = t1der2th-(thtem)*(mder2th/mp1-mder1th*mder1th/(mp1*mp1))
! !   t1der2dl = -(thtem)*(mder2dl/mp1-mder1dl*mder1dl/(mp1*mp1))
! !   t1derthd = mder1dl/(thsq*(mp1))-(thtem)*(mderthd/(mp1)-mder1th*mder1dl/(mp1*mp1))
! !   
! !   tem=theta*(delta-1.0d0)+(theta*delta+1.d0)*m
! !   tem1u=(theta*delta+1.0d0)*mder1u
! !   tem1v=(theta*delta+1.0d0)*mder1v
! !   tem1th=delta-1.0d0+mder1th+delta*(m+mder1th*theta)
! !   tem1dl=theta+mder1dl+theta*(m+mder1dl*delta)
! !   tem2uth=delta*mder1u+mder2uth*(theta*delta+1.0d0)
! !   tem2vth=delta*mder1v+mder2vth*(theta*delta+1.0d0)
! !   print*,tem
! !   print*,tem1th
! !   print*,tem2uth
! !   dltem=1.d0-2.d0*delta; dl1=delta-1.d0
! !   t20 = (dltem)*log(m)
! !   t2der1th = (dltem)*mder1th/m
! !   t2der1dl = -2.d0*log(m)+(dltem)*mder1dl/m
! !   t2der2th = (dltem)*(mder2th/m-mder1th*mder1th/msq)
! !   t2der2dl = -4.d0*mder1dl/m+(dltem)*(mder2dl/m-mder1dl*mder1dl/msq)
! !   t2derthd = -2.d0*mder1th/m+(dltem)*(mderthd/m-mder1th*mder1dl/msq)
! ! 
! !   coef = theta*delta+1.d0
! !   den = theta*(dl1)+(coef)*m
! !   den2=den*den
! !   t30 = log(theta*(dl1)+coef*m)
! !   t3der1th = (dl1+delta*m+coef*mder1th)/den
! !   t3der1dl = (theta+theta*m+coef*mder1dl)/den
! !   t3der2th = (2.d0*delta*mder1th+coef*mder2th)/den
! !   t3der2th = t3der2th-(dl1+delta*m+coef*mder1th)**2/den2
! !   t3der2dl = (2.d0*theta*mder1dl+coef*mder2dl)/den
! !   t3der2dl = t3der2dl-(theta+theta*m+coef*mder1dl)**2/den2;
! !   t3derthd = (1.d0+m+theta*mder1th+delta*mder1dl+coef*mderthd)/den;
! !   t3derthd = t3derthd - (theta+theta*m+coef*mder1dl)*(dl1+delta*m+coef*mder1th)/den2
! ! 
! !   t1 = u1**(-theta); t2 = u2**(-theta)
! !   t1a=t1-1.d0; t2a=t2-1.d0; smlog=log(t1a)+log(t2a)
! !   tu1 = -log(u1); tu2 = -log(u2)
! !   t40 = (dl1)*(smlog)+(theta+1.d0)*(tu1+tu2)
! !   t4der1th = (dl1)*(t1*tu1/t1a+t2*tu2/t2a)+tu1+tu2  
! !   t4der1dl = smlog
! !   t4der2th = -(dl1)*(t1*tu1*tu1/t1a**2+t2*tu2*tu2/t2a**2)
! !   t4der2dl = 0.d0
! !   t4derthd = t1*tu1/t1a+t2*tu2/t2a
! ! 
! !   lpdf = t10+t20+t30+t40
! !   der1th = t1der1th+t2der1th+t3der1th+t4der1th
! !   der1dl = t1der1dl+t2der1dl+t3der1dl+t4der1dl
! !   der2th = t1der2th+t2der2th+t3der2th+t4der2th
! !   der2dl = t1der2dl+t2der2dl+t3der2dl+t4der2dl
! !   derthd = t1derthd+t2derthd+t3derthd+t4derthd
! !   lder11(1)=der1th; lder11(2)=der1dl  ! gradient of log pdf
! !   lder22(1)=der2th; lder22(2)=derthd; lder22(3)=der2dl; ! Hessian terms
! !   
! !   th2=(-1.0d0/theta)-2.0d0
! !   dl2=1.0d0-2.0d0*delta
! !   th1=-theta-1.0d0
! !   
! !   tem1=th2*mder1u/(m+1.0d0)+(dl2)*mder1u/m
! !   tem2=coef*mder1u/tem+dl1*(-theta*u1**(-theta-1.0d0)/(u1**(-theta)-1.0d0))&
! !   -(theta+1.0d0)/u1
! !   lpdf1u=tem1+tem2 !correct
! !   
! !   tem1=th2*mder1v/(m+1.0d0)+(dl2)*mder1v/m
! !   tem2=coef*mder1v/tem+dl1*(-theta*u2**(-theta-1.0d0)/(u2**(-theta)-1.0d0))&
! !   -(theta+1.0d0)/u2
! !   lpdf1v=tem1+tem2 !correct
! !   
! !   tem1=th2*(-mder1u**2.0d0/(m+1.0d0)**2.0d0+mder2u/(m+1.0d0))&
! !   +dl2*(-mder1u**2.0d0/m**2+mder2u/m)
! !   tem2=coef*(-tem1u*mder1u/tem**2.0d0+mder2u/tem)+(theta+1.0d0)/u1**2.0d0
! !   tem3=(dl1/(u1**(-theta)-1)**2.d0)*(-theta*th1*u1**(th1-1.0d0)*(u1**(-theta)-1.0d0)&
! !   -(theta*u1**th1)**2.0d0)
! !   lpdf2u=tem1+tem2+tem3 !correct
! !   
! !   
! !   tem1=th2*(-mder1v**2.0d0/mp1**2.0d0+mder2v/mp1)+dl2*(-mder1v**2.0d0/m**2+mder2v/m)
! !   tem2=coef*(-tem1v*mder1v/tem**2.0d0+mder2v/tem)+(theta+1.0d0)/u2**2.0d0
! !   tem3=(dl1/(u2**(-theta)-1)**2.d0)*(-theta*th1*u2**(th1-1.0d0)*(u2**(-theta)-1.0d0)&
! !   -(theta*u2**th1)**2.0d0)
! !   lpdf2v=tem1+tem2+tem3 !correct
! !   
! !   tem1=th2*(-mder1u*mder1v/mp1**2+mder2uv/mp1)+dl2*(-mder1u*mder1v/m**2+mder2uv/m)
! !   tem2=coef*(-tem1v*mder1u/tem**2+mder2uv/tem)
! !   lpdf2uv=tem1+tem2!correct
! !   
! !   tem1=mder1u/(mp1*thsq)+(-1.0d0/theta-2.0d0)*(-mder1th*mder1u/mp1**2+mder2uth/mp1)&
! !   +dltem*(-mder1th*mder1u/m**2+mder2uth/m)
! !   tem2=(-tem1th*tem1u/tem/tem)+tem2uth/tem
! !   tem3=(-t1/u1+t1*log(u1)*theta/u1)*t1a-t1*log(u1)*theta*t1/u1
! !   lpdf2uth=tem1+tem2+dl1*tem3/t1a**2-1.0d0/u1
! !   
! !   tem1=mder1v/(mp1*thsq)+(-1.0d0/theta-2.0d0)*(-mder1th*mder1v/mp1**2+mder2vth/mp1)&
! !   +dltem*(-mder1th*mder1v/m**2+mder2vth/m)
! !   tem2=(-tem1th*tem1v/tem/tem)+tem2vth/tem
! !   tem3=(-t2/u2+t2*log(u2)*theta/u2)*t2a-t2*log(u2)*theta*t2/u2
! !   lpdf2vth=tem1+tem2+dl1*tem3/t2a**2-1.0d0/u2
! !   
! !   tem1=th2*(-mder1dl*mder1u/mp1**2+mder2udl/mp1)-2.0d0*mder1u/m
! !   tem2=dl2*(-mder1dl*mder1u/m**2+mder2udl/m)-tem1dl*coef*mder1u/tem**2.0d0
! !   tem3=(theta*mder1u+coef*mder2udl)/tem
! !   lpdf2udl=tem1+tem2+tem3-theta*(t1/u1)/(t1-1.0d0)!correct
! !   
! !   tem1=th2*(-mder1dl*mder1v/mp1**2+mder2vdl/mp1)-2.0d0*mder1v/m
! !   tem2=dl2*(-mder1dl*mder1v/m**2+mder2vdl/m)-tem1dl*coef*mder1v/tem**2.0d0
! !   tem3=(theta*mder1v+coef*mder2vdl)/tem
! !   lpdf2vdl=tem1+tem2+tem3-theta*(t2/u2)/(t2-1.0d0)!correct
! !   
! !   lderu(1)=lpdf1u;lderu(2)=lpdf2u;lderu(3)=lpdf2uth;lderu(4)=lpdf2udl;
! !   lderv(1)=lpdf1v;lderv(2)=lpdf2v;lderv(3)=lpdf2vth;lderv(4)=lpdf2vdl;
! !   
! !   lderuv=lpdf2uv
! !   return
! ! end
! ! 
! ! 
! ! Function with Gumbel log pdf c(u1,u2,cpar) for tree 2, and 
! !   lder1, lder2 (partial wrt cpar, 1st and 2nd order)
! !   also lder1u lder2u (partial wrt u1, 1st and 2nd order)
! ! ***  also lder1v (partial wrt u2, 1st order) and others
! !   and ldermixu  (partial wrt cpar and u1) 
! ! inputs
! !   u1,u2 = values in (0,1)
! !   cpar = scalar in (1,100)
! ! outputs 
! !   lpdf = log pdf, 
! !   lder1 = \p lpdf/\p cpar
! !   lder2 = \p^2 lpdf/\p cpar^2
! !   lder1u = \p lpdf/\p u1
! !   lder2u = \p^2 lpdf/\p u1^2
! !   ldermixu = \p^2 lpdf/\p u1 \p cpar
! !   lder1v = \p lpdf/\p u2
! !   lder2v = \p^2 lpdf/\p u2^2
! !   lder2uv = \p^2 lpdf/\p u1 \p u2
! !   ldermixv = \p^2 lpdf/\p u2 \p cpar
! !  some derivs added; function name with ending of 't'
! ! subroutine  lgum2derivt(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,&
! ! ldermixu,lder1v,lder2v,lder2uv,ldermixv)
! !   implicit none
! !   double precision u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,lder2uv
! !   double precision x,y,tx,ty,xd,yd,s,m,logs,logm,msq,dlsq,dlcu
! !   double precision sder1,sder2,mder1,mder2,den,den2
! !   double precision mu,m2u,u1sq,muder1,u2sq,lder1v,lder2v,ldermixv,ldermixu
! !   double precision mv,m2v,mvder1,m2uv
! !   x = -log(u1); y = -log(u2);
! !   tx = log(x); ty = log(y); 
! !   xd = x**cpar; yd = y**cpar; s = xd+yd; m = s**(1.d0/cpar);
! !   logs=log(s);
! !   dlsq=cpar*cpar; dlcu=dlsq*cpar
! !   for 1-factor and 2-factor models
! !   sder1 = xd*tx+yd*ty;
! !   sder2 = xd*tx*tx+yd*ty*ty;
! !   mder1 = m*sder1/(s*cpar)-m*logs/dlsq;
! !   mder2 = -mder1*logs/dlsq-2.d0*m*sder1/(s*dlsq)+2.d0*m*logs/dlcu;
! !   mder2 = mder2+sder2*m/(s*cpar)+(mder1/s-m*sder1/(s*s))*sder1/cpar;
! !   den = m+cpar-1.d0; den2=den*den
! !   logm=log(m); msq=m*m
! !   lpdf = -m+log(den)+(1.d0-2.d0*cpar)*logm+(cpar-1.d0)*(tx+ty)+x+y;
! !   lder1 = -mder1+(mder1+1.d0)/den-2.d0*logm+(1.d0-2.d0*cpar)*mder1/m+tx+ty;
! !   lder2 = -mder2+mder2/den-(mder1+1.d0)**2/den2-4.d0*mder1/m&
! !   +(1.d0-2.d0*cpar)*(mder2/m-mder1*mder1/msq);
! !   for 2-factor model
! !   u1sq=u1*u1;
! !   mu = -m*xd/(u1*s*x); 
! !   m2u = (1.d0-cpar)*m*xd*xd/(u1*s*x)**2+(cpar-1.d0)*m*xd/(s*x*x*u1sq)+m*xd/(s*x*u1sq);
! !   muder1 = -(mder1/s-m*sder1/(s*s))*xd/(x*u1)-m*xd*tx/(s*u1*x);
! !   lder1u = -mu+mu/den+(1.d0-2.d0*cpar)*mu/m-(cpar-1.d0)/(u1*x)-1.d0/u1;
! !   lder2u = -m2u+m2u/den-mu*mu/den2+(1.d0-2.d0*cpar)*(m2u/m-mu*mu/msq)+(cpar-1.d0)/(u1sq*x); 
! !   lder2u = lder2u-(cpar-1.d0)/(x*x*u1sq)+1.d0/u1sq;
! !   ldermixu = -muder1+muder1/den-mu*(mder1+1.d0)/den2-2.d0*mu/m&
! !   +(1.d0-2.d0*cpar)*(muder1/m-mu*mder1/msq)-1.d0/(x*u1);
! !   
! !   u2sq=u2*u2;
! !   mv = -m*yd/(u2*s*y); 
! !   m2v = (1.d0-cpar)*m*yd*yd/(u2*s*y)**2+(cpar-1.d0)*m*yd/(s*y*y*u2sq)+m*yd/(s*y*u2sq);
! !   mvder1 = -(mder1/s-m*sder1/(s*s))*yd/(y*u2)-m*yd*ty/(s*u2*y);
! !   lder1v = -mu+mu/den+(1.d0-2.d0*cpar)*mu/m-(cpar-1.d0)/(u2*y)-1.d0/u2;
! !   lder2v = -m2u+m2u/den-mu*mu/den2+(1.d0-2.d0*cpar)*(m2u/m-mu*mu/msq)+(cpar-1.d0)/(u2sq*y); 
! !   lder2v = lder2v-(cpar-1.d0)/(y*y*u2sq)+1.d0/u2sq;
! !   ldermixv = -muder1+muder1/den-mu*(mder1+1.d0)/den2-2.d0*mu/m&
! !   +(1.d0-2.d0*cpar)*(muder1/m-mu*mder1/msq)-1.d0/(y*u2);
! !   
! !   m2uv=(x**(cpar-1.0d0)/u1/u2)*(cpar*(1.0d0/cpar-1.0d0)*m*y**(cpar-1.0d0)/s/s)
! !   lder2uv=-m2uv+(-mu*mv/den**2+m2uv/den)+(1.0d0-2*cpar)*(m2uv/m-mu*mv/m**2)
! !   return
! !   end
! !   
! ! 
! !   Function with Frank log pdf c(u1,u2,cpar) for tree 2, and 
! !   lder1, lder2 (partial wrt cpar, 1st and 2nd order)
! !   also lder1u lder2u (partial wrt u1, 1st and 2nd order)
! ! ***  also lder1v (partial wrt u2, 1st order) and others
! !   and ldermixu  (partial wrt cpar and u1) 
! ! inputs
! !   u1,u2 = values in (0,1)
! !   cpar = scalar in (-oo,oo)
! ! outputs 
! !   lpdf = log pdf, 
! !   lder1 = \p lpdf/\p cpar
! !   lder2 = \p^2 lpdf/\p cpar^2
! !   lder1u = \p lpdf/\p u1
! !   lder2u = \p^2 lpdf/\p u1^2
! !   ldermixu = \p^2 lpdf/\p u1 \p cpar
! !   lder1v = \p lpdf/\p u2
! !   lder2v = \p^2 lpdf/\p u2^2
! !   lder2uv = \p^2 lpdf/\p u1 \p u2
! !   ldermixv = \p^2 lpdf/\p u2 \p cpar
! !  some derivs added; function name with ending of 't'
! ! subroutine lfrk2derivt(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu,&
! !  lder1v,lder2v,lder2uv,ldermixv)
! !   implicit none
! !   double precision u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermixu
! !   double precision lder1v,lder2v,lder2uv,ldermixv
! !   double precision den,den1,den2,den1u,den2u,denmixu,t0,t1,t2
! !   double precision den1v,den2v,denmixv
! ! 
! !   t0 = exp(-cpar);
! !   t1 = exp(-cpar*u1);
! !   t2 = exp(-cpar*u2);
! !   den = t1+t2-t0-t1*t2;
! !   den1 = -u1*t1-u2*t2+t0+(u1+u2)*t1*t2;
! !   den2 = u1*u1*t1+u2*u2*t2-t0-(u1+u2)*(u1+u2)*t1*t2;
! !   den1u = -cpar*t1*(1.d0-t2);  
! !   den1v = -cpar*t2*(1.d0-t1);   ! added
! !   den2u = cpar*cpar*t1*(1.d0-t2);
! !   den2v = cpar*cpar*t2*(1.d0-t1); ! added
! !   denmixu = t1*(-1.d0+cpar*u1+t2-(u1+u2)*cpar*t2); 
! !   denmixv = t2*(-1.d0+cpar*u2+t1-(u1+u2)*cpar*t1); ! added 
! !    
! !   lpdf = log(abs(cpar))+log(abs(1-t0))-cpar*(u1+u2)-2*log(abs(den));
! !   121023 maybe later add the limits as cpar->0
! !   pdf = cpar*(1-t0)/den^2 where
! !      1-t0 has same sign as cpar,  den has same sign as cpar
! !   lpdf = log(cpar*(1.d0-t0))-cpar*(u1+u2)-2.d0*log(abs(den));
! !   lder1 = 1.d0/cpar+t0/(1.d0-t0)-(u1+u2)-2.d0*den1/den;
! !   lder2 = -1.d0/(cpar*cpar)-t0/((1.d0-t0)*(1.d0-t0))-2.d0*den2/den&
! !   +2.d0*den1*den1/(den*den); 
! !   lder1u = -cpar-2.d0*den1u/den;
! !   lder1v = -cpar-2.d0*den1v/den;  ! added
! !   lder2u = -2.d0*den2u/den+2.d0*den1u*den1u/(den*den);
! !   lder2v = -2.d0*den2v/den+2.d0*den1v*den1v/(den*den); ! added
! !   lder2uv = 2.d0*cpar*cpar*t1*t2/den+2.d0*den1u*den1v/(den*den); ! added
! !   ldermixu = -1.d0-2.d0*denmixu/den+2.d0*den1u*den1/(den*den);
! !   ldermixv = -1.d0-2.d0*denmixv/den+2.d0*den1v*den1/(den*den); ! added
! !   return
! !   end
! ! 
! ! 
! !   Function with bb8 logpdf =log f(u1,u2,nu,delta) 
! ! lder1, lder2 (partial wrt rho, 1st and 2nd order)
! !  also lderu (partial wrt u1, 1st and 2nd order and also wrt parmaters)
! ! ***  also lderv (partial wrt u2, 1st order and 2nd order and also wrt parameters) and others
! ! and lderuv  (partial wrt u and v) 
! ! outputs 
! ! lpdf = log pdf, 
! ! lder11 = \p lpdf/\p nu  \p lpdf/\p delta !2-dimensional vector
! ! lder22 = \p^2 lpdf/\p^2 nu^2 \p^2 \p^2 lpdf/\p^2 nu delta lpdf/\p^2 delta^2  !3-dimensional vector
! ! lderu = \p lpdf/\p u, \p^2 lpdf/\p u^2, \p^2 lpdf/\p u nu, \p^2 lpdf/\p u delta 
! ! lderv = \p lpdf/\p v, \p^2 lpdf/\p v^2, \p^2 lpdf/\p v nu, \p^2 lpdf/\p v delta 
! ! lderuv= \p^2 lpdf/\p uv
! ! subroutine lbb8derivs(u1,u2,cparv,lpdf,lder11,lder22,lderu,lderv,lderuv)
! !  implicit none
! !  double precision u1,u2,cparv(2),lpdf,lder11(2),lder22(3)
! !  double precision eta,x,y,delta,dl1,dlu1,dlu2,nu,nu1,nu2
! !  double precision eta1nu,eta1del,x1nu,x1del,y1nu,y1del,eta2nu,eta2del,eta2nudel
! !  double precision tem,temm,tem1,tem2,tem3,tem4,lpdf1nu,lpdf1del,tem1nu,tem1dl,eta1dl
! !  double precision x2nu,y2nu,y2nudel,y2del,x2nudel,x2del,lpdf2nu,lpdf2del,lpdf2nudel
! !  double precision xy1mixdel,tem2nu,tem2del,tem2nudel
! !  double precision lderu(4),lderv(4),lderuv,lpdf1u,lpdf2u,lpdf1v,lpdf2v
! !  double precision lpdf2udl,lpdf2unu,lpdf2vdl,lpdf2vnu
! !  double precision tem1u,tem1v,x1u,y1v,tem2u,tem2v,x2u,y2v,x2unu,y2vnu
! !  double precision x2udl,y2vdl,tem2unu,tem2udl,tem2vnu,tem2vdl
! !  double precision tem2uv
! !  nu=cparv(1)
! !  delta=cparv(2)
! !  dl1=1.0d0-delta
! !  nu1=nu-1.0d0
! !  nu2=1.0d0/nu-2.0d0
! !  dlu1=delta*u1
! !  dlu2=delta*u2
! !  eta=1.0d0-dl1**nu
! !  x=1.0d0-(1.0d0-delta*u1)**nu
! !  y=1.0d0-(1.0d0-delta*u2)**nu
! !  tem=x*y/eta !x*y/eta
! !  temm=1.0d0-tem
! !  lpdf=-log(eta)+log(delta)+(nu2)*log(temm)+log(nu-tem)+(nu1)*log(1.0d0-dlu1)+&
! !  (nu1)*log(1-dlu2)
! !  
! !  eta1nu=-log(dl1)*((dl1)**nu)!\p eta/\p nu
! !  eta1del=nu*((dl1)**(nu1)) !\p eta/\p delta
! !  x1nu=-log(1.0d0-dlu1)*((1-dlu1)**nu) !\p x/\p nu
! !  x1del=(nu*(1.0d0-dlu1)**(nu1))*u1 !\p x/\p delta
! !  y1nu=-log(1.0d0-dlu2)*((1-delta*u2)**nu)  !\p y/\p nu
! !  y1del=nu*((1.0d0-dlu2)**(nu1))*u2 !\p y/\p delta
! !  tem1nu=-(x*y/eta**2)*eta1nu+(x1nu*y+y1nu*x)/eta !\p tem/\p nu
! !  tem1dl=(x1del*y+y1del*x)/eta-eta1del*x*y/eta**2!\p tem/\p delta
! !  
! !  tem1=-eta1nu/eta-log(temm)/nu**2+nu2*(1.0d0/(temm))*(-tem1nu)
! !  tem2=(1.0d0/(nu-tem))*(1.0d0-tem1nu)+log(1.0d0-dlu1)+log(1.0d0-dlu2)
! !  lpdf1nu=tem1+tem2!\p lpdf/\p nu
! !  
! !  tem1=-eta1del/eta+1.0d0/delta+nu2*(1.0d0/temm)*(-tem1dl)
! !  tem2=-(1.0d0/(nu-tem))*(tem1dl)-(nu1)*u1/(1.0d0-delta*u1)-nu1*u2/(1.0d0-delta*u2)
! !  lpdf1del=tem1+tem2 !\p lpdf/\p delta
! !  
! !  lder11(1)=lpdf1nu;lder11(2)=lpdf1del;
! !  
! !  eta2nu=-((log(dl1))**2)*(dl1**nu) !\p^2 eta/\p nu^2
! !  eta2del=-nu*(nu1)*(dl1)**(nu-2.0d0) !\p^2 eta/\p delta^2
! !  eta2nudel=nu*log(dl1)*(dl1**nu1)+dl1**nu1!\p^2 eta/\p nu delta
! !  x2nu=-(log(1.0d0-dlu1)**2)*((1.0d0-dlu1)**nu) !\x^2 eta/\p nu^2
! !  y2nu=-(log(1.0d0-delta*u2)**2)*((1.0d0-dlu2)**nu)!\y^2 eta/\p nu^2
! !  x2del=-u1*u1*nu*nu1*((1-dlu1)**(nu-2.0d0))!\x^2 eta/\p delta^2
! !  y2del=-u2*u2*nu*nu1*((1-dlu2)**(nu-2.0d0))!\y^2 eta/\p delta^2
! !  x2nudel=u1*nu*(1-dlu1)**(nu1)*log(1-dlu1)+u1*((1-dlu1)**nu)/(1-dlu1) !\x^2 eta/\p nu delta
! !  y2nudel=u2*nu*(1-dlu2)**(nu1)*log(1-dlu2)+u2*((1-dlu2)**nu)/(1-dlu2)!\y^2 eta/\p nu delta
! !  
! !  tem1=(2.0d0/eta**3)*x*y*((eta1nu)**2)-(eta2nu*x*y+(x1nu*y+y1nu*x)*eta1nu)/eta**2
! !  tem2=-eta1nu*x1nu*y/eta**2+(x2nu*y+x1nu*y1nu)/eta
! !  tem3=-eta1nu*y1nu*x/eta**2+(y2nu*x+x1nu*y1nu)/eta
! !  tem2nu=tem1+tem2+tem3 !\p^2 tem/\p nu^2
! !  
! !  xy1mixdel=x1del*y1del
! !  tem1=(2.0d0*x*y*((eta1del)**2))/eta**3-(eta2del*x*y+(x1del*y+y1del*x)*eta1del)/eta**2
! !  tem2=-eta1del*y*x1del/eta**2+(xy1mixdel+x2del*y)/eta
! !  tem3=-eta1del*x*y1del/eta**2+(xy1mixdel+y2del*x)/eta
! !  tem2del=tem1+tem2+tem3!\p^2 tem/\p delta^2
! !  
! !  
! !  tem1=2.0d0*eta1nu*eta1del*x*y/eta**3.0d0
! !  tem2=-(eta2nudel*x*y+(x1del*y+y1del*x)*eta1nu)/eta**2.0d0
! !  tem3=-eta1del*x1nu*y/eta**2+(x2nudel*y+y1del*x1nu)/eta
! !  tem4=-eta1del*y1nu*x/eta**2+(y2nudel*x+x1del*y1nu)/eta
! !  tem2nudel=tem1+tem2+tem3+tem4 !\p^2 tem/\p nu delta
! !  
! !  tem1=eta1nu**2/eta**2-eta2nu/eta+2.0d0*log(temm)/nu**3.0d0
! !  tem2=tem1nu/(temm)/nu**2.0d0+tem1nu/(temm)/nu**2.0d0
! !  tem3=(1.0d0/nu-2.0d0)*(-(tem1nu**2.0d0/(temm)**2.0d0)-tem2nu/(temm))
! !  tem4=-(1.0d0-tem1nu)**2.0d0/(nu-tem)**2.0d0-tem2nu/(nu-tem)
! !  lpdf2nu=tem1+tem2+tem3+tem4!\p^2 lpdf/\p nu^2 
! !  
! !  tem1=eta1del**2/eta**2-eta2del/eta-1.0d0/delta**2
! !  tem2=(-(1.0d0/(1.0d0-tem)**2)*(tem1del)**2+(-tem2del)/(1.0d0-tem)))*(1.0d0/nu-2.0d0)
! !  tem2=nu2*(-tem1dl**2/(temm)**2-tem2del/temm)
! !  tem3=-(tem1dl**2/(nu-tem)**2+tem2del/(nu-tem))
! !  tem4=-(nu1*u1**2.0d0/(1-dlu1)**2.0d0+nu1*u2**2.0d0/(1-dlu2)**2.0d0)
! !  lpdf2del=tem1+tem2+tem3+tem4!\p^2 lpdf/\p delta^2 
! !  
! !  tem1=eta1nu*eta1del/eta**2-eta2nudel/eta+tem1dl/nu**2/(1.0d0-tem)
! !  tem2=nu2*(-(tem1nu*tem1dl/(1-tem)**2)-tem2nudel/(1.0d0-tem))
! !  tem3=(1.0d0-tem1nu)*tem1dl/(nu-tem)**2-tem2nudel/(nu-tem)
! !  tem4=-u1/(1.0d0-dlu1)-u2/(1.0d0-dlu2)
! !  lpdf2nudel=tem1+tem2+tem3+tem4
! !  
! !  lder22(1)=lpdf2nu
! !  lder22(2)=lpdf2nudel
! !  lder22(3)=lpdf2del
! !  
! !  x1u=delta*nu*(1-dlu1)**nu1
! !  y1v=delta*nu*(1-dlu2)**nu1
! !  tem1u=y*x1u/eta
! !  tem1v=x*y1v/eta
! !  
! !  lpdf1u=-nu2*tem1u/temm-tem1u/(nu-tem)-nu1*delta/(1-dlu1)
! !  lpdf1v=-nu2*tem1v/temm-tem1v/(nu-tem)-nu1*delta/(1-dlu2)
! !  
! !  x2u=-(delta**2)*nu*nu1*(1-dlu1)**(nu-2.0d0)
! !  y2v=-(delta**2)*nu*nu1*(1-dlu2)**(nu-2.0d0)
! !  tem2u=y*x2u/eta
! !  tem2v=x*y2v/eta
! !  
! !  lpdf2u=-nu2*tem1u**2/temm**2-nu2*(tem2u/temm)-tem1u**2/(nu-tem)**2&
! !  -tem2u/(nu-tem)-nu1*delta**2/(1-dlu1)**2
! !  lpdf2v=-nu2*tem1v**2/temm**2-nu2*(tem2v/temm)-tem1v**2/(nu-tem)**2&
! !  -tem2v/(nu-tem)-nu1*delta**2/(1-dlu2)**2
! !  
! !  x2unu=delta*(1-dlu1)**nu1+delta*nu*log(1-dlu1)*(1-dlu1)**nu1
! !  tem2unu=-eta1nu*y*x1u/eta**2+(y1nu*x1u+x2unu*y)/eta 
! !  tem1=tem1u/nu**2/temm+(-tem1nu*tem1u/temm**2-tem2unu/temm)*nu2
! !  tem2=tem1u*(1-tem1nu)/(nu-tem)**2-tem2unu/(nu-tem)-delta/(1-dlu1)
! !  lpdf2unu=tem1+tem2
! !  
! !  eta1dl=nu*dl1**nu1
! !  x2udl=nu*(1.0d0-dlu1)**nu1-delta*nu*nu1*u1*(1.0d0-dlu1)**(nu1-1.0d0)
! !  tem2udl=-eta1dl*y*x1u/eta**2+(y1del*x1u+x2udl*y)/eta
! !  tem1=-nu2*tem1dl*tem1u/temm**2-nu2*tem2udl/temm
! !  tem2=-tem1dl*tem1u/(nu-tem)**2-tem2udl/(nu-tem)-nu1/(1-dlu1)**2
! !  lpdf2udl=tem1+tem2
! !  
! !  y2vnu=delta*(1-dlu2)**nu1+delta*nu*log(1-dlu2)*(1-dlu2)**nu1
! !  tem2vnu=-eta1nu*x*y1v/eta**2+(x1nu*y1v+y2vnu*x)/eta 
! !  tem1=tem1v/nu**2/temm+(-tem1nu*tem1v/temm**2-tem2vnu/temm)*nu2
! !  tem2=tem1v*(1-tem1nu)/(nu-tem)**2-tem2vnu/(nu-tem)-delta/(1-dlu2)
! !  lpdf2vnu=tem1+tem2
! !  
! !  y2vdl=nu*(1.0d0-dlu2)**nu1-delta*nu*nu1*u2*(1.0d0-dlu2)**(nu1-1.0d0)
! !  tem2vdl=-eta1dl*x*y1v/eta**2+(x1del*y1v+y2vdl*x)/eta
! !  tem1=-nu2*tem1dl*tem1v/temm**2-nu2*tem2vdl/temm
! !  tem2=-tem1dl*tem1v/(nu-tem)**2-tem2vdl/(nu-tem)-nu1/(1-dlu2)**2
! !  lpdf2vdl=tem1+tem2
! !  
! !  tem2uv=y1v*x1u/eta
! !  tem1=-nu2*tem1u*tem1v/temm**2-tem2uv*nu2/temm
! !  tem2=-tem1u*tem1v/(nu-tem)**2-tem2uv/(nu-tem)
! !  lderuv=tem1+tem2
! !  
! !  
! !  lderu(1)=lpdf1u;lderu(2)=lpdf2u;lderu(3)=lpdf2unu;lderu(4)=lpdf2udl;
! !  lderv(1)=lpdf1v;lderv(2)=lpdf2v;lderv(3)=lpdf2vnu;lderv(4)=lpdf2vdl;
! ! return
! ! end
! !  
! ! 
! !  compute the density of frank copula  
! ! subroutine dfrk(x1,x2,cpar,pdf)
! !   implicit none
! !   double precision:: x1,x2
! !   double precision :: cpar,tem1,tem2,tem,t1,tem3,pdf
! ! 
! ! 
! !   t1=1.d0-exp(-cpar);
! !   tem1=exp(-cpar*x1); tem2=exp(-cpar*x2);
! !   tem=t1-(1.d0-tem1)*(1.d0-tem2);
! !   tem3=cpar*tem1*tem2*t1
! !   pdf=tem3/(tem*tem);
! ! return
! ! end subroutine dfrk
! ! 
! ! 
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
! ! subroutine mderivs(u1,u2,cparv,m,mder1th,mder1dl,mder2th,mder2dl,mderthd,&
! ! mderu,mderv,mder2uv)
! !   implicit none
! !   double precision u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd
! !   double precision t1,t2,tu1,tu2,ttu1,ttu2,td01,td02,td11,td12,td21,td22
! !   double precision s,sder1th,sder1dl,sder2th,sder2dl,sderthd,m1,ts,m1der1dl 
! !   double precision t1a,t2a,dlsq,dlcu
! !   double precision xder1u,yder1v,xder2u,mder1u,mder1v,mder2u,mder2uv,mder2v,del1,th1
! !   double precision yder2v,xder2uth,yder2vth,xder1th,xder2udl,yder2vdl,yder1th
! !   double precision mder2vdl,mder2vth,mder2udl,mder2uth
! !   double precision mderu(4),mderv(4),cparv(2)
! !   theta=cparv(1);delta=cparv(2)
! !   del1=1.0d0/delta-1.0d0
! !   th1=-theta-1.0d0
! !   
! !   t1 = u1**(-theta); t2 = u2**(-theta)
! !   t1a=t1-1.d0; t2a=t2-1.d0
! !   tu1 = -log(u1); tu2 = -log(u2)
! !   ttu1 = log(t1a); ttu2 = log(t2a)
! !   td01 = (t1a)**delta; td02 = (t2a)**delta
! !   td11=td01/t1a; td12=td02/t2a;
! !   td21=td11/t1a; td22=td12/t2a;
! !  
! !   s = td01+td02
! !   sder1th = delta*(td11*t1*tu1+td12*t2*tu2)
! !   sder1dl = td01*ttu1+td02*ttu2
! !   sder2th = delta*(delta-1.d0)*(td21*t1*t1*tu1*tu1+td22*t2*t2*tu2*tu2)
! !   sder2th = sder2th+delta*(td11*t1*tu1*tu1+td12*t2*tu2*tu2)
! !   sder2dl = td01*ttu1*ttu1+td02*ttu2*ttu2
! !   sderthd = sder1th/delta+delta*(td11*ttu1*tu1*t1+td12*ttu2*tu2*t2)
! ! 
! !   m = s**(1.d0/delta); m1 = m/s
! !   ts = log(s)
! !   dlsq=delta*delta; dlcu=delta*dlsq
! !   mder1th = m1*sder1th/delta
! !   mder1dl = m1*sder1dl/delta - m*ts/dlsq
! !   m1der1dl = mder1dl/s - m*sder1dl/s**2
! !   mder2th = (1.d0-delta)*m1*sder1th**2/(dlsq*s)+m1*sder2th/delta
! !   mder2dl = 2.d0*m*ts/dlcu-mder1dl*ts/dlsq-2.d0*m1*sder1dl/dlsq
! !   mder2dl = mder2dl+sder2dl*m1/delta+sder1dl*m1der1dl/delta
! !   mderthd = -m1*sder1th/dlsq+sder1th*m1der1dl/delta+m1*sderthd/delta
! !   
! !   xder1u=-delta*theta*td11*u1**(th1) !correct
! !   yder1v=-delta*theta*td12*u2**(th1) !correct
! !   xder2u=-delta*theta*(-theta*(delta-1.0d0)*td21*(u1**(2*th1))+&
! !   th1*td11*u1**(th1-1.0d0))!correct
! !   yder2v=-delta*theta*(-theta*(delta-1.0d0)*td22*(u2**(2*th1))+&
! !   th1*td12*u2**(th1-1.0d0))!correct
! !   xder2uth=-delta*td11*u1**th1&
! !   -delta*theta*(-(delta-1.0d0)*td21*u1**(-theta)*log(u1)*u1**th1&
! !   -u1**th1*log(u1)*td11)!correct
! !   yder2vth=-delta*td12*u2**th1-&
! !   delta*theta*(-(delta-1.0d0)*td22*u2**(-theta)*log(u2)*u2**th1&
! !   -u2**th1*log(u2)*td12)!correct
! !   xder1th=delta*(td11*t1*tu1)
! !   yder1th=delta*(td12*t2*tu2)
! !   xder2udl=-theta*u1**th1*(td11+delta*td11*ttu1)
! !   yder2vdl=-theta*u2**th1*(td12+delta*td12*ttu2)
! !   
! !   mder1u=xder1u*s**del1/delta !correct
! !   mder1v=yder1v*s**del1/delta !correct
! !   mder2uv=del1*(s**(del1-1))*yder1v*xder1u/delta!correct
! !   mder2u=del1*s**(del1-1.0d0)*xder1u**2/delta+xder2u*s**del1/delta!correct
! !   mder2v=del1*(s**(del1-1.0d0))*(yder1v)**2/delta+yder2v*s**(del1)/delta!correct
! !   mder2uth=(xder2uth*s**del1+del1*s**(del1-1.0d0)*sder1th*xder1u)/delta !correct
! !   mder2udl=xder2udl*s**del1/delta+(-s**del1/delta**2+&
! !   (1.0d0/delta)*s**(del1)*(-log(s)/delta**2+sder1dl*del1/s))*xder1u
! !   mder2vth=(yder2vth*s**del1+del1*s**(del1-1.0d0)*sder1th*yder1v)/delta!correct
! !   mder2vdl=yder2vdl*s**del1/delta+(-s**del1/delta**2+&
! !   (1.0d0/delta)*s**(del1)*(-log(s)/delta**2+sder1dl*del1/s))*yder1v
! !   
! !   mderu(1)=mder1u;mderu(2)=mder2u;mderu(3)=mder2uth;mderu(4)=mder2udl;
! !   mderv(1)=mder1v;mderv(2)=mder2v;mderv(3)=mder2vth;mderv(4)=mder2vdl;
! ! 
! !   return
! !   
! !   end
! ! 
! ! 
! ! 
! !   this routine is in bb1facts.f90
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
! subroutine mderivs2(u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd)
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