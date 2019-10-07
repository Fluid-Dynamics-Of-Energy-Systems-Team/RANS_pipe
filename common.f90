      


!     ****************      
      integer         istep,nsteps,inlet,centerBC,numDomain
      common /intvar/ istep,nsteps,inlet,centerBC,numDomain
      save   /intvar/




!     ****************
      real*8          pi,dpdz
      common /numvar/ pi,dpdz
      save   /numvar/



!     ****************
      real*8          fmu(0:i1,0:k1),f1(0:i1,0:k1),f2(0:i1,0:k1),ReT(0:i1,0:k1),yp(0:i1,0:k1),ReTauS(0:i1,0:k1),&
                      ypt(0:i1,0:k1),dterm(0:i1,0:k1),eterm(0:i1,0:k1),fv2(imax,kmax),Lh(imax,kmax), & 
                                 at(imax,kmax),Lht(imax,kmax),sigmat,sigmak,sigmae,cmu,ce1,ce2,sigmah2,atmp(0:i1,0:k1), &
                     bF1(0:i1,0:k1),bF2(imax,kmax),cdKOM(imax,kmax) 
      common /turb/   fmu,f1,f2,ReT,yp,ypt,dterm,eterm,fv2,Lh,at,Lht,sigmat,sigmak,sigmae, &
                     cmu,ce1,ce2,sigmah2,atmp,ReTauS, bF1, bF2, cdKOM
      save   /turb/




!     ****************
      real*8          Ru(0:i1),Rp(0:i1),y_fa(0:i1),y_cv(0:i1),dru(0:i1),drp(0:i1),dz,z1(0:k1),z2(0:k1),wallDist(1:imax)
      common /phsgrd/ Ru,Rp,y_fa,y_cv,dru,drp,dz,z1,z2,wallDist
      save   /phsgrd/



!     ****************
      real*8          ekm(0:i1,0:k1),ekmt(0:i1,0:k1),ekme(0:i1,0:k1),ekh(0:i1,0:k1),cp(0:i1,0:k1),&
                     temp(0:i1,0:k1),tco(0:i1,0:k1),Pk(0:i1,0:k1),Gk(0:i1,0:k1),peclet(0:i1,0:k1),beta(0:i1,0:k1),&
                     ekhi(0:i1,0:k1),ekhk(0:i1,0:k1),ekmi(0:i1,0:k1),ekmk(0:i1,0:k1),cpi(0:i1,0:k1),cpk(0:i1,0:k1), &
                     Tt(0:i1,0:k1),enth_wall
      common /energy/ ekm,ekmt,ekme,ekh,cp,temp,tco,Pk,Gk,peclet,beta,ekhi,ekhk,ekmi,ekmk,cpi,cpk,Tt,enth_wall
      save   /energy/
      
      

!     ****************
      real*8          dt,dtmax
      common /timvar/ dt,dtmax
      save   /timvar/





!     ****************
      real*8          Unew(0:i1,0:k1),Vnew(0:i1,0:k1),Wnew(0:i1,0:k1),v2new(0:i1,0:k1),h2new(0:i1,0:k1), &
                      Cnew(0:i1,0:k1),knew(0:i1,0:k1),enew(0:i1,0:k1),qcrit(0:i1,0:k1),nuSAnew(0:i1,0:k1),omNew(0:i1,0:k1)
      common /new/    Unew,Vnew,Wnew,v2new,h2new,Cnew,knew,enew,qcrit,nuSAnew,omNew
      save   /new/





!     ****************
      real*8          Uold(0:i1,0:k1),Vold(0:i1,0:k1),eold(0:i1,0:k1),v2old(0:i1,0:k1), &
                      Wold(0:i1,0:k1),Cold(0:i1,0:k1),kold(0:i1,0:k1),h2old(0:i1,0:k1),nuSAold(0:i1,0:k1),omOld(0:i1,0:k1)
      common /old/    Uold,Vold,eold,v2old,Wold,Cold,kold,h2old,nuSAold,omOld
      save   /old/




!     ****************  
      real*8          dUdt(0:i1,0:k1),dVdt(0:i1,0:k1),dWdt(0:i1,0:k1)
      common /derive/ dUdt,dVdt,dWdt
      save   /derive/
      


!     ****************
      real*8          rold(0:i1,0:k1),rnew(0:i1,0:k1)
      common /density/rold,rnew
      save   /density/




!      ****************
      real*8          work(iwork),save(isave),bdrs(kmax), &
                      bdts(imax,kmax),bdzs(imax), &
                      p(imax,kmax),phi(imax,kmax)
      common /press/  work,save,bdrs,bdts,bdzs,p,phi
      save   /press/




!     ****************
      real*8          tempTab(1:nTab),rhoTab(1:nTab),betaTab(1:nTab), &
                      muTab(1:nTab),lamTab(1:nTab), &
                      cpTab(1:nTab),enthTab(1:nTab),lamocpTab(1:nTab), &
                      temp2Tab(1:nTab),rho2Tab(1:nTab),beta2Tab(1:nTab), &
                      mu2Tab(1:nTab),lam2Tab(1:nTab), &
                      cp2Tab(1:nTab),enth2Tab(1:nTab),lamocp2Tab(1:nTab)
      common /table/  tempTab,rhoTab,betaTab,muTab,lamTab,cpTab,enthTab,lamocpTab, &
                     temp2Tab,rho2Tab,beta2Tab,mu2Tab,lam2Tab,cp2Tab,enth2Tab,lamocp2Tab
      save   /table/


           
      