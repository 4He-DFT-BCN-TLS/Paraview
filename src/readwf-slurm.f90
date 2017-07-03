!Use Para_DerivnD
Implicit Real*8(A-H,O-Z)
Complex   (Kind=8), Allocatable :: psi0(:,:,:), psi02D(:,:), psi2D(:,:), Dxpsi2D(:,:), Dypsi2D(:,:)
Complex   (Kind=8), Allocatable :: Sto1c(:,:,:), Sto2c(:,:,:),Sto3c(:,:,:)
Complex   (Kind=8), Allocatable :: DDxpsi2D(:,:), DDypsi2D(:,:), DDxypsi2D(:,:)
Real      (Kind=8), Allocatable :: den0(:,:,:), x0(:), y0(:), z0(:), den2D(:,:),x(:), y(:), z(:), Wx(:), Wy(:), &  
                                   x0v(:), y0v(:), xv(:), yv(:)
                               
Integer   (Kind=4) :: nn(2), Icon=13, nnn(3), i_minloc(2)
Complex   (Kind=8) :: ci=(0.d0,1.d0), caux, caux1
Character (len =1) :: cchar="#", Vortex_axis='Z'
Logical            :: limp, L2Dplot=.false.,Ldensity=.true., Ldynamic=.false., L3Dplot=.false.,                 &
                      Lvortex=.false., Lvortex_position=.false.,Lxyv(2),Lreal_wave_function=.false.,            &
                      Lvortex_position_as_a_function_of_z=.false.,L_3D_filter=.false.
Character (len=80) :: File5='readwf.dat', File6='readwf.res', File31='position.dat', File32='velocity.dat'
Character (len=80) :: File7='denxy.dat', File8='current.out',File10='den.inp',File9='vortex_position.out',      &
                      File11='denxz.dat',File12='denyz.dat',File13='angular_momentum.dat', File15='Q2.out'
Character (len=80) :: File20='denxy.vtk',File21='denxz.vtk',File22='denyz.vtk',File23='current.vtk',File41='denxyz.vtk'
Data nx/436/, ny/436/, nz/436/, hx/-1.d0/, hy/-1.d0/, hz/-1.d0/, npd/13/,Km1/4/, ndmax/2/, nthread/4/, npi/4/
Data fac/158.66624d0/,epsrho/1.d-6/,drop_radius/15.6d0/  !drop_radius = 2.2*N_He**(1/3)/Sqrt(2)
Data xi/-40.d0/,xf/40.d0/,yi/-40.d0/,yf/40.d0/,xlv/-20.d0/,xrv/+20.d0/, ylv/-20.d0/, yrv/+20.d0/, nv/2/
Data zv_min/-20.d0/, zv_max/20.d0/
Data X_P_filter/100.d0/,Y_P_filter/100.d0/,Z_P_filter/100.d0/, X_N_filter/-100.d0/,Y_N_filter/-100.d0/,Z_N_filter/-100.d0/
!
! xlv, xrv, ylv, ylv:   valors per acotar la búsqueda de la posicoó dels vortex 
! nv numero de vortex que busquerem
!
Namelist/Input/File4,File6,File7,File8,File10,File11,File12,File13,L2Dplot,Ldynamic,xlv,xrv,ylv,yrv,nv,Lvortex_position,   &
File20,File21,File22,File23,nx,ny,nz,hx,hy,hz,npd,npi,Km1,icon,epsrho,xi,xf,yi,yf,nthread,Vortex_axis,Lvortex,File9,       &
Lreal_wave_function,drop_radius, File15, File31, File32, L3Dplot, File41,Lvortex_position_as_a_function_of_z,zv_min,zv_max,&
L_3D_filter,X_P_filter,Y_P_filter,Z_P_filter,X_N_filter,Y_N_filter,Z_N_filter
read(5,nml=Input)
K=npd    ! Km1 will be the number of derivatives for the Taylor expansion
!Open(Unit=6,File=File6)
!Write(6,Input)
Open(10,File=File10, Status='Old', buffered='yes')
Open(31,File=File31, buffered='yes')
!Open(32,File=File32)
!Allocate(x(nx)); Allocate(y(ny)); Allocate(z(nz))
!Allocate(Dxpsi2D(nx,ny)); Allocate(Dypsi2D(nx,ny))
Go to 20
10 Continue
  Ldensity=.false.
20 Continue  
Rewind(10)
call titols(10,cchar,isalto)
If(Ldynamic)Then
  Read(10,*)xmax0,ymax0,zmax0,hx0,hy0,hz0,nx0,ny0,nz0
  Read(10,*)ximp,yimp,zimp
  Read(10,*)vximp,vyimp,vzimp
Else        
  Read(10,*)xmax0,ymax0,zmax0,hx0,hy0,hz0,nx0,ny0,nz0,limp,ximp,yimp,zimp
  vximp=0.d0; vyimo=0.d0; vzimp=0.d0
EndiF  
!Write(6,'("xmax0,ymax0,zmax0,hx0,hy0,hz0,nx0,ny0,nz0...:",/,1p6E15.6,0p,3I5)') &
!           xmax0,ymax0,zmax0,hx0,hy0,hz0,nx0,ny0,nz0
If(.Not.Allocated(x0)    )Allocate(x0(nx0)          )
If(.Not.Allocated(y0)    )Allocate(y0(ny0)          )
If(.Not.Allocated(z0)    )Allocate(z0(nz0)          )
If(.Not.Allocated(psi0)  )Allocate(psi0(nx0,ny0,nz0))
!If(.Not.Allocated(Sto1c)  )Allocate(Sto1c(nx0,ny0,nz0))
!If(.Not.Allocated(Sto2c)  )Allocate(Sto2c(nx0,ny0,nz0))
!If(.Not.Allocated(Sto3c)  )Allocate(Sto3c(nx0,ny0,nz0))
If(.Not.Allocated(den0)  )Allocate(den0(nx0,ny0,nz0))
!If(.Not.Allocated(psi02D))Allocate(psi02D(nx0,ny0)  )
!If(.Not.Allocated(psi2D) )Allocate(psi2D(nx,ny)     )
!If(.Not.Allocated(den2D) )Allocate(den2D(nx,ny)     )
If(Ldensity)Then
  If(Lreal_wave_function)Then
    Read(10,*,Err=10)den0
!    Write(6,'("We have read a real w.f.")')
    psi0=Cmplx(den0)
    den0 = Abs(den0)**2
  Else  
    Read(10,*,Err=10)den0
!    Write(6,'("We have read a density")')
    psi0=Cmplx(Sqrt(den0))
  Endif  
Else  
  Read(10,*)psi0
!  Write(6,'("We have read a complex w.f.")')
  den0=Abs(psi0)**2
EndIf
Write(31,'(1p,3E18.10)')ximp,yimp,zimp
Close(31)
!Write(32,'(1p,3E18.10)')vximp,vyimp,vzimp
Close(10)

 !   Write(6,'(" We write the value of MaxLoc(den0(:,:,:)(x,y,z)....:")')
 !   Write(6,*)MaxLoc(den0(:,:,:))
 !   Write(6,'(" We write the value of MaxVal(den0(:,:,:)(x,y,z)....:")')
 !   Write(6,*)MaxVal(den0(:,:,:))

!If(hx.Lt.0.0d0)Then
!  hx=xmax0/nx*2.
!!  Write(6,'("New hx...:",1p,E15.6)')hx
!EndIf  
!If(hy.Lt.0.0d0)Then
!  hy=ymax0/ny*2.
!!  Write(6,'("New hy...:",1p,E15.6)')hy
!EndIf  
!If(hz.Lt.0.0d0)Then
!  hz=zmax0/nz*2.
!!  Write(6,'("New hz...:",1p,E15.6)')hz
!EndIf  

dxyz0=hx0*hy0*hz0
dxyz=dxyz0

!Write(6,'("den(0,0,0)....:",1p,E15.6)')den0(nx0/2+1,ny0/2+1,nz0/2+1)
n4 = (sum(den0)*dxyz0 + 0.5)
!Write(6,'("Nombre de atoms...:",I10)')n4
!Write(6,'("Norma(den)........:",1p,E15.6)')sum(den0)*dxyz0
!Write(6,'("Norma(psi)........:",1p,E15.6)')sum((Abs(psi0)**2))*dxyz0
Do ix=1,nx0; x0(ix)=-xmax0+(ix-1)*hx0; EndDo
Do iy=1,ny0; y0(iy)=-ymax0+(iy-1)*hy0; EndDo
Do iz=1,nz0; z0(iz)=-zmax0+(iz-1)*hz0; EndDo

!Open(Unit=1,File="den-i-x.dat")
!Do ix=1,nx0; Write(1,'(1p,2E15.6)')x0(ix),den0(ix,ny0/2+1,nz0/2+1); EndDo
!Close (1)
!
!Open(Unit=1,File="den-i-y.dat")
!Do iy=1,ny0; Write(1,'(1p,2E15.6)')y0(iy),den0(nx0/2+1,iy,nz0/2+1); EndDo
!Close (1)
!
!Open(Unit=1,File="den-i-z.dat")
!Do iz=1,nz0; Write(1,'(1p,2E15.6)')z0(iz),den0(nx0/2+1,ny0/2+1,iz); EndDo
!Close (1)
!
!xmax=(nx/2)*hx; ymax=(ny/2)*hy; zmax=(nz/2)*hz
!Do ix=1,nx; x(ix)=-xmax+(ix-1)*hx; EndDo
!Do iy=1,ny; y(iy)=-ymax+(iy-1)*hy; EndDo
!Do iz=1,nz; z(iz)=-zmax+(iz-1)*hz; EndDo

!
!       We store the density corresponding to the z=0 plane
!

!Psi02D = Psi0(:,:,nz0/2+1)
!
!Call INTERxy(nx,ny,x,y,Psi2D,nx0,ny0,hx0,hy0,Psi02D,K,KM1)
!
!den2D = Max(Abs(psi2D)**2,1.d-99)
!    Write(6,'(" We write the value of MaxLoc(den2Di(:,:)(x,y,0)....:")')
!    Write(6,*)MaxLoc(den2D(:,:))
!    Write(6,'(" We write the value of MaxVal(den2Di(:,:)(x,y,0)....:")')
!    Write(6,*)MaxVal(den2D(:,:))
!
!Open(Unit=1,File="den-o-x.dat")
!Do ix=1,nx; Write(1,'(1p,2E15.6)')x(ix),den2D(ix,ny/2+1); EndDo
!Close (1)
!
!Open(Unit=1,File="den-o-y.dat")
!Do iy=1,ny; Write(1,'(1p,2E15.6)')y(iy),den2D(nx/2+1,iy); EndDo
!Close (1)
!
!Open(Unit=7,File=File7)
!Do ix=1, nx
!  Do iy=1,ny
!    Write(7,'(1p,3E18.10)')x(ix),y(iy),den2D(ix,iy)
!  EndDo
!  Write(7,*)
!EndDo
!Close(7)

!Open (Unit=15, File=File15)
!x2 = 0.d0; y2 = 0.d0
!Do iz = 1 , nz0
!  Do iy = 1 , ny0
!    Do ix = 1, nx0
!      x2 = x2 + den0(ix,iy,iz)*x0(ix)**2
!      y2 = y2 + den0(ix,iy,iz)*y0(iy)**2
!    EndDo
!  EndDo  
!EndDo  
!x2 = x2*dxyz/n4
!y2 = y2*dxyz/n4
!Write(15,'(1p,3E15.6)')x2,y2,(x2-y2)
!
!If(.Not.Ldensity.And.Lvortex)Then
!  Call Init_deriv_p(npd,ndmax,nthread)
!  nn(1)=nx;nn(2)=ny
!  nnn(1)=nx0;nnn(2)=ny0;nnn(3)=nz0
!!
!!  Calcularem la posició del centre de masses
!!
!          xcm=0.0d0
!          ycm=0.0d0
!          zcm=0.0d0
!          Do ix=1, nx0
!            xcm = xcm + Sum(den0(ix,:,:)*x0(ix))
!          EndDo 
!          Do iy=1, ny0
!            ycm = ycm + Sum(den0(:,iy,:)*y0(iy))
!          EndDo 
!          Do iz=1, nz0
!            zcm = zcm + Sum(den0(:,:,iz)*z0(iz))
!          EndDo 
!          xcm = xcm*dxyz0/n4
!          ycm = ycm*dxyz0/n4
!          zcm = zcm*dxyz0/n4
!
!          Write(6,'("Posició del C.M....: x_cm, y_cm, z_cm...",1p,3E15.6)')xcm,ycm,zcm
!!
!!  Calcularem <Lx>, <Ly> i <Lz>
!!
!          Call DerivnD(1,nnn,hx0,1,psi0,Sto1c,Icon)    ! Derivada respecte de x
!          Call DerivnD(1,nnn,hy0,2,psi0,Sto2c,Icon)    ! Derivada respecte de y
!          caux = 0.d0
!          Do iz=1, nz0
!            Do iy=1, ny0
!              Do ix=1, nx0
!                caux1 = ci*((y0(iy)-ycm)*sto1c(ix,iy,iz) - (x0(ix)-xcm)*sto2c(ix,iy,iz)) 
!                caux = caux + Conjg(Psi0(ix,iy,iz))*caux1
!                sto3c(ix,iy,iz) = caux1
!              EndDo
!            EndDo
!          EndDo
!          caux = caux*dxyz0
!          Write(6,'(T5,"<L_z>..............:",8x,1p,2E15.6," \hbar")')caux
!
!          Call DerivnD(1,nnn,hx0,1,Sto3c,Sto1c,Icon)    ! Derivada respecte de x
!          Call DerivnD(1,nnn,hy0,2,Sto3c,Sto2c,Icon)    ! Derivada respecte de y
!          caux = 0.d0
!          Do iz=1, nz0
!            Do iy=1, ny0
!              Do ix=1, nx0
!                caux1 = ci*((y0(iy)-ycm)*sto1c(ix,iy,iz) - (x0(ix)-xcm)*sto2c(ix,iy,iz)) 
!                caux = caux + Conjg(Psi0(ix,iy,iz))*caux1
!              EndDo
!            EndDo
!          EndDo
!          caux = caux*dxyz0
!          Write(6,'(T5,"<L_z2>.............:",8x,1p,2E15.6," \hbar")')caux
!
!          Call derivnD(1,nnn,hx0,1,psi0,sto1c,Icon)    ! Derivada respecte de x
!          Call derivnD(1,nnn,hz0,3,psi0,sto2c,Icon)    ! Derivada respecte de z
!          caux = 0.d0
!          Do iz=1, nz0
!            Do iy=1, ny0
!              Do ix=1, nx0
!                caux = caux + Ci*Conjg(Psi0(ix,iy,iz))*                  &
!                ((x0(ix)-xcm)*sto2c(ix,iy,iz) - (z0(iz)-zcm)*sto1c(ix,iy,iz)) 
!              EndDo
!            EndDo
!          EndDo
!          caux = caux*dxyz0
!          Write(6,'(T5,"<L_y>..............:",8x,1p,E15.6," \hbar")')Real(caux)
!
!          Call derivnD(1,nnn,hy0,2,psi0,sto1c,Icon)    ! Derivada respecte de y
!          Call derivnD(1,nnn,hz0,3,psi0,sto2c,Icon)    ! Derivada respecte de z
!          caux = 0.d0
!          Do iz=1, nz0
!            Do iy=1, ny0
!              Do ix=1, nx0
!                caux = caux + Ci*Conjg(Psi0(ix,iy,iz))*                  &
!                ((z0(iz)-zcm)*sto1c(ix,iy,iz) - (y0(iy)-ycm)*sto2c(ix,iy,iz)) 
!              EndDo
!            EndDo
!          EndDo
!          caux = caux*dxyz0
!          Write(6,'(T5,"<L_x>..............:",8x,1p,E15.6," \hbar")')Real(caux)
!
!  Call DerivnD(1,nn,hx,1,psi2D,Dxpsi2D,Icon)
!  Call DerivnD(1,nn,hy,2,psi2D,Dypsi2D,Icon)
!
!
!!
!! Busquem la posició dels vortex
!!
!  If(Lvortex_position)Then
!    Open(Unit=9,File=File9)      
!    If(.Not.Lvortex_position_as_a_function_of_z)Then
!
!!    ixmin = (-drop_radius + xmax)/hx + 1.5 
!!    ixmax = ( drop_radius + xmax)/hx + 1.5 
!!    iymin = (-drop_radius + ymax)/hy + 1.5 
!!    iymax = ( drop_radius + ymax)/hy + 1.5 
!
!      ixmin = (xlv + xmax)/hx + 1.5 
!      ixmax = (xrv + xmax)/hx + 1.5 
!      iymin = (ylv + ymax)/hy + 1.5 
!      iymax = (yrv + ymax)/hy + 1.5 
!
!
!
!      Write(6,'("ixmin, ixmax, iymin, iymax...:",4I10)')ixmin, ixmax, iymin, iymax
!      Write(6,'(" We write the value of MinLoc(den2Di(ixmin:ixmax,iymin:iymax)....:")')
!      Write(6,*)MinLoc(den2D(ixmin:ixmax,iymin:iymax))
!      I_MinLoc = MinLoc(den2D(ixmin:ixmax,iymin:iymax))
!      Write(6,'(" I_MinLoc(1,2)....:",2I10)')I_MinLoC
!      Write(6,'(" We write the value of MinVal(den2Di(ixmin:ixmax,iymin:iymax)....:")')
!      Write(6,*)MinVal(den2D(ixmin:ixmax,iymin:iymax))
!
!      ix1 = (I_MinLoc(1) - 1) + ixmin
!      iy1 = (I_MinLoc(2) - 1) + iymin
!
!      Write(6,'("Index del primer vortex:ix1,iy1............:",2I10)')ix1,iy1
!      Write(6,'("Posició del primer vortex:x(ix1),y(iy1)....:",1p,2E15.6)')x(ix1),y(iy1)
!
!      Write(6,'("den2D(ix1, iy1)...:",1p,E15.6)')den2D(ix1,iy1)
!      If(nv.Gt.1)Then
!        den_min = den2D(ix1,iy1)
!
!        Write(6,'(" We write the value of MinLoc(den2Di(ixmin:ixmax,iymin:iymax),Mask=(den2D(ixmin:ixmax,iymin:iymax)>den_min))....:")')
!        Write(6,*)MinLoc(den2D(ixmin:ixmax,iymin:iymax),Mask=(den2D(ixmin:ixmax,iymin:iymax)>den_min))
!        I_MinLoc = MinLoc(den2D(ixmin:ixmax,iymin:iymax),Mask=(den2D(ixmin:ixmax,iymin:iymax)>den_min))
!        Write(6,'(" I_MinLoc(1,2)....:",2I10)')I_MinLoC
!        Write(6,'(" We write the value of MinVal(den2Di(ixmin:ixmax,iymin:iymax),Mask=(den2D(ixmin:ixmax,iymin:iymax)>den_min))....:")')
!        Write(6,*)MinVal(den2D(ixmin:ixmax,iymin:iymax),Mask=(den2D(ixmin:ixmax,iymin:iymax)>den_min))
!
!        ix2 = (I_MinLoc(1) - 1) + ixmin
!        iy2 = (I_MinLoc(2) - 1) + iymin
!
!        Write(6,'("Index del segon vortex:ix2,iy2............:",2I10)')ix2,iy2
!        Write(6,'("Posició del segon vortex:x(ix2),y(iy2)....:",1p,2E15.6)')x(ix2),y(iy2)
!
!        Write(6,'("den2D(ix2, iy2)...:",1p,E15.6)')den2D(ix2,iy2)
!      Endif
!!
!! Refinem la posició dels vortex
!!
!
!      Allocate(xv(2))      
!      Allocate(yv(2))
!
!      det   = Dreal(Dxpsi2D(ix1,iy1))*Dimag(Dypsi2D(ix1,iy1)) -    &
!            Dreal(Dypsi2D(ix1,iy1))*Dimag(Dxpsi2D(ix1,iy1)) 
!      Stepx = (Dimag(psi2D(ix1,iy1))*Dreal(Dypsi2D(ix1,iy1)) -     &
!               Dreal(psi2D(ix1,iy1))*Dimag(Dypsi2D(ix1,iy1)))/det 
!      Stepy = (Dreal(psi2D(ix1,iy1))*Dimag(Dxpsi2D(ix1,iy1)) -     &
!               Dimag(psi2D(ix1,iy1))*Dreal(Dxpsi2D(ix1,iy1)))/det 
!      xv(1) = x(ix1) + Stepx
!      yv(1) = y(iy1) + Stepy
!
!      Write(6,'(" Stepx, Stepy....:",1p,2E15.6)')Stepx, Stepy
!
!      Write(6,'(" Posició final del 1er vortex...:",1p,2E15.6)')xv(1), yv(1)
!      If(Nv.Gt.1)Then
!        det   = Dreal(Dxpsi2D(ix2,iy2))*Dimag(Dypsi2D(ix2,iy2)) -    &
!              Dreal(Dypsi2D(ix2,iy2))*Dimag(Dxpsi2D(ix2,iy2)) 
!        Stepx = (Dimag(psi2D(ix2,iy2))*Dreal(Dypsi2D(ix2,iy2)) -     &
!               Dreal(psi2D(ix2,iy2))*Dimag(Dypsi2D(ix2,iy2)))/det 
!        Stepy = (Dreal(psi2D(ix2,iy2))*Dimag(Dxpsi2D(ix2,iy2)) -     &
!               Dimag(psi2D(ix2,iy2))*Dreal(Dxpsi2D(ix2,iy2)))/det 
!        xv(2) = x(ix2) + Stepx
!        yv(2) = y(iy2) + Stepy
!
!        Write(6,'(" Stepx, Stepy....:",1p,2E15.6)')Stepx, Stepy
! 
!        Write(6,'(" Posició final del 2on vortex...:",1p,2E15.6)')xv(2), yv(2)
!
!        Write(6,'("Posicio final dels vortex.....:  (1)-->:",1p,2E15.6,"  (2)-->:",2E15.6)')xv(1),yv(1), xv(2), yv(2)
!        Write(9,'(1p,4E15.6)')xv(1),yv(1), xv(2), yv(2)
!      Else  
!        Write(9,'(1p,4E15.6)')xv(1),yv(1)
!      EndIF
!    Else
!!
!!  Here we will compute the vortex position as a function of z
!!
!      Allocate(xv(2))      
!      Allocate(yv(2))
!      ixmin = (xlv + xmax)/hx + 1.5 
!      ixmax = (xrv + xmax)/hx + 1.5 
!      iymin = (ylv + ymax)/hy + 1.5 
!      iymax = (yrv + ymax)/hy + 1.5 
!      izmin = (zv_min + zmax0)/hz0 + 1.5 
!      izmax = (zv_max + zmax0)/hz0 + 1.5 
!      Write(6,'("izmin,izmzx,zv_min,zvmzx...:",2I8,1p,2E15.6)')izmin,izmax,zv_min,zv_max
!      Do iz = izmin, izmax
!        Psi02D = Psi0(:,:,iz)
!        Call INTERxy(nx,ny,x,y,Psi2D,nx0,ny0,hx0,hy0,Psi02D,K,KM1)
!        den2D = Max(Abs(psi2D)**2,1.d-99)
!        Write(6,'(" We write the value of MinLoc(den2Di(ixmin:ixmax,iymin:iymax)....:")')
!        Write(6,*)MinLoc(den2D(ixmin:ixmax,iymin:iymax))
!        I_MinLoc = MinLoc(den2D(ixmin:ixmax,iymin:iymax))
!        Write(6,'(" I_MinLoc(1,2)....:",2I10)')I_MinLoC
!        Write(6,'(" We write the value of MinVal(den2Di(ixmin:ixmax,iymin:iymax)....:")')
!        Write(6,*)MinVal(den2D(ixmin:ixmax,iymin:iymax))
!
!        ix1 = (I_MinLoc(1) - 1) + ixmin
!        iy1 = (I_MinLoc(2) - 1) + iymin
!
!        Write(6,'("Index del vortex:ix1,iy1............:",2I10)')ix1,iy1
!        Write(6,'("Posició del vortex:x(ix1),y(iy1)....:",1p,2E15.6)')x(ix1),y(iy1)
!
!        Write(6,'("den2D(ix1, iy1)...:",1p,E15.6)')den2D(ix1,iy1)
!!
!! Refinem la posició dels vortex
!!
!
!        det   = Dreal(Dxpsi2D(ix1,iy1))*Dimag(Dypsi2D(ix1,iy1)) -    &
!                Dreal(Dypsi2D(ix1,iy1))*Dimag(Dxpsi2D(ix1,iy1)) 
!        Stepx = (Dimag(psi2D(ix1,iy1))*Dreal(Dypsi2D(ix1,iy1)) -     &
!                 Dreal(psi2D(ix1,iy1))*Dimag(Dypsi2D(ix1,iy1)))/det 
!        Stepy = (Dreal(psi2D(ix1,iy1))*Dimag(Dxpsi2D(ix1,iy1)) -     &
!                 Dimag(psi2D(ix1,iy1))*Dreal(Dxpsi2D(ix1,iy1)))/det 
!        xv(1) = x(ix1) + Stepx
!        yv(1) = y(iy1) + Stepy
!
!        Write(6,'(" Stepx, Stepy....:",1p,2E15.6)')Stepx, Stepy
!
!        Write(6,'(" Posició final del 1er vortex...:",1p,2E15.6)')xv(1), yv(1)
!        Write(9,'(1p,4E15.6)')z0(iz),xv(1),yv(1)
!      EndDo  
!    EndIf        
!  EndIF
!
!  Pi=4.0d0*Datan(1.0d0)
!  TwoPi=2.0d0*Pi
!
!  ixi = (xi+xmax)/hx + 1.5
!  ixf = (xf+xmax)/hx + 1.5
!  iyi = (yi+ymax)/hy + 1.5
!  iyf = (yf+ymax)/hy + 1.5
!
!  write(6,'("ixi, ixf, iyi,iyf....:",4I5)')ixi,ixf,iyi,iyf
!  write(6,'("xi, xf, yi,yf....:",1p,4E15.6)')x(ixi),x(ixf),y(iyi),y(iyf)
!
!!
!! Calculo de la circulacion
!!
!
!  ynorm1 = 0.d0
!  ynorm2 = 0.d0
!  nyi = iyf - iyi +1
!  Allocate(Wy(nyi))
!  nxi = ixf - ixi +1
!  Allocate(Wx(nxi))
!  If(nyi.Lt.npi)Then
!    npi = 2*(nyi/2)
!    Write(6,'("I had changed npi..:",I4)')npi
!  Endif  
!  If(nxi.Lt.npi)Then
!    npi = 2*(nxi/2)      
!    Write(6,'("I had changed npi..:",I4)')npi
!  Endif  
!  Call Simps(npi,2,hx,x,Wx,nxi)
!  Call Simps(npi,2,hy,y,Wy,nyi)
!  j=0
!  Do iy=iyi,iyf
!    j = j +1
!!    if(iy.Gt.1.And.iy.Lt.iyf)Then
!      aux1 =  Dimag(Dypsi2D(ixi,iy)/Psi2D(ixi,iy))
!      aux2 =  Dimag(Dypsi2D(ixf,iy)/Psi2D(ixf,iy))
!!    Else
!!      aux1 =  0.5*Dimag(Dypsi2D(ixi,iy)/Psi2D(ixi,iy))
!!      aux2 =  0.5*Dimag(Dypsi2D(ixf,iy)/Psi2D(ixf,iy))
!!    Endif
!    ynorm1 = ynorm1 + Aux1*Wy(j)
!    ynorm2 = ynorm2 + Aux2*Wy(j)
!  EndDo
!
!  xnorm1 = 0.d0
!  xnorm2 = 0.d0
!  j = 0
!  Do ix=ixi,ixf
!    j = j +1
!!    If(ix.Gt.1.And.ix.Lt.ixf)Then
!      aux1 =  Dimag(Dxpsi2D(ix,iyi)/Psi2D(ix,iyi))
!      aux2 =  Dimag(Dxpsi2D(ix,iyf)/Psi2D(ix,iyf))
!!    Else
!!      aux1 =  0.5*Dimag(Dxpsi2D(ix,iyi)/Psi2D(ix,iyi))
!!      aux2 =  0.5*Dimag(Dxpsi2D(ix,iyf)/Psi2D(ix,iyf))
!!    Endif
!    xnorm1 = xnorm1 + aux1*Wx(j)
!    xnorm2 = xnorm2 + aux2*Wx(j)
!  EndDo
!  xnorm1=xnorm1*hx
!  xnorm2=xnorm2*hx
!
!  ynorm1=ynorm1*hy
!  ynorm2=ynorm2*hy
!
!  Write(6,'("xnorm1,xnorm2...:",1p,2E15.6)')xnorm1,xnorm2
!  Write(6,'("ynorm1,ynorm2...:",1p,2E15.6)')ynorm1,ynorm2
!
!  xnorm=-xnorm1+xnorm2
!  ynorm=ynorm1-ynorm2
!
!  circ = (xnorm + ynorm)/Twopi
!
!  Write(6,'("Circulacion.....:",1p,E15.6)')circ
!!
!!     Here we will compute <L_axis>
!!
!      If(Lvortex)Then
!        Open(Unit=13,File=File13)
!        xcm=0.0d0
!        ycm=0.0d0
!        zcm=0.0d0
!        Do ix=1, nx0
!          xcm = xcm + Sum(den0(ix,:,:)*x0(ix))
!        EndDo 
!        Do iy=1, ny0
!          ycm = ycm + Sum(den0(:,iy,:)*y0(iy))
!        EndDo 
!        Do iz=1, nz0
!          zcm = zcm + Sum(den0(:,:,iz)*z0(iz))
!        EndDo 
!        xcm = xcm*dxyz0/n4
!        ycm = ycm*dxyz0/n4
!        zcm = zcm*dxyz0/n4
!        Write(6,'("C.M. position: X,Y,Z..:",1p,3E15.6)')xcm,ycm,zcm
!        caux = (0.d0, 0.d0)      
!        If(Vortex_axis.Eq.'Z')Then
!          Write(13,'("#  <L_z> ")')      
!          Call derivnD(1,nnn,hx0,1,psi0,sto1c,Icon)
!          Call derivnD(1,nnn,hy0,2,psi0,sto2c,Icon)
!          Do iz=1, nz0
!            Do iy=1, ny0
!              Do ix=1, nx0
!                caux = caux + Ci*Conjg(Psi0(ix,iy,iz))*                  &
!                ((y0(iy)-ycm)*sto1c(ix,iy,iz) - (x0(ix)-xcm)*sto2c(ix,iy,iz)) 
!              EndDo
!            EndDo
!          EndDo
!          caux = caux*dxyz0
!          Write(6,'(T5,"<L_z>.........................:",8x,1p,E15.6," \hbar +-",E13.4," \hbar")')Real(caux)
!          Write(13,'(1p,E15.6)')Real(caux)
!        Endif
!        If(Vortex_axis.Eq.'Y')Then
!          Write(13,'("#  <L_y> ")')      
!          Call derivnD(1,nnn,hx0,1,psi0,sto1c,Icon)
!          Call derivnD(1,nnn,hz0,3,psi0,sto2c,Icon)
!          Do iz=1, nz0
!            Do iy=1, ny0
!              Do ix=1, nx0
!                caux = caux + Ci*Conjg(Psi0(ix,iy,iz))*                  &
!                ((x0(ix)-xcm)*sto2c(ix,iy,iz) - (z0(iz)-zcm)*sto1c(ix,iy,iz)) 
!              EndDo
!            EndDo
!          EndDo
!          caux = caux*dxyz0
!          Write(6,'(T5,"<L_y>.........................:",8x,1p,E15.6," \hbar +-",E13.4," \hbar")')Real(caux)
!          Write(13,'(1p,E15.6)')Real(caux)
!        Endif
!        If(Vortex_axis.Eq.'X')Then
!          Write(13,'("#  <L_y> ")')      
!          Call derivnD(1,nnn,hy0,2,psi0,sto1c,Icon)
!          Call derivnD(1,nnn,hz0,3,psi0,sto2c,Icon)
!          Do iz=1, nz0
!            Do iy=1, ny0
!              Do ix=1, nx0
!                caux = caux + Ci*Conjg(Psi0(ix,iy,iz))*                  &
!                ((z0(iz)-zcm)*sto1c(ix,iy,iz) - (y0(iy)-ycm)*sto2c(ix,iy,iz)) 
!              EndDo
!            EndDo
!          EndDo
!          caux = caux*dxyz0
!          Write(6,'(T5,"<L_x>.........................:",8x,1p,E15.6," \hbar +-",E13.4," \hbar")')Real(caux)
!          Write(13,'(1p,E15.6)')Real(caux)
!        Endif
!      Endif        
!
! We write a file with a format for the paraview program 
!
!EndIf
!Open(Unit=23,File=File23)
!Write(23,'("# vtk DataFile Version 4.2.0")')
!Write(23,'("stream")')
!Write(23,'("ASCII")')
!Write(23,'("DATASET RECTILINEAR_GRID")')
!Write(23,'("DIMENSIONS", 3I7)')nx, ny, 1
!Write(23,'("X_COORDINATES", I7," float")')nx
!Write(23,'(1p,(E18.10))')(x(ix), ix=1,nx)
!Write(23,'("Y_COORDINATES", I7," float")')ny
!Write(23,'(1p,(E18.10))')(y(iy), iy=1,ny)
!Write(23,'("Z_COORDINATES", I7," float")')1
!Write(23,'(1p,(E18.10))')0.d0
!Write(23,'("POINT_DATA", I17)')(nx*ny*1)
!!
!! We write a file with a format for the paraview program 
!!
!Open(Unit=20,File=File20)
!Write(20,'("# vtk DataFile Version 4.2.0")')
!Write(20,'("# Density and velocity field in the z=0 plane")')
!Write(20,'("ASCII")')
!Write(20,'("DATASET STRUCTURED_GRID")')
!Write(20,'("DIMENSIONS", 3I7)')nx, ny, 1
!Write(20,'("POINTS", I17," float")')(nx*ny*1)
!Do iy=1, ny
!  Do ix=1, nx
!    Write(20,'(1p,3E18.10)')x(ix),y(iy),0.d0
!  EndDo  
!EndDo
!Write(20,'("POINT_DATA", I17)')(nx*ny*1)
!Write(20,'("SCALARS dens float 1")')
!Write(20,'("LOOKUP_TABLE default")')
!Do iy=1, ny
!  Do ix=1, nx
!    Write(20,'(1p,E18.10)')den2D(ix,iy)
!  EndDo  
!EndDo
!Write(20,'("VECTORS veloc float")')
!
!Write(23,'("SCALARS Xvelocity float 1")')
!Write(23,'("LOOKUP_TABLE default")')
!Do iy=1, ny
!  Do ix=1,nx
!    If(den2D(ix,iy).Gt.epsrho)Then
!!      Stox = fac*dimag(Dxpsi2D(ix,iy)/psi2D(ix,iy))
!!      Stoy = fac*dimag(Dypsi2D(ix,iy)/psi2D(ix,iy))
!      Stox = fac*dimag(Dxpsi2D(ix,iy)/psi2D(ix,iy))*den2D(ix,iy)
!      Stoy = fac*dimag(Dypsi2D(ix,iy)/psi2D(ix,iy))*den2D(ix,iy)
!      If(Stox.Gt.vxmax)vxmax=Stox
!      If(Stoy.Gt.vymax)vymax=Stoy
!      If(Stox.Lt.vxmin)vxmin=Stox
!      If(Stoy.Lt.vymin)vymin=Stoy
!      Write(20,'(1p,3E18.10)')Stox, Stoy,0.d0
!      Write(23,'(1p,3E18.10)')Stox
!    Else  
!      Write(20,'(1p,3E18.10)')0.0d0,0.0d0,0.d0
!      Write(23,'(1p,3E18.10)')0.0d0
!    Endif
!  EndDo
!EndDo
!Close(20)
!Write(23,'("SCALARS Yvelocity float 1")')
!Write(23,'("LOOKUP_TABLE default")')
!Do iy=1, ny
!  Do ix=1,nx
!    If(den2D(ix,iy).Gt.epsrho)Then
!!      Stox = fac*dimag(Dxpsi2D(ix,iy)/psi2D(ix,iy))
!!      Stoy = fac*dimag(Dypsi2D(ix,iy)/psi2D(ix,iy))
!      Stox = fac*dimag(Dxpsi2D(ix,iy)/psi2D(ix,iy))*den2D(ix,iy)
!      Stoy = fac*dimag(Dypsi2D(ix,iy)/psi2D(ix,iy))*den2D(ix,iy)
!      If(Stox.Gt.vxmax)vxmax=Stox
!      If(Stoy.Gt.vymax)vymax=Stoy
!      If(Stox.Lt.vxmin)vxmin=Stox
!      If(Stoy.Lt.vymin)vymin=Stoy
!      Write(23,'(1p,3E18.10)')Stoy
!    Else  
!      Write(23,'(1p,3E18.10)')0.0d0
!    Endif
!  EndDo
!EndDo
!Write(23,'("VECTORS VecVelocity float")')
!Do iy=1, ny
!  Do ix=1,nx
!    If(den2D(ix,iy).Gt.epsrho)Then
!!      Stox = fac*dimag(Dxpsi2D(ix,iy)/psi2D(ix,iy))
!!      Stoy = fac*dimag(Dypsi2D(ix,iy)/psi2D(ix,iy))
!      Stox = fac*dimag(Dxpsi2D(ix,iy)/psi2D(ix,iy))*den2D(ix,iy)
!      Stoy = fac*dimag(Dypsi2D(ix,iy)/psi2D(ix,iy))*den2D(ix,iy)
!      If(Stox.Gt.vxmax)vxmax=Stox
!      If(Stoy.Gt.vymax)vymax=Stoy
!      If(Stox.Lt.vxmin)vxmin=Stox
!      If(Stoy.Lt.vymin)vymin=Stoy
!      Write(23,'(1p,3E18.10)')Stox,Stoy,0.d0
!    Else  
!      Write(23,'(1p,3E18.10)')0.0d0,0.d0,0.d0
!    Endif
!  EndDo
!EndDo
!Close(23)
!vxmax = -1.d10
!vxmin =  1.d10
!vymax = -1.d10
!vymin =  1.d10
!Open(Unit=8,File=File8)
!Do ix=1, nx
!  Do iy=1,ny
!    If(den2D(ix,iy).Gt.epsrho)Then
!!      Stox = fac*dimag(Dxpsi2D(ix,iy)/psi2D(ix,iy))
!!      Stoy = fac*dimag(Dypsi2D(ix,iy)/psi2D(ix,iy))
!      Stox = fac*dimag(Dxpsi2D(ix,iy)/psi2D(ix,iy))*den2D(ix,iy)
!      Stoy = fac*dimag(Dypsi2D(ix,iy)/psi2D(ix,iy))*den2D(ix,iy)
!      If(Stox.Gt.vxmax)vxmax=Stox
!      If(Stoy.Gt.vymax)vymax=Stoy
!      If(Stox.Lt.vxmin)vxmin=Stox
!      If(Stoy.Lt.vymin)vymin=Stoy
!      Write(8,'(1p,5E18.10)')x(ix), y(iy), Stox, Stoy, den2D(ix,iy)
!    Else
!      Write(8,'(1p,5E18.10)')x(ix),y(iy),0.0d0,0.0d0,den2D(ix,iy)
!    Endif
!  EndDo
!  Write(8,*)
!EndDo
!Close(8)
!Write(6,'("VxMin, VxMax....:",1p,2E15.6)')vxmin,vxmax
!Write(6,'("VyMin, VyMax....:",1p,2E15.6)')vymin,vymax
!!Close(6)
!
!If(.Not.L2Dplot)Stop
!Deallocate(psi2D); Deallocate(den2D);Deallocate(psi02D)
!Allocate(psi2D(nx,nz)); Allocate(den2D(nx,nz));Allocate(psi02D(nx0,nz0)  )
!
!!
!!       We store the density corresponding to the y=0 plane
!!
!
!Psi02D = Psi0(:,ny0/2+1,:)
!
!Call INTERxy(nx,nz,x,z,Psi2D,nx0,nz0,hx0,hz0,Psi02D,K,KM1)
!
!den2D = Max(Abs(psi2D)**2,1.d-99)
!    Write(6,'(" We write the value of MaxLoc(den2Di(:,:)(x,0,z)....:")')
!    Write(6,*)MaxLoc(den2D(:,:))
!    Write(6,'(" We write the value of MaxVal(den2Di(:,:)(x,0,z)....:")')
!    Write(6,*)MaxVal(den2D(:,:))
!
!Open(Unit=1,File="den-o-z.dat")
!Do iz=1,nz; Write(1,'(1p,2E15.6)')z(iz),den2D(nx/2+1,iz); EndDo
!Close (1)
!
!!
!! We write a file with a format for the paraview program 
!!
!Open(Unit=21,File=File21)
!Write(21,'("# vtk DataFile Version 4.2.0")')
!Write(21,'("# Density at y=0 plane")')
!Write(21,'("ASCII")')
!Write(21,'("DATASET STRUCTURED_GRID")')
!Write(21,'("DIMENSIONS", 3I7)')nx, 1, nz
!Write(21,'("POINTS", I17," float")')(nx*1*nz)
!Do iz=1, nz
!  Do ix=1, nx
!    Write(21,'(1p,3E18.10)')x(ix),0.d0,z(iz)
!  EndDo  
!EndDo
!Write(21,'("POINT_DATA", I17)')(nx*1*nz)
!Write(21,'("SCALARS dens float 1")')
!Write(21,'("LOOKUP_TABLE default")')
!Do iz=1, nz
!  Do ix=1, nx
!    Write(21,'(1p,E18.10)')den2D(ix,iz)
!  EndDo  
!EndDo
!Close(21)
!
!Open(Unit=11,File=File11)
!Do ix=1, nx
!  Do iz=1,nz
!    Write(11,'(1p,3E18.10)')x(ix),z(iz),den2D(ix,iz)
!  EndDo
!  Write(11,*)
!EndDo
!Close(11)
!
!Deallocate(psi2D); Deallocate(den2D);Deallocate(psi02D)
!Allocate(psi2D(ny,nz)); Allocate(den2D(ny,nz));Allocate(psi02D(ny0,nz0)  )
!
!!
!!       We store the density corresponding to the x=0 plane
!!
!
!Psi02D = Psi0(nx0/2+1,:,:)
!
!Call INTERxy(ny,nz,y,z,Psi2D,ny0,nz0,hy0,hz0,Psi02D,K,KM1)
!
!den2D = Max(Abs(psi2D)**2,1.d-99)
!    Write(6,'(" We write the value of MaxLoc(den2Di(:,:)(0,y,z)....:")')
!    Write(6,*)MaxLoc(den2D(:,:))
!    Write(6,'(" We write the value of MaxVal(den2Di(:,:)(0,y,z)....:")')
!    Write(6,*)MaxVal(den2D(:,:))
!!Close(6)
!
!!
!! We write a file with a format for the paraview program 
!!
!Open(Unit=22,File=File22)
!Write(22,'("# vtk DataFile Version 4.2.0")')
!Write(22,'("# Density in the x=0 plane")')
!Write(22,'("ASCII")')
!Write(22,'("DATASET STRUCTURED_GRID")')
!Write(22,'("DIMENSIONS", 3I7)')1, ny, nz
!Write(22,'("POINTS", I17," float")')(1*ny*nz)
!Do iz=1, nz
!  Do iy=1, ny
!    Write(22,'(1p,3E18.10)')0.d0,y(iy),z(iz)
!  EndDo  
!EndDo
!Write(22,'("POINT_DATA", I17)')(1*ny*nz)
!Write(22,'("SCALARS dens float 1")')
!Write(22,'("LOOKUP_TABLE default")')
!Do iz=1, nz
!  Do iy=1, ny
!    Write(22,'(1p,E18.10)')den2D(iy,iz)
!  EndDo  
!EndDo
!Close(22)
!
!Open(Unit=12,File=File12)
!Do iy=1, ny
!  Do iz=1,nz
!    Write(12,'(1p,3E18.10)')y(iy),z(iz),den2D(iy,iz)
!  EndDo
!  Write(12,*)
!EndDo
!Close(12)
If(L3Dplot)Then
  Open(Unit=41,File=File41, buffered='yes')      
  If(L_3D_filter)Then
!
!  Calcularem la posició del centre de masses
!
    xcm=0.0d0
    ycm=0.0d0
    zcm=0.0d0
    Do ix=1, nx0
      xcm = xcm + Sum(den0(ix,:,:)*x0(ix))
    EndDo 
    Do iy=1, ny0
      ycm = ycm + Sum(den0(:,iy,:)*y0(iy))
    EndDo 
    Do iz=1, nz0
      zcm = zcm + Sum(den0(:,:,iz)*z0(iz))
    EndDo 
    xcm = xcm*dxyz0/n4
    ycm = ycm*dxyz0/n4
    zcm = zcm*dxyz0/n4
!    Write(6,'("Posició del C.M....: xcm, ycm, zcm...",1p,3E15.6)')xcm,ycm,zcm
    x_ini = xcm + x_N_filter; ix_ini = Nint((x_ini + xmax0)/hx0) + 1; ix_ini = max(ix_ini,1)
    x_fin = xcm + x_P_filter; ix_fin = Nint((x_fin + xmax0)/hx0) + 1; ix_fin = min(ix_fin,nx0)          
    y_ini = ycm + y_N_filter; iy_ini = Nint((y_ini + ymax0)/hy0) + 1; iy_ini = max(iy_ini,1)
    y_fin = ycm + y_P_filter; iy_fin = Nint((y_fin + ymax0)/hy0) + 1; iy_fin = min(iy_fin,ny0)          
    z_ini = zcm + z_N_filter; iz_ini = Nint((z_ini + zmax0)/hz0) + 1; iz_ini = max(iz_ini,1)
    z_fin = zcm + z_P_filter; iz_fin = Nint((z_fin + zmax0)/hz0) + 1; iz_fin = min(iz_fin,nz0)
  Else
    ix_ini = 1; ix_fin = nx0 
    iy_ini = 1; iy_fin = ny0 
    iz_ini = 1; iz_fin = nz0 
  EndIf        
  nxx = ix_fin - ix_ini +1
  nyy = iy_fin - iy_ini +1
  nzz = iz_fin - iz_ini +1
!  Write(6,'("Dades del filtre en direccio X: ix_ini, ix_fin, nxx...:",3I8)')ix_ini, ix_fin, nxx
!  Write(6,'("Dades del filtre en direccio Y: iy_ini, iy_fin, nyy...:",3I8)')iy_ini, iy_fin, nyy
!  Write(6,'("Dades del filtre en direccio Z: iz_ini, iz_fin, nzz...:",3I8)')iz_ini, iz_fin, nzz
!
! Cabecera para las coordenadas
!
  Write(41,'("# vtk DataFile Version 4.2.0")')
  Write(41,'("# 3D Density ")')
  Write(41,'("ASCII")')
  Write(41,'("DATASET STRUCTURED_GRID")')
  Write(41,'("DIMENSIONS", 3I7)') nxx, nyy, nzz
  Write(41,'("POINTS", I17," float")') (nxx*nyy*nzz)
!
! Imprimir coordenadas
!
  DO iz = iz_ini, iz_fin
    DO iy = iy_ini, iy_fin
      DO ix = ix_ini, ix_fin
        Write(41,'(1p,3E18.10)') x0(ix),y0(iy),z0(iz)
      END DO
    END DO  
  END DO

! Cabecera para los datos

  Write(41,'("POINT_DATA", I17)') (nxx*nyy*nzz)
  Write(41,'("SCALARS dens float 1")')
  Write(41,'("LOOKUP_TABLE default")')

! Bucle para los datos
  DO iz = iz_ini, iz_fin
    DO iy = iy_ini, iy_fin
      DO ix = ix_ini, ix_fin
        Write(41,'(1p,E18.10)') den0(ix,iy,iz)
      END DO
    END DO  
  END DO
  Close(41)
Endif
Stop
End
