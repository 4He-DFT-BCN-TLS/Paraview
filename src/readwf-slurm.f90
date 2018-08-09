!Use Para_DerivnD
Implicit Real*8(A-H,O-Z)
Complex   (Kind=8), Allocatable :: psi0(:,:,:), psi02D(:,:), psi2D(:,:), Dxpsi2D(:,:), Dypsi2D(:,:)
complex(kind = 8), ALLOCATABLE	:: invar(:)			! Lambda (internal variables)
Complex   (Kind=8), Allocatable :: Sto1c(:,:,:), Sto2c(:,:,:),Sto3c(:,:,:)
Complex   (Kind=8), Allocatable :: DDxpsi2D(:,:), DDypsi2D(:,:), DDxypsi2D(:,:)
Real      (Kind=8), Allocatable :: den0(:,:,:), x0(:), y0(:), z0(:), den2D(:,:),x(:), y(:), z(:), Wx(:), Wy(:), &  
                                   x0v(:), y0v(:), xv(:), yv(:), rimp(:,:), vimp(:,:)
                               
Integer   (Kind=4) :: nn(2), Icon=13, nnn(3), i_minloc(2)
integer (kind=4)	:: denmode = 42
integer (kind = 4)	:: ninvar, nimp			! Number of components of lambda
Complex   (Kind=8) :: ci=(0.d0,1.d0), caux, caux1
Character (len =1) :: cchar="#", Vortex_axis='Z'
Logical            :: limp, L2Dplot=.false.,Ldensity=.true., Ldynamic=.false., L3Dplot=.false.,                 &
                      Lvortex=.false., Lvortex_position=.false.,Lxyv(2),Lreal_wave_function=.false.,            &
                      Lvortex_position_as_a_function_of_z=.false.,L_3D_filter=.false.
Character (len=80) :: File5='readwf.dat', File6='readwf.res', File31='position.dat', File32='velocity.dat'
Character (len=80) :: File7='denxy.dat', File8='current.out',File10='den.inp',File9='vortex_position.out',      &
                      File11='denxz.dat',File12='denyz.dat',File13='angular_momentum.dat', File15='Q2.out', nimpFile='nimp.dat'
Character (len=80) :: File20='denxy.vtk',File21='denxz.vtk',File22='denyz.vtk',File23='current.vtk',File41='denxyz.vtk'
Data nx/436/, ny/436/, nz/436/, hx/-1.d0/, hy/-1.d0/, hz/-1.d0/, npd/13/,Km1/4/, ndmax/2/, nthread/4/, npi/4/
Data fac/158.66624d0/,epsrho/1.d-6/,drop_radius/15.6d0/  !drop_radius = 2.2*N_He**(1/3)/Sqrt(2)
Data xi/-40.d0/,xf/40.d0/,yi/-40.d0/,yf/40.d0/,xlv/-20.d0/,xrv/+20.d0/, ylv/-20.d0/, yrv/+20.d0/, nv/2/
Data zv_min/-20.d0/, zv_max/20.d0/
Data X_P_filter/100.d0/,Y_P_filter/100.d0/,Z_P_filter/100.d0/, X_N_filter/-100.d0/,Y_N_filter/-100.d0/,Z_N_filter/-100.d0/

Namelist/Input/File4,File6,File7,File8,File10,File11,File12,File13,L2Dplot,Ldynamic,xlv,xrv,ylv,yrv,nv,Lvortex_position,   &
File20,File21,File22,File23,nx,ny,nz,hx,hy,hz,npd,npi,Km1,icon,epsrho,xi,xf,yi,yf,nthread,Vortex_axis,Lvortex,File9,       &
Lreal_wave_function,drop_radius, File15, File31, File32, L3Dplot, File41,Lvortex_position_as_a_function_of_z,zv_min,zv_max,&
L_3D_filter,X_P_filter,Y_P_filter,Z_P_filter,X_N_filter,Y_N_filter,Z_N_filter,denmode,nimpFile
read(5,nml=Input)
K=npd    ! Km1 will be the number of derivatives for the Taylor expansion

Open(31,File=File31, buffered='yes')

!! HERE STARTS THE WAVE FUNCTION READING
Open(10,File=File10, Status='Old', buffered='yes')
call titols(10,cchar,isalto)
select case (denmode)
	case (1)
		read(10,*) xmax0,ymax0,zmax0,hx0,hy0,hz0,nx0,ny0,nz0,limp,ximp,yimp,zimp
		Allocate(x0(nx0))
		Allocate(y0(ny0))
		Allocate(z0(nz0))
		Allocate(den0(nx0,ny0,nz0))
		Allocate(psi0(nx0,ny0,nz0))
		read(10,*) den0
		psi0 = sqrt(den0)
	case (2)
		read(10,*) xmax0,ymax0,zmax0,hx0,hy0,hz0,nx0,ny0,nz0,limp,ximp,yimp,zimp
		Allocate(x0(nx0))
		Allocate(y0(ny0))
		Allocate(z0(nz0))
		Allocate(den0(nx0,ny0,nz0))
		Allocate(psi0(nx0,ny0,nz0))
		read(10,*) psi0
		den0 = Conjg(psi0) * psi0
	case (3)
		read(10,*) xmax0,ymax0,zmax0,hx0,hy0,hz0,nx0,ny0,nz0,nimp
		write(*,*) "Number of impurities: ", nimp
		Allocate(rimp(nimp,3))
		Allocate(vimp(nimp,3))
		Allocate(x0(nx0))
		Allocate(y0(ny0))
		Allocate(z0(nz0))
		Allocate(den0(nx0,ny0,nz0))
		Allocate(psi0(nx0,ny0,nz0))
		read(10,*) rimp
		read(10,*) vimp
		read(10,*) psi0
		den0 = Conjg(psi0) * psi0
		do ix=1,nimp
			Write(31,'(F9.5,1X,F9.5,1X,F9.5)') rimp(ix,1), rimp(ix,2), rimp(ix,3)
		enddo
		Close(31)
		Open(314,File=nimpFile, buffered='yes')
		write(314,*) nimp
		close(314)
	case (4)
		read(10,*) xmax0,ymax0,zmax0,hx0,hy0,hz0,nx0,ny0,nz0,ximp,yimp,zimp,vximp,vyimp,vzimp,ninvar
		Allocate(x0(nx0))
		Allocate(y0(ny0))
		Allocate(z0(nz0))
		ALLOCATE(invar(ninvar))
		Allocate(den0(nx0,ny0,nz0))
		Allocate(psi0(nx0,ny0,nz0))
		read(10,*) invar
		read(10,*) psi0
		den0 = Conjg(psi0) * psi0
	case default
		write(*,*)
		write(*,*) "You have chosen a 'denmode' unequal to {1,2,3,4}. Please modify '2Dden.settings' and choose one of:"
		write(*,*)
		write(*,*) "denmode = 1:	STATIC helium REAL density"
		write(*,*) "denmode = 2:	STATIC helium COMPLEX density supporting vorticity"
		write(*,*) "denmode = 3:	DYNAMIC helium density"
		write(*,*) "denmode = 4:	DYNAMIC helium density + electronic state of impurity"
		call EXIT(10)
end select
close(10)

!! HERE STOPS THE WAVE FUNCTION READING


dxyz0=hx0*hy0*hz0
dxyz=dxyz0

n4 = (sum(den0)*dxyz0 + 0.5)
Do ix=1,nx0; x0(ix)=-xmax0+(ix-1)*hx0; EndDo
Do iy=1,ny0; y0(iy)=-ymax0+(iy-1)*hy0; EndDo
Do iz=1,nz0; z0(iz)=-zmax0+(iz-1)*hz0; EndDo

If(L3Dplot)Then
  Open(Unit=41,File=File41, buffered='yes')      
  If(L_3D_filter)Then
!
!  Calcularem la posici? del centre de masses
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

