      implicit real*8(a-h,o-z)
	  real*8  M_b, M_int, nu, nuh
	  	  
	  open(unit=7, file= 'tang_vectors.dat')
	  open(unit=8, file= 'acc.dat')
	  
C   Constants; units are pc, My, M-solar
	   G = 0.0045d0
	   a0 = 3.88d0
	   
	   pi = acos(-1.d0)
	  
c     data dwarf spheroidal

c--- Segue I

c      xlum = 346.7d0
c	  r_half = 28.d0
c	  R_mw = 28.d3
c	  xML = 10.d0

c--- Crater II

c      xlum = 1.6d5
c	  r_half = 1421d0   !1066d0 this is half-light radius
c	  R_mw = 120.d3
c      xM_mw = 6.5d10  
c	  xML = 2.d0
	  
c--- Andromeda XIX

c      xlum = 3.31d5
c	  r_half = 1799.d0
c	  R_mw = 116.d3   
c	  xM_mw = 2.d0 * xM_mw  !need baryonic mass of Andromeda, for now take twice M_mw
c	  xML = 2.d0
	  
	  
c--- Leo I (old)

c      xlum = 3.4d6
c	  r_half = 328.d0
c	  R_mw = 258.d3
c	  xML = 2.d0

c--- Leo I (new data)

c      xlum = 5.5d6
c	  r_half = 298.d0
c	  R_mw = 258.d3
c	  xML = 2.d0

c      Fornax (old, from our `17 paper)

c       xlum = 1.4d7   !2.04d7
c       r_half = 891d0   !792.d0
c	   xML = 1.d0
c       xM_mw = 6.5d10  
c	   R_mw = 149.d3
	   
c      Fornax (new, from Lelli et al.)

c       xlum = 2.04d7
c       r_half = 792.d0
c	   xML = 1.d0
c       xM_mw = 6.5d10  
c	   R_mw = 149.d3
	  
c--- Carina
c      xlum = 2.4d5
c	  r_half = 321.d0
c	  R_mw = 107.d3
c	  xML = 5.d0
	  

c      Leo II (new, from Lelli et al.)

c      xlum = 7.41d5
c      r_half = 219.d0
c      xML = 2.d0
c      xM_mw = 6.5d10  
c      R_mw = 236.d3	

	
c--- Sculptor
c      xlum = 1.4d6
c	  r_half = 347.d0
c	  R_mw = 86.d3
c	  xML = 4.d0
	  
c--- Draco
c      xlum = 2.7d5
c	  r_half = 261.d0
c	  R_mw = 76.d3
c	  xML = 20.d0
	  
c--- Sextans (old)
c      xlum = 4.1d5
c	  r_half = 909.d0
c	  R_mw = 89.d3
c	  xML = 6.d0

c--- Sextans (new)
c      xlum = 4.37d5
c	  r_half = 748.d0
c	  R_mw = 89.d3
c	  xML = 2.d0

c--- Ursa Minor
c      xlum = 2.0d5
c	  r_half = 373.d0
c	  R_mw = 78.d3
c	  xML = 31.d0

c--- NGC 2419 Globular Cluster
c      xlum = 9.d5
c	  r_half = 23.d0
c	  R_mw = 87.5d3
c	  xML = 1.d0

c--- NGC1052-DF2
c      xlum = 1.1d8
c	  r_half = 2.2d3
c	  xML = 2.d0
c     xM_mw = 1.d11  ! mass of NGC1052
c      R_mw = 80.d3

c--- Toy Model
c      xlum = 1.d7
c	  r_half = 500.d0
c	  xML = 1.d0
c	  R_mw = 200.d3
c	  xM_mw = 1.d15 !xM_mw ! * 1.d+4

	  
c--- Baryonic mass of dSph
c	  M_b = xML * xlum
	  M_b = 2.d8

      write(*,*)'Enter the desired LOG of the internal', 
     &' NEWTONIAN acceleration ratio at the half-mass radius'
	  read(*,*) accel_Rat_int_log ! log(a_N/a0)
	  accel_Rat_int = 10.d0**accel_Rat_int_log
	  r_half = sqrt(G * M_b / (2.d0 *a0*accel_Rat_int))
	  a_s = r_half * dsqrt(2.d0**(2./3.)-1.d0)     !  Plummer Model
	  
	  write(*,*) 'Half Mass Radius (pc) = ', r_half
	  
	  R_mw = 1.d6
	    
	  write(*,*)'Enter the desired', 
     &' ACTUAL external acceleration ratio (g_e/a_0)'
	  read(*,*) accel_Rat_ext
	  
      accel_Rat_ext_newt = (accel_Rat_ext**2/
     &                    (1.d0 + accel_Rat_ext))

      
	  
	  xM_mw = (accel_Rat_ext_newt * a0) * (R_mw)**2 / G
	  
	  write(*,555)xM_mw
 555   format('Host Mass =', 1p,e12.3,2x, 'Solar Masses')
	

c --- polar coordinates of the host galaxy: fixed
 	  
	  rh = R_mw
	  theta_h = 0.d0
	  
c --- initialize coordinates
	  rad_fact = 1.d0
	  
 17	  continue   ! loop that changes r
	  
	  r = r_half * rad_fact
	  theta = 0.d0
	  dtheta = 2*pi / (50 * rad_fact)
	  
c--- interior mass

      M_int = M_b * r**3/((r**2 + a_s**2)**(1.5d0))   ! Plummer Model
	  
 16   continue   ! start of loop over theta values
 
c ---  Cartesian coordinates

      x = r * cos(theta)
      y = r * sin(theta)
 
c--- internal accelerations
c--- Newtonian, x & y components & magnitude

      gni_x = (-G*M_int/r**2) * cos(theta)
      gni_y = (-G*M_int/r**2) * sin(theta)
	  gni = sqrt(gni_x**2 + gni_y**2)
	  
c--- Polar components

	  gni_r = gni_x * cos(theta) + gni_y * sin(theta)
	  gni_t = gni_x * (-sin(theta)) + gni_y * cos(theta)
	  
c---- Isolated MOND

C   Lelli Nu
      yi = gni/a0
	  nui = 1.d0/(1.d0 - exp(-sqrt(yi)))
c      Simple Nu
c      nui = (1.d0 + sqrt(1.d0 + 4.d0/yi))/2.d0

	  gi_iso_x = nui * gni_x
	  gi_iso_y = nui * gni_y

      gi_iso_r = nui * gni_r
	  gi_iso_t = nui * gni_t
	  
	  gi_iso = sqrt(gi_iso_x**2 + gi_iso_y**2)
	  
c--- external accelerations
c--- Newtonian, x & y components & magnitude

      re_x = x - R_mw
	  re_y = y 
	  
	  re = sqrt(re_x**2 + re_y**2)
	  re3 = re**3
	  factor_1 = G * xM_mw /re3
	  
      gne_x = -re_x * factor_1  
	  gne_y = -re_y * factor_1
	  gne = sqrt(gne_x**2 + gne_y**2)
	  
c--- Polar components

	  gne_r = gne_x * cos(theta) + gne_y * sin(theta)
	  gne_t = gne_x * (-sin(theta)) + gne_y * cos(theta)
	  
c--- MOND correction

c   Lelli Nu
      yh = gne/a0
	  nuh = 1.d0/(1.d0 - exp(-sqrt(yh)))
c    Simple Nu
c      nuh = (1.d0 + sqrt(1.d0 + 4.d0/yh))/2.d0	  

	  ge_x = nuh * gne_x
	  ge_y = nuh * gne_y
	  
	  ge_r = nuh * gne_r
	  ge_t = nuh * gne_t
	  
	  ge = sqrt(ge_x**2 + ge_y**2)
	  
c--- Total Newtonian accelerations

	  gntot_x = gni_x + gne_x
      gntot_y = gni_y + gne_y
      gntot = sqrt(gntot_x**2 + gntot_y**2)
	  
	  gntot_r = gntot_x * cos(theta) + gntot_y * sin(theta)
	  gntot_t = gntot_x * (-sin(theta)) + gntot_y * cos(theta)
	  
c--- External Field Effect gi = nu*gt - nu*ge

c   Lelli Nu
      yn = (gntot)/a0
	  nu = 1.d0/(1.d0 - exp(-sqrt(yn)))
c      Simple Nu
c      nu = (1.d0 + sqrt(1.d0 + 4.d0/yn))/2.d0

c--- Total Acceleration

      gtot_x = nu * gntot_x
	  gtot_y = nu * gntot_y
	  
	  gtot = sqrt(gtot_x**2 + gtot_y**2)
	  
	  gi_x = gtot_x - ge_x
	  gi_y = gtot_y - ge_y
	  
	  gi_r = gi_x * cos(theta) + gi_y * sin(theta)
	  gi_t = gi_x * (-sin(theta)) + gi_y * cos(theta)
	  
	  gi = sqrt(gi_x**2 + gi_y**2)
	  
	  
c---- Normalize all accelerations to a0
      gni_xa  = gni_x/a0				! Newtonian internal
	  gi_iso_xa    = gi_iso_x/a0		! MOND internal isolated
	  gnea  = gne/a0				! Newtonian external
	  ge_xa    = ge_x/a0				! MOND external
	  gntota = gntot/a0				! Newtonian total
      gtota = gtot/a0				! MOND total	  
	  gia    = gi/a0				! MOND internal (w/ EFE)
	  
	  gi_ra  = gi_r/a0
	  gi_ta  = gi_t/a0				! Tangential comp of integrated internal acc
	  ge_ra  = ge_r/a0
	  ge_ta  = ge_t/a0
	  gi_iso_ra = gi_iso_r/a0
	  gi_iso_ta = gi_iso_t/a0
	  gne_ra = gne_r/a0
	  gne_ta = gne_t/a0
	  
	  diff_x = gtot_x - gi_iso_x - ge_x
	  diff_y = gtot_y - gi_iso_y - ge_y
	  diff_r = diff_x * cos(theta) + diff_y * sin(theta)
	  diff_t = diff_x * (-sin(theta)) + diff_y * cos(theta)
	  
	  
 1978 format(9e14.6)
      write(8,1978)theta, gi_ra, gi_ta, ge_ra, ge_ta, gi_iso_ra
     &             , gi_iso_ta, diff_r, diff_t			!gi_xa,ge_xa,gni_xa,gi_iso_xa
	 
	  
 2024 format(4e14.6)
      write(7,2024) x, y, theta + (pi/2.0), 1000 * diff_t
	  write(7,2024) x, y, theta, 1000 * diff_r
	  
	  theta = theta + dtheta
	  
	  if(theta .le. 2*pi)then
	     go to 16
c	  else if(rad_fact .gt. 0.5)then
c	     rad_fact = rad_fact - 0.5
c		 go to 17
	  endif
	  
	  stop
	  end
	  
	  
	  
	  


 
	  

	  
