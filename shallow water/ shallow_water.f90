program shallow_water
  implicit none

  integer, parameter :: n = 100
  integer, parameter :: na = n + 1, nb = n, nt = 1800
  real(8), parameter :: km = 15000.0d0
  real(8), parameter :: dx=30000.0d0, dy=30000.0d0
  real(8), parameter :: dt=60.0d0, g=9.81d0
  real(8), parameter :: omega=7.2921150d-5
  real(8), parameter :: pi=3.141592653589793d0
  real(8), parameter :: phi0=35.0d0*pi/180.0d0
  real(8), parameter :: f=2.0d0*omega*sin(phi0)
  
  real(8), parameter :: tend=30.0d0*3600.0d0
  real(8), parameter :: z0=3000.0d0

  real(8) :: radius, depth_max
  real(8) :: zmin, vmax, vrmin
  

  write(*,'(a)') '=============================================='
  write(*,'(a)') ' SHALLOW-WATER MODEL: A-GRID / C-GRID'
  write(*,'(a)') '=============================================='

  radius = 300000.0d0
  depth_max = 150.0d0

  write(*,'(a)') 'R = 300 km, depth = 150 m'
  write(*,'(a)') 'Grid          Zmin(m)       Vmax(m/s)     vrmin(m/s)'

  call run_a_case(radius,depth_max,zmin,vmax,vrmin)
  write(*,'(a8,3(2x,es13.5))') 'A-grid',zmin,vmax,vrmin

  call run_c_case(radius,depth_max,zmin,vmax,vrmin)
  write(*,'(a8,3(2x,es13.5))') 'C-grid',zmin,vmax,vrmin

  write(*,'(a)') '=============================================='
  write(*,'(a)') 'Simulations completed.'

contains

  subroutine run_a_case(radius,depth_max,zmin,vmax,vrmin)
    implicit none
   
    real(8), intent(in) :: radius,depth_max
    real(8), intent(out) :: zmin,vmax,vrmin
    real(8) :: u(na,na),v(na,na),z(na,na)
    real(8) :: ut(na,na),vt(na,na),zs(na,na)
    real(8) :: k1u(na,na),k2u(na,na),k3u(na,na),k4u(na,na)
    real(8) :: k1v(na,na),k2v(na,na),k3v(na,na),k4v(na,na)
    real(8) :: time
    integer :: n,unit_out
    character(40) :: fname

    fname='unit3_Agrid.txt'
    open(newunit=unit_out,file=trim(fname),status='replace',action='write')
    write(unit_out,'(a)') 'time_h x_km y_km u v z'

    u=0.0d0
    v=0.0d0
    call calc_z_a(z,0.0d0,radius,depth_max)
    call output_a(unit_out,0.0d0,u,v,z)
    
    do n=0,nt-1
      time=dble(n)*dt
      call calc_z_a(zs,time,radius,depth_max)
      call rhs_a(u,v,zs,k1u,k1v)

      ut=u+0.5d0*dt*k1u
      vt=v+0.5d0*dt*k1v
      call calc_z_a(zs,time+0.5d0*dt,radius,depth_max)
      call rhs_a(ut,vt,zs,k2u,k2v)

      ut=u+0.5d0*dt*k2u
      vt=v+0.5d0*dt*k2v
      call calc_z_a(zs,time+0.5d0*dt,radius,depth_max)
      call rhs_a(ut,vt,zs,k3u,k3v)

      ut=u+dt*k3u
      vt=v+dt*k3v
      call calc_z_a(zs,time+dt,radius,depth_max)
      call rhs_a(ut,vt,zs,k4u,k4v)

      u=u+dt*(k1u+2.0d0*k2u+2.0d0*k3u+k4u)/6.0d0
      v=v+dt*(k1v+2.0d0*k2v+2.0d0*k3v+k4v)/6.0d0

      if (mod(n+1,120)==0) then
        call calc_z_a(z,dble(n+1)*dt,radius,depth_max)
        call output_a(unit_out,dble(n+1)*dt,u,v,z)
      end if
    end do
    close(unit_out)

    call calc_z_a(z,tend,radius,depth_max)
    call diag_a(u,v,z,zmin,vmax,vrmin)
    call output_profile_a(z)
   call output_2d_a(z)
  end subroutine run_a_case


    subroutine run_c_case(radius,depth_max,zmin,vmax,vrmin)
    implicit none
    
    real(8), intent(in) :: radius,depth_max
    real(8), intent(out) :: zmin,vmax,vrmin
    real(8) :: u(nb+1,nb),v(nb,nb+1),z(nb,nb)
    real(8) :: ut(nb+1,nb),vt(nb,nb+1),zs(nb,nb)
    real(8) :: k1u(nb+1,nb),k2u(nb+1,nb),k3u(nb+1,nb),k4u(nb+1,nb)
    real(8) :: k1v(nb,nb+1),k2v(nb,nb+1),k3v(nb,nb+1),k4v(nb,nb+1)
    integer :: n,unit_out
    real(8) :: time
    character(40) :: fname


    fname='unit3_Cgrid.txt'
    open(newunit=unit_out,file=trim(fname),status='replace',action='write')
    write(unit_out,'(a)') 'time_h x_km y_km u v z'

    u=0.0d0
    v=0.0d0
    call calc_z_c(z,0.0d0,radius,depth_max)
    call output_c(unit_out,0.0d0,u,v,z)

    do n=0,nt-1
      time=dble(n)*dt
      call calc_z_c(zs,time,radius,depth_max)
      call rhs_c(u,v,zs,k1u,k1v)

      ut=u+0.5d0*dt*k1u
      vt=v+0.5d0*dt*k1v
      call calc_z_c(zs,time+0.5d0*dt,radius,depth_max)
      call rhs_c(ut,vt,zs,k2u,k2v)

      ut=u+0.5d0*dt*k2u
      vt=v+0.5d0*dt*k2v
      call calc_z_c(zs,time+0.5d0*dt,radius,depth_max)
      call rhs_c(ut,vt,zs,k3u,k3v)

      ut=u+dt*k3u
      vt=v+dt*k3v
      call calc_z_c(zs,time+dt,radius,depth_max)
      call rhs_c(ut,vt,zs,k4u,k4v)

      u=u+dt*(k1u+2.0d0*k2u+2.0d0*k3u+k4u)/6.0d0
      v=v+dt*(k1v+2.0d0*k2v+2.0d0*k3v+k4v)/6.0d0

      if (mod(n+1,120)==0) then
        call calc_z_c(z,dble(n+1)*dt,radius,depth_max)
        call output_c(unit_out,dble(n+1)*dt,u,v,z)
      end if
    end do
    close(unit_out)

    call calc_z_c(z,tend,radius,depth_max)
    call diag_c(u,v,z,zmin,vmax,vrmin)
    call output_profile_c(z)
    call output_2d_c(z)
  end subroutine run_c_case


  subroutine calc_z_a(z,time,radius,depth_max)
    implicit none
    real(8), intent(out) :: z(na,na)
    real(8), intent(in) :: time,radius,depth_max
    real(8) :: x,y,r,tau,ramp,bowl
    integer :: i,j

    tau = min(1.0d0,max(0.0d0,time/tend))
    ramp = 3.0d0*tau**2 - 2.0d0*tau**3

    do j=1,na
      y = -1500000.0d0 + dble(j-1)*dy

      do i=1,na
        x = -1500000.0d0 + dble(i-1)*dx
        r = sqrt(x*x+y*y)

        bowl = 0.0d0
        if (r <= radius) then
          bowl = 1.0d0 - r/radius
        end if

        z(i,j) = z0-depth_max*ramp*bowl
      end do
    end do
  end subroutine calc_z_a


  subroutine calc_z_c(z,time,radius,depth_max)
    implicit none
    real(8), intent(out) :: z(nb,nb)
    real(8), intent(in) :: time,radius,depth_max
    real(8) :: x,y,r,tau,ramp,bowl
    integer :: i,j

    tau = min(1.0d0,max(0.0d0,time/tend))
    ramp = 3.0d0*tau**2 - 2.0d0*tau**3

    do j=1,nb
      y = dble(j-(nb+1)/2)*dy

      do i=1,nb
        x = dble(i-(nb+1)/2)*dx
        r = sqrt(x*x+y*y)

        bowl = 0.0d0
        if (r <= radius) then
          bowl = 1.0d0 - r/radius
        end if

        z(i,j) = z0-depth_max*ramp*bowl
      end do
    end do
  end subroutine calc_z_c

  subroutine rhs_a(u,v,z,du,dv)
    implicit none
    real(8), intent(in) :: u(na,na),v(na,na),z(na,na)
    real(8), intent(out) :: du(na,na),dv(na,na)
    real(8) :: dzdx,dzdy,lapu,lapv,advu,advv
    integer :: i,j

    du=0.0d0
    dv=0.0d0

    do j=2,na-1
      do i=2,na-1
        advu=upwind_a(u(i,j),u(i-1,j),u(i+1,j), &
          u(i,j-1),u(i,j+1),u(i,j),v(i,j))
        advv=upwind_a(v(i,j),v(i-1,j),v(i+1,j), &
          v(i,j-1),v(i,j+1),u(i,j),v(i,j))

        dzdx=(z(i+1,j)-z(i-1,j))/(2.0d0*dx)
        dzdy=(z(i,j+1)-z(i,j-1))/(2.0d0*dy)
        lapu=(u(i+1,j)-2.0d0*u(i,j)+u(i-1,j))/dx**2
        lapv=(v(i,j+1)-2.0d0*v(i,j)+v(i,j-1))/dy**2

        du(i,j)=advu-g*dzdx+f*v(i,j)+km*lapu
        dv(i,j)=advv-g*dzdy-f*u(i,j)+km*lapv
      end do
    end do

    do j=2,na-1
      du(1,j)=boundary_tendency(-u(1,j),u(1,j),u(2,j),dx)
      du(na,j)=boundary_tendency(u(na,j),u(na,j),u(na-1,j),dx)
      dv(1,j)=boundary_tendency(-u(1,j),v(1,j),v(2,j),dx)
      dv(na,j)=boundary_tendency(u(na,j),v(na,j),v(na-1,j),dx)
    end do

    do i=2,na-1
      du(i,1)=boundary_tendency(-v(i,1),u(i,1),u(i,2),dy)
      du(i,na)=boundary_tendency(v(i,na),u(i,na),u(i,na-1),dy)
      dv(i,1)=boundary_tendency(-v(i,1),v(i,1),v(i,2),dy)
      dv(i,na)=boundary_tendency(v(i,na),v(i,na),v(i,na-1),dy)
    end do

    du(1,1)=boundary_tendency(-u(1,1),u(1,1),u(2,1),dx)+ &
            boundary_tendency(-v(1,1),u(1,1),u(1,2),dy)
    du(na,1)=boundary_tendency(u(na,1),u(na,1),u(na-1,1),dx)+ &
             boundary_tendency(-v(na,1),u(na,1),u(na,2),dy)
    du(1,na)=boundary_tendency(-u(1,na),u(1,na),u(2,na),dx)+ &
             boundary_tendency(v(1,na),u(1,na),u(1,na-1),dy)
    du(na,na)=boundary_tendency(u(na,na),u(na,na),u(na-1,na),dx)+ &
              boundary_tendency(v(na,na),u(na,na),u(na,na-1),dy)

    dv(1,1)=boundary_tendency(-u(1,1),v(1,1),v(2,1),dx)+ &
            boundary_tendency(-v(1,1),v(1,1),v(1,2),dy)
    dv(na,1)=boundary_tendency(u(na,1),v(na,1),v(na-1,1),dx)+ &
             boundary_tendency(-v(na,1),v(na,1),v(na,2),dy)
    dv(1,na)=boundary_tendency(-u(1,na),v(1,na),v(2,na),dx)+ &
             boundary_tendency(v(1,na),v(1,na),v(1,na-1),dy)
    dv(na,na)=boundary_tendency(u(na,na),v(na,na),v(na-1,na),dx)+ &
              boundary_tendency(v(na,na),v(na,na),v(na,na-1),dy)
  end subroutine rhs_a


  real(8) function upwind_a(qc,qw,qe,qs,qn,a,b)
    implicit none
    real(8), intent(in) :: qc,qw,qe,qs,qn,a,b
    real(8) :: ap,am,bp,bm
    ap=max(a,0.0d0)
    am=min(a,0.0d0)
    bp=max(b,0.0d0)
    bm=min(b,0.0d0)
    upwind_a=-ap*(qc-qw)/dx-am*(qe-qc)/dx &
             -bp*(qc-qs)/dy-bm*(qn-qc)/dy
  end function upwind_a


  real(8) function boundary_tendency(un,qb,qin,dn)
    implicit none
    real(8), intent(in) :: un,qb,qin,dn
    real(8) :: un_plus,un_minus
    real(8), parameter :: qext=0.0d0
    un_plus=max(un,0.0d0)
    un_minus=min(un,0.0d0)
    boundary_tendency=-un_plus*(qb-qin)/dn+un_minus*(qb-qext)/dn
  end function boundary_tendency


  subroutine rhs_c(u,v,z,du,dv)
    implicit none
    real(8), intent(in) :: u(nb+1,nb),v(nb,nb+1),z(nb,nb)
    real(8), intent(out) :: du(nb+1,nb),dv(nb,nb+1)
    real(8) :: vint_u,uint_v,dzdx_u,dzdy_v,lapu,lapv
    real(8) :: ap,am,bp,bm,adv,un_e,un_w,un_n,un_s
    integer :: i,j

    du=0.0d0
    dv=0.0d0

    do j=2,nb-1
      do i=2,nb
        vint_u=0.25d0*(v(i-1,j)+v(i,j)+v(i-1,j+1)+v(i,j+1))
        ap=max(u(i,j),0.0d0)
        am=min(u(i,j),0.0d0)
        bp=max(vint_u,0.0d0)
        bm=min(vint_u,0.0d0)
        adv=-ap*(u(i,j)-u(i-1,j))/dx &
            -am*(u(i+1,j)-u(i,j))/dx &
            -bp*(u(i,j)-u(i,j-1))/dy &
            -bm*(u(i,j+1)-u(i,j))/dy
        dzdx_u=(z(i,j)-z(i-1,j))/dx
        lapu=(u(i+1,j)-2.0d0*u(i,j)+u(i-1,j))/dx**2
        du(i,j)=adv-g*dzdx_u+f*vint_u+km*lapu
      end do
    end do

    do j=2,nb
      do i=2,nb-1
        uint_v=0.25d0*(u(i,j-1)+u(i+1,j-1)+u(i,j)+u(i+1,j))
        ap=max(uint_v,0.0d0)
        am=min(uint_v,0.0d0)
        bp=max(v(i,j),0.0d0)
        bm=min(v(i,j),0.0d0)
        adv=-ap*(v(i,j)-v(i-1,j))/dx &
            -am*(v(i+1,j)-v(i,j))/dx &
            -bp*(v(i,j)-v(i,j-1))/dy &
            -bm*(v(i,j+1)-v(i,j))/dy
        dzdy_v=(z(i,j)-z(i,j-1))/dy
        lapv=(v(i,j+1)-2.0d0*v(i,j)+v(i,j-1))/dy**2
        dv(i,j)=adv-g*dzdy_v-f*uint_v+km*lapv
      end do
    end do

    do j=1,nb
      un_w=-u(1,j)
      du(1,j)=boundary_tendency(un_w,u(1,j),u(2,j),dx)
      un_e=u(nb+1,j)
      du(nb+1,j)=boundary_tendency(un_e,u(nb+1,j),u(nb,j),dx)
    end do

    do i=2,nb
      vint_u=0.25d0*(v(i-1,1)+v(i,1)+v(i-1,2)+v(i,2))
      du(i,1)=boundary_tendency(-vint_u,u(i,1),u(i,2),dy)
      vint_u=0.25d0*(v(i-1,nb)+v(i,nb)+v(i-1,nb+1)+v(i,nb+1))
      du(i,nb)=boundary_tendency(vint_u,u(i,nb),u(i,nb-1),dy)
    end do

    du(1,1)=boundary_tendency(-u(1,1),u(1,1),u(2,1),dx)+ &
            boundary_tendency(-v(1,1),u(1,1),u(1,2),dy)
    du(nb+1,1)=boundary_tendency(u(nb+1,1),u(nb+1,1),u(nb,1),dx)+ &
               boundary_tendency(-v(nb,1),u(nb+1,1),u(nb+1,2),dy)
    du(1,nb)=boundary_tendency(-u(1,nb),u(1,nb),u(2,nb),dx)+ &
             boundary_tendency(v(1,nb+1),u(1,nb),u(1,nb-1),dy)
    du(nb+1,nb)=boundary_tendency(u(nb+1,nb),u(nb+1,nb),u(nb,nb),dx)+ &
                boundary_tendency(v(nb,nb+1),u(nb+1,nb),u(nb+1,nb-1),dy)

    do i=1,nb
      dv(i,1)=boundary_tendency(-v(i,1),v(i,1),v(i,2),dy)
      dv(i,nb+1)=boundary_tendency(v(i,nb+1),v(i,nb+1),v(i,nb),dy)
    end do

    do j=2,nb
      uint_v=0.25d0*(u(1,j-1)+u(2,j-1)+u(1,j)+u(2,j))
      dv(1,j)=boundary_tendency(-uint_v,v(1,j),v(2,j),dx)
      uint_v=0.25d0*(u(nb,j-1)+u(nb+1,j-1)+u(nb,j)+u(nb+1,j))
      dv(nb,j)=boundary_tendency(uint_v,v(nb,j),v(nb-1,j),dx)
    end do

    dv(1,1)=boundary_tendency(-u(1,1),v(1,1),v(2,1),dx)+ &
            boundary_tendency(-v(1,1),v(1,1),v(1,2),dy)
    dv(nb,1)=boundary_tendency(u(nb+1,1),v(nb,1),v(nb-1,1),dx)+ &
             boundary_tendency(-v(nb,1),v(nb,1),v(nb,2),dy)
    dv(1,nb+1)=boundary_tendency(-u(1,nb),v(1,nb+1),v(2,nb+1),dx)+ &
               boundary_tendency(v(1,nb+1),v(1,nb+1),v(1,nb),dy)
    dv(nb,nb+1)=boundary_tendency(u(nb+1,nb),v(nb,nb+1),v(nb-1,nb+1),dx)+ &
                boundary_tendency(v(nb,nb+1),v(nb,nb+1),v(nb,nb),dy)
  end subroutine rhs_c

  subroutine output_a(unit_out,time,u,v,z)
    implicit none
    integer, intent(in) :: unit_out
    real(8), intent(in) :: time,u(na,na),v(na,na),z(na,na)
    integer :: i,j
    real(8) :: x,y

    do j=1,na
      y=dble(j-(na+1)/2)*dy/1000.0d0
      do i=1,na
        x=dble(i-(na+1)/2)*dx/1000.0d0
        write(unit_out,'(f7.2,1x,f9.2,1x,f9.2,1x,es14.6,1x,es14.6,1x,f10.3)') &
          time/3600.0d0,x,y,u(i,j),v(i,j),z(i,j)
      end do
    end do
    write(unit_out,*)
  end subroutine output_a


  subroutine output_c(unit_out,time,u,v,z)
    implicit none
    integer, intent(in) :: unit_out
    real(8), intent(in) :: time,u(nb+1,nb),v(nb,nb+1),z(nb,nb)
    integer :: i,j
    real(8) :: x,y,uc,vc

    do j=1,nb
      y = (dble(j-(nb+1)/2)*dy)/1000.0d0
      do i=1,nb
       x = (dble(i-(nb+1)/2)*dx)/1000.0d0
        uc=0.5d0*(u(i,j)+u(i+1,j))
        vc=0.5d0*(v(i,j)+v(i,j+1))
        write(unit_out,'(f7.2,1x,f9.2,1x,f9.2,1x,es14.6,1x,es14.6,1x,f10.3)') &
          time/3600.0d0,x,y,uc,vc,z(i,j)
      end do
    end do
    write(unit_out,*)
  end subroutine output_c


  subroutine output_profile_a(z)
    implicit none
    real(8), intent(in) :: z(na,na)
    integer :: i,j0,un
    real(8) :: x
    character(80) :: filename

     filename='Z_profile_Agrid.txt'
    open(newunit=un,file=trim(filename),status='replace',action='write')
    write(un,'(a)') 'x_km z_m'
    j0=(na+1)/2
    do i=1,na
      x=dble(i-(na+1)/2)*dx/1000.0d0
      write(un,'(f12.3,1x,f12.5)') x,z(i,j0)
    end do
    close(un)
  end subroutine output_profile_a


  subroutine output_profile_c(z)
    implicit none
    real(8), intent(in) :: z(nb,nb)
    integer :: i,j0,un
    real(8) :: x
    character(80) :: filename

     filename='Z_profile_Cgrid.txt'
    open(newunit=un,file=trim(filename),status='replace',action='write')
    write(un,'(a)') 'x_km z_m'
    j0=(nb+1)/2
    do i=1,nb
      x = (dble(i-(nb+1)/2)*dx)/1000.0d0
      write(un,'(f12.3,1x,f12.5)') x,z(i,j0)
    end do
    close(un)
  end subroutine output_profile_c


  subroutine output_2d_a(z)
    implicit none
    real(8), intent(in) :: z(na,na)
    integer :: i,j,un
    real(8) :: x,y
    character(80) :: filename

   filename='Z_final_Agrid.txt'
    open(newunit=un,file=trim(filename),status='replace',action='write')
    write(un,'(a)') 'x_km y_km z_m'
    do j=1,na
      y=dble(j-(na+1)/2)*dy/1000.0d0
      do i=1,na
        x=dble(i-(na+1)/2)*dx/1000.0d0
        write(un,'(f12.3,1x,f12.3,1x,f12.5)') x,y,z(i,j)
      end do
    end do
    close(un)
  end subroutine output_2d_a


  subroutine output_2d_c(z)
    implicit none
    real(8), intent(in) :: z(nb,nb)
    integer :: i,j,un
    real(8) :: x,y
    character(80) :: filename

    filename='Z_final_Cgrid.txt'
    open(newunit=un,file=trim(filename),status='replace',action='write')
    write(un,'(a)') 'x_km y_km z_m'
    do j=1,nb
      y = (dble(j-(nb+1)/2)*dy)/1000.0d0
      do i=1,nb
        x = (dble(i-(nb+1)/2)*dx)/1000.0d0
        write(un,'(f12.3,1x,f12.3,1x,f12.5)') x,y,z(i,j)
      end do
    end do
    close(un)
  end subroutine output_2d_c


  subroutine diag_a(u,v,z,zmin,vmax,vrmin)
    implicit none
    real(8), intent(in) :: u(na,na),v(na,na),z(na,na)
    real(8), intent(out) :: zmin,vmax,vrmin
    integer :: i,j
    real(8) :: x,y,r,speed,vr

    zmin=minval(z)
    vmax=0.0d0
    vrmin=huge(1.0d0)
    do j=1,na
      y=dble(j-(na+1)/2)*dy
      do i=1,na
        x=dble(i-(na+1)/2)*dx
        speed=sqrt(u(i,j)**2+v(i,j)**2)
        vmax=max(vmax,speed)
        r=sqrt(x*x+y*y)
        if (r>0.0d0) then
          vr=(x*u(i,j)+y*v(i,j))/r
          vrmin=min(vrmin,vr)
        end if
      end do
    end do
  end subroutine diag_a


  subroutine diag_c(u,v,z,zmin,vmax,vrmin)
    implicit none
    real(8), intent(in) :: u(nb+1,nb),v(nb,nb+1),z(nb,nb)
    real(8), intent(out) :: zmin,vmax,vrmin
    integer :: i,j
    real(8) :: x,y,r,uc,vc,speed,vr

    zmin=minval(z)
    vmax=0.0d0
    vrmin=huge(1.0d0)
    do j=1,nb
      y = dble(j-(nb+1)/2)*dy
      do i=1,nb
        x = dble(i-(nb+1)/2)*dx
        uc=0.5d0*(u(i,j)+u(i+1,j))
        vc=0.5d0*(v(i,j)+v(i,j+1))
        speed=sqrt(uc**2+vc**2)
        vmax=max(vmax,speed)
        r=sqrt(x*x+y*y)
        if (r>0.0d0) then
          vr=(x*uc+y*vc)/r
          vrmin=min(vrmin,vr)
        end if
      end do
    end do
  end subroutine diag_c

end program shallow_water

