program unit2
    implicit none

    !------------------------------------------------------------
    ! Kich thuoc luoi
    !------------------------------------------------------------
    integer, parameter :: nx = 1668
    ! nx: so diem luoi khong gian

    integer, parameter :: nhour = 24
    ! nhour: so gio luu ket qua

    !------------------------------------------------------------
    ! Tham so mo hinh
    !------------------------------------------------------------
    real(kind=8), parameter :: dx = 300.0d0
    ! dx: buoc luoi khong gian (m)

    real(kind=8), parameter :: dt = 10.0d0
    ! dt: buoc thoi gian (s)

    real(kind=8), parameter :: Cf = 5.0d0
    ! Cf: van toc binh luu (m/s)

    real(kind=8), parameter :: C = 5.0d0
    ! C: van toc dung trong dieu kien bien Orlanski (m/s)

    real(kind=8), parameter :: t_end = 24.0d0 * 3600.0d0
    ! t_end: thoi gian tinh toan 24 gio, doi sang giay

    !------------------------------------------------------------
    ! Khai bao bien
    !------------------------------------------------------------
    real(kind=8) :: x(nx)
    ! x: toa do khong gian (m)

    real(kind=8) :: U(nx), Unew(nx)
    ! U: nghiem U tai thoi diem hien tai
    ! Unew: nghiem U tai thoi diem moi

    real(kind=8) :: Uhour(nx,0:nhour)
    ! Uhour: luu U sau moi 1 gio
    ! Uhour(:,0) la nghiem ban dau

    real(kind=8) :: time
    ! time: thoi gian hien tai (s)

    integer :: i, istep, nsteps, hour, noutput
    ! i: chi so diem luoi
    ! istep: chi so buoc thoi gian
    ! nsteps: tong so buoc thoi gian
    ! hour: so gio hien tai
    ! noutput: so buoc tinh trong 1 gio

    ! Tong so buoc thoi gian
    nsteps = int(t_end / dt)

    ! So buoc tinh trong 1 gio
    noutput = int(3600.0d0 / dt)

    !------------------------------------------------------------
    ! Khoi tao song tam giac
    !------------------------------------------------------------
    call ini(U, x, dx, nx)

    ! Luu nghiem tai t = 0
    Uhour(:,0) = U

    !------------------------------------------------------------
    ! Vong lap tinh phuong trinh binh luu
    !------------------------------------------------------------
    do istep = 1, nsteps

        ! Thoi gian hien tai
        time = real(istep, kind=8) * dt

        ! Tinh U moi bang phuong phap RK4
        call forecast(U, Unew, nx, dx, Cf, dt)

        ! Ap dung dieu kien bien Orlanski
        call boundary(Unew, U, nx, dx, dt, C)

        ! Cap nhat nghiem
        U = Unew

        ! Luu ket qua sau moi 1 gio
        if (mod(istep, noutput) == 0) then

            ! Doi thoi gian tu giay sang gio
            hour = int(time / 3600.0d0)

            if (hour <= nhour) then
                Uhour(:,hour) = U
            end if

            ! In thoi gian ra man hinh
            write(*,'(A,I2,A)') 'Time = ', hour, ' h'

        end if

    end do

    !------------------------------------------------------------
    ! Ghi ket qua ra file TXT
    !------------------------------------------------------------
    open(unit=10, file='result_unit2.txt', &
         status='replace', action='write')

    ! Ghi ten cot x
    write(10,'(A)',advance='no') 'x(km)'

    ! Ghi ten cac cot U theo tung gio
    do hour = 0, nhour
        write(10,'(A,I2,A)',advance='no') &
            achar(9)//'U(', hour, 'h)'
    end do

    write(10,*)

    ! Ghi du lieu
    do i = 1, nx

        ! Doi x tu m sang km
        write(10,'(F10.3)',advance='no') x(i)/1000.0d0

        ! Ghi U tu 0 gio den 24 gio
        do hour = 0, nhour
            write(10,'(A,ES12.4)',advance='no') &
                achar(9), Uhour(i,hour)
        end do

        write(10,*)

    end do

    close(10)

    print *
    print *, 'Done! result_unit2.txt generated.'


contains


    !============================================================
    ! HAM KHOI TAO SONG TAM GIAC
    !============================================================
    !
    ! Tao nghiem ban dau U(x,0) dang hinh tam giac.
    !
    ! xA: chan trai cua song
    ! xC: dinh song
    ! xB: chan phai cua song
    ! Umax: gia tri lon nhat cua U
    !
    subroutine ini(U, x, dx, nx)

        implicit none

        integer, intent(in) :: nx
        ! nx: so diem luoi

        real(kind=8), intent(out) :: U(nx)
        ! U: gia tri song ban dau

        real(kind=8), intent(out) :: x(nx)
        ! x: toa do cac diem luoi

        real(kind=8), intent(in) :: dx
        ! dx: buoc luoi khong gian

        integer :: i
        ! i: chi so diem luoi

        real(kind=8) :: xA, xC, xB, Umax
        ! xA: vi tri chan trai
        ! xC: vi tri dinh
        ! xB: vi tri chan phai
        ! Umax: gia tri cuc dai cua song

        Umax = 15.0d0
        xA = 6000.0d0
        xC = 9000.0d0
        xB = 12000.0d0

        do i = 1, nx

            ! Tao toa do x
            x(i) = real(i-1, kind=8) * dx

            ! Ben trai song: U = 0
            if (x(i) < xA) then

                U(i) = 0.0d0

            ! Canh len cua tam giac
            else if (x(i) <= xC) then

                U(i) = Umax * (x(i)-xA) / (xC-xA)

            ! Canh xuong cua tam giac
            else if (x(i) <= xB) then

                U(i) = Umax * (xB-x(i)) / (xB-xC)

            ! Ben phai song: U = 0
            else

                U(i) = 0.0d0

            end if

        end do

    end subroutine ini


    !============================================================
    ! HAM TINH DAO HAM KHONG GIAN
    !============================================================
    !
    ! Phuong trinh:
    !
    ! dU/dt = -Cf * dU/dx
    !
    ! Dao ham dU/dx duoc tinh bang sai phan trung tam:
    !
    ! dU/dx = [U(i+1) - U(i-1)] / (2*dx)
    !
    subroutine rhs(U, dUdx, nx, dx)

        implicit none

        integer, intent(in) :: nx
        ! nx: so diem luoi

        real(kind=8), intent(in) :: U(nx)
        ! U: nghiem hien tai

        real(kind=8), intent(out) :: dUdx(nx)
        ! dUdx: dao ham khong gian cua U

        real(kind=8), intent(in) :: dx
        ! dx: buoc luoi khong gian

        integer :: i
        ! i: chi so diem luoi

        ! Khoi tao dao ham bang 0
        dUdx = 0.0d0

        ! Sai phan trung tam trong mien tinh
        do i = 2, nx-1

            dUdx(i) = (U(i+1) - U(i-1)) / (2.0d0*dx)

        end do

        ! Hai bien duoc xu ly trong ham boundary
        dUdx(1)  = 0.0d0
        dUdx(nx) = 0.0d0

    end subroutine rhs


    !============================================================
    ! HAM TINH RK4
    !============================================================
    !
    ! Ap dung phuong phap Runge-Kutta bac 4 (RK4)
    ! cho phuong trinh binh luu:
    !
    ! dU/dt = -Cf * dU/dx
    !
    subroutine forecast(U, Unew, nx, dx, Cf, dt)

        implicit none

        integer, intent(in) :: nx
        ! nx: so diem luoi

        real(kind=8), intent(in) :: U(nx)
        ! U: nghiem tai thoi diem hien tai

        real(kind=8), intent(out) :: Unew(nx)
        ! Unew: nghiem sau mot buoc dt

        real(kind=8), intent(in) :: dx
        ! dx: buoc luoi khong gian

        real(kind=8), intent(in) :: Cf
        ! Cf: van toc binh luu

        real(kind=8), intent(in) :: dt
        ! dt: buoc thoi gian

        real(kind=8) :: k1(nx), k2(nx), k3(nx), k4(nx)
        ! k1, k2, k3, k4: cac buoc trung gian cua RK4

        real(kind=8) :: Utmp(nx)
        ! Utmp: nghiem trung gian trong RK4

        !-------------------------
        ! Buoc 1
        !-------------------------
        call rhs(U, k1, nx, dx)
        k1 = -Cf * k1

        !-------------------------
        ! Buoc 2
        !-------------------------
        Utmp = U + 0.5d0 * dt * k1

        call rhs(Utmp, k2, nx, dx)
        k2 = -Cf * k2

        !-------------------------
        ! Buoc 3
        !-------------------------
        Utmp = U + 0.5d0 * dt * k2

        call rhs(Utmp, k3, nx, dx)
        k3 = -Cf * k3

        !-------------------------
        ! Buoc 4
        !-------------------------
        Utmp = U + dt * k3

        call rhs(Utmp, k4, nx, dx)
        k4 = -Cf * k4

        !-------------------------
        ! Cong thuc RK4
        !-------------------------
        !
        ! U(n+1) = U(n)
        !        + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        !
        Unew = U + (dt/6.0d0) * &
               (k1 + 2.0d0*k2 + 2.0d0*k3 + k4)

    end subroutine forecast


    !============================================================
    ! HAM DIEU KIEN BIEN ORLANSKI
    !============================================================
    !
    ! Xu ly gia tri U tai hai bien cua mien tinh.
    !
    subroutine boundary(Unew, Uold, nx, dx, dt, C)

        implicit none

        integer, intent(in) :: nx
        ! nx: so diem luoi

        real(kind=8), intent(inout) :: Unew(nx)
        ! Unew: nghiem moi, duoc cap nhat tai bien

        real(kind=8), intent(in) :: Uold(nx)
        ! Uold: nghiem tai buoc thoi gian truoc

        real(kind=8), intent(in) :: dx
        ! dx: buoc luoi khong gian

        real(kind=8), intent(in) :: dt
        ! dt: buoc thoi gian

        real(kind=8), intent(in) :: C
        ! C: van toc dung trong dieu kien Orlanski

        real(kind=8) :: alpha
        ! alpha: he so C*dt/dx

        ! He so trong cong thuc Orlanski
        alpha = C * dt / dx

        ! Bien trai: U = 0
        Unew(1) = 0.0d0

        ! Bien phai: dieu kien bien Orlanski
        !
        ! U(nx,t+dt) = U(nx,t)
        !               - alpha*[U(nx,t)-U(nx-1,t)]
        !
        Unew(nx) = Uold(nx) - &
                   alpha * (Uold(nx) - Uold(nx-1))

    end subroutine boundary


end program unit2