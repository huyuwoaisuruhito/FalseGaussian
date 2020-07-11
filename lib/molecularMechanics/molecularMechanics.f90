! module mm_f
! contains

subroutine cross_product (v1, v2, ans)
    implicit none
    REAL(8), intent(in), dimension(0:2):: v1
    REAL(8), intent(in), dimension(0:2):: v2
    REAL(8), intent(out), dimension(0:2):: ans

    ans(0) = v1(1) * v2(2) - v1(2) * v2(1)
    ans(1) = v1(2) * v2(0) - v1(0) * v2(2)
    ans(2) = v1(0) * v2(1) - v1(1) * v2(0)
end subroutine cross_product

subroutine dot_product (v1, v2, ans)
    implicit none
    REAL(8), intent(in), dimension(0:2):: v1
    REAL(8), intent(in), dimension(0:2):: v2
    REAL(8), intent(out):: ans

    ans = v1(0) * v2(0) + v1(1) * v2(1) + v1(2) * v2(2);
end subroutine dot_product

subroutine potential_dihedral (a1, a2, a3, a4, dihedraldata, len, potential)
    REAL(8), intent(in), dimension(0:2):: a1
    REAL(8), intent(in), dimension(0:2):: a2
    REAL(8), intent(in), dimension(0:2):: a3
    REAL(8), intent(in), dimension(0:2):: a4
    REAL(8), intent(in), dimension(0:3*len-1):: dihedraldata
    integer, intent(in):: len
    REAL(8), intent(out):: potential

    REAL(8), dimension(0:2):: d1, d2, d3
    REAL(8), dimension(0:2):: da1, da2
    REAL(8):: cosa, phi, dotp, lda1, lda2

    do i = 0, 2
        d1(i) = a2(1) - a1(i)
        d2(i) = a3(1) - a2(i)
        d3(i) = a4(1) - a3(i)
    end do

    call cross_product(d1, d2, da1)
    call cross_product(d2, d3, da2)
    call dot_product(da1, da2, dotp)
    call dot_product(da1,da1, lda1)
    call dot_product(da2,da2, lda2)
    cosa = dotp / sqrt(lda1 * lda2)
    if (-1 < cosa .AND. cosa < 1) then
        phi = acos(cosa)
    else
        if (cosa > 0) then
            phi = 0
        else
            phi = dacos(-1.D0)
        end if
    end if

    potential = 0
    do i = 0, len-1
        potential = potential + dihedraldata(3*i+2) * cos(dihedraldata(3*i+0) * phi - dihedraldata(3*i+1))
    end do
end subroutine potential_dihedral

function fdihedral_c(a1, a2, a3, a4, dihedraldata, len, force)
    !DEC$ ATTRIBUTES DLLEXPORT::fdihedral_c
    REAL(8), intent(in), dimension(0:2):: a1
    REAL(8), intent(in), dimension(0:2):: a2
    REAL(8), intent(in), dimension(0:2):: a3
    REAL(8), intent(in), dimension(0:2):: a4
    REAL(8), intent(in), dimension(0:3*len-1):: dihedraldata
    integer, intent(in):: len
    REAL(8), intent(out), dimension(0:2):: force

    REAL(8):: temp1, temp2
    REAL(8), dimension(0:2):: da2
    logical :: fdihedral_c

    call potential_dihedral(a1, a2, a3, a4, dihedraldata, len, temp1)
    do i=0, 2
        da2(i) = a2(i) + 0.00001
        call potential_dihedral(a1, da2, a3, a4, dihedraldata, len, temp2)
        force(i) = ( temp1-temp2 ) / 0.00001
        da2(i) = a2(i)
    end do 

    fdihedral_c = .True.
end function fdihedral_c

function fdihedral_s(a1, a2, a3, a4, dihedraldata, len, force)
    !DEC$ ATTRIBUTES DLLEXPORT::fdihedral_s
    REAL(8), intent(in), dimension(0:2):: a1
    REAL(8), intent(in), dimension(0:2):: a2
    REAL(8), intent(in), dimension(0:2):: a3
    REAL(8), intent(in), dimension(0:2):: a4
    REAL(8), intent(in), dimension(0:3*len-1):: dihedraldata
    integer, intent(in):: len
    REAL(8), intent(out), dimension(0:2):: force

    REAL(8):: temp1, temp2
    REAL(8), dimension(0:2):: da1
    logical :: fdihedral_s

    do i=0, 2
        da1(i) = a1(i)
    end do 

    call potential_dihedral(a1, a2, a3, a4, dihedraldata, len, temp1)
    if (temp1 > 100000 .OR. temp1 < -100000) then
        write(*,*) 'Error:', temp1
    end if
    do i=0, 2
        da1(i) = a1(i) + 0.00001
        call potential_dihedral(da1, a2, a3, a4, dihedraldata, len, temp2)
        force(i) = ( temp1-temp2 ) / 0.00001
        da1(i) = a1(i)
    end do

    fdihedral_s = .True.
end function fdihedral_s

! end module mm_f


! program test
!     use mm_f

!     integer:: len = 1
!     REAL(8), dimension(0:2):: a1 = [3.1011916589473683, -2.9756210204210527, -1.5546421468421052]
!     REAL(8), dimension(0:2):: a2 = [2.101631658947368, -2.945955020421053, -1.5546421468421052]
!     REAL(8), dimension(0:2):: a3 = [1.4146776589473684, -1.6701430204210528, -1.5546421468421052]
!     REAL(8), dimension(0:2):: a4 = [ 0.78984466, -1.59334002, -2.44446515]
!     REAL(8), dimension(0:2):: dihedraldata = [3., 0., 0.65084]
!     REAL(8), dimension(0:2):: ans

!     write(*, *) fdihedral_s(a1, a2, a3, a4, dihedraldata, len, ans)

!     write(*, *) ans
!     READ *, pause

! end program test