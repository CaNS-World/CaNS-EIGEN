! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
program gen_grid
!
! this program generates the grid files to be read by the xdmf file
! the output files from CaNS 'grid_x.bin', 'grid_y.bin' and 'grid_z.bin'
! must be in this folder
!
! Pedro Costa (p.simoes.costa@gmail.com)
!
implicit none
include 'param.h90'
!
integer :: iunit
!
iunit = 99
call write_centered_grid(iunit,'grid_x.bin','x.bin',nx)
call write_centered_grid(iunit,'grid_y.bin','y.bin',ny)
call write_centered_grid(iunit,'grid_z.bin','z.bin',nz)
contains
  subroutine write_centered_grid(iunit,gridfile,outfile,n)
    implicit none
    integer         , intent(in) :: iunit,n
    character(len=*), intent(in) :: gridfile,outfile
    real(8), allocatable, dimension(:) :: dummy,rc
    !
    allocate(dummy(n),rc(n))
    open(iunit,file=gridfile,action='read',form='unformatted',access='stream',status='old')
    read(iunit) dummy(1:n),dummy(1:n),rc(1:n),dummy(1:n)
    close(iunit)
    open(iunit,file=outfile,action='write',form='unformatted',access='stream',status='replace')
    write(iunit) rc(1:n)
    close(iunit)
    deallocate(dummy,rc)
  end subroutine write_centered_grid
end program gen_grid
