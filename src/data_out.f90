! This file is a part of DUDI-heliocentric, the Fortran-90 implementation 
! of the two-body model for the dynamics of dust ejected from an atmosphereless
! body moving around the Sun
! Version 1.0.1
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 

! File: data_out.f90
! Description: Contains a subroutine used for output of the result.

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module data_out
    use const
    implicit none
        
         contains

                subroutine matrix_out(fname, image, nt1, nt2)
                    use const
                    implicit none
                    integer, intent(in) :: nt1, nt2
                    real, intent(in) :: image(nt1,nt2)
                    integer, parameter :: outchannel = 115
                    integer i 
                    character(*) fname
                    character(len = 20) outformat
                    
                    write(outformat,'(A1,I3,A1,I3,A11)') '(', nt1, 'x',nt2, '(ES12.4E2))'
                    
                    open(outchannel, file = fname, status = 'replace')
                    
                        do i = 1, nt1
                            write(outchannel,outformat) image(i,:)
                        enddo
                    
                    close(outchannel)
                
                end subroutine matrix_out

            
end module data_out
