module m
        use iso_fortran_env

        contains
        subroutine cop_carnot()
        	implicit none
      		
      		real(kind=8) :: t, q_c, q_h, min, max
      		character :: category
      		integer :: i, n, io
      		real(kind=8), dimension(:,:), allocatable :: m
        	
        	
       		open(unit=100, file='Qc.dat')
        	open(unit=200, file='Qh.dat')
        	open(unit=300, file='result.dat')
        
        	print*, "number of lines"
        	n=0
        	do 
        		read(100,*,iostat=io) t, category, q_c!min, max, q_c
        		if (io==iostat_end) exit
    			n=n+1
        	end do
        	
        	print*, n
        	allocate (m(n,4))
        	
        	rewind(100)
        	
        	print*, "reading and writing qc"
        	
        	do i=1,n
        		read(100,*,iostat=io) t, category, q_c!min, max, q_c
        		if (io==iostat_end) exit
        		m(i,1)= t 
        		m(i,2)=q_c
        	end do
        	
        	print*, "reading and writing qh"
        	
        	do i=1,n
        		read(200,*,iostat=io) t, category, q_h!min, max, q_h
        		if (io==iostat_end) exit
        		m(i,3)=q_h
        	end do
        	
        	do i=1,n
        		m(i,4)= abs(m(i, 2))/(abs(m(i, 3))-abs(m(i, 2)))
        		write(300,*) m(i, :)
        	end do
        	
        	close(100)
        	close(200)
        	close(300)
        
        end subroutine
        subroutine work()
        	implicit none
        	
        	integer :: io, i, j, k, len_modified
        	real, parameter :: T=0.005, p_m=100000
        	real(kind=8), dimension(14) :: t_2, dt, sum_p, sum_u
        	real(kind=8), dimension(14, 40) :: mat_p, mat_pt
        	real(kind=8), dimension(14, 126) :: mat_u, mat_ux, mat_ut
        	character(len=4000) :: line, modified_string, line_p
        	character(len=20) :: line_t, t_1
        	!real(kind=8) :: t_1
        	real(kind=8) :: sum_w
		
		
        	
        	open(unit=400, file='p_end')
        	open(unit=500, file='u_end')
        	open(unit=600, file='p.txt')
        	open(unit=700, file='velocity.txt')
        	
        	print*, "time "
        	print*, ""
        	
        	do
        		read(400, '(F20.15)', iostat=io) t_2
        		if (io==iostat_end) exit
    			print*,  t_2
        	end do
        	
        	call delta(t_2, dt)
        	
        	print*, ""
        	print*, "time steps"
        	print*, ""
        	
        	do i=1, 14
        		print*, dt(i)
        	end do
        	
        	rewind (400)
        	
        	
        	!print*, ""
        	!print*, "p x dt"
        	!print*, ""
        		
        	!do i=1, 14
        	!	read(400, *, iostat=io) mat_p(i,:)
        	!	if (io==iostat_end) exit
        		!mat_pt(i,:)=(mat_p(i,:)-100000)*dt(i)
        		!print*, mat_p(i,:)-100000, " X ", dt(i), " = ", mat_pt(i,:)
        		!print*, mat_p(i,:)
        	!end do
        	!'(25x, 40F13.7)'
        	
        	
        	print*, ""
        	print*, "p matrix as characters"
        	print*, ""
        	
        	do i=1, 14
        		read(400,'(a20,A2000)',iostat=io)line_t, line_p
        		if (io==iostat_end) exit
        		!print*, line_p
			write(600,*) trim(line_p(1:2000))
        	end do
        	
        	rewind (600)
        	
        	print*, ""
        	print*, "p x dt"
        	print*, ""
        	
        	do i=1, 14
        		read(600, *, iostat=io) mat_p(i, :)
        		if (io==iostat_end) exit
        		!print*, mat_p(i,:)
        		mat_pt(i,:)=(mat_p(i,:)-100000)*dt(i)
        		print*, mat_p(i,:)-100000, " X ", dt(i), " = ", mat_pt(i,:)
        	end do
        	
        	print*, ""
        	print*, "sum p by number of samples"
        	print*, ""
        	
        	do i=1, 14
        		sum_p(i)=sum(mat_p(i,:)-100000)/40
        		print*, sum_p(i)
        	end do
   	
        	print*, ""
        	print*, "velocity as a character"
        	print*, ""
        	
        	do i=1, 14
        		read(500,'(a20,A4000)',iostat=io)t_1, line
        		if (io==iostat_end) exit
        		!print*, line
        	
        		len_modified=0	
        		
        		do k=1,4000
        			if (line(k:k) /= '(' .and. line(k:k) /= ')') then
        				len_modified = len_modified + 1
            				modified_string(len_modified:len_modified) = line(k:k)
        			end if
			end do
			!print*, trim(modified_string(1:len_modified))
			write(700,*) trim(modified_string(1:len_modified))
        	end do
        	
        	rewind (700)
        	
        	print*, ""
        	print*, "velocity X t"
        	print*, ""
        	
        	do i=1, 14
        		read(700, *, iostat=io) (mat_u(i, j), j=1, 126)
        		if (io==iostat_end) exit
        		!print*, mat_u(i,:)
        		mat_ut(i,:)=mat_u(i,:)*dt(i)
        		print*, mat_u(i,:), " X ", dt(i), " = ", mat_ut(i,:)
        	end do
        	
        	sum_u=0

        	do i=1, 14
        		do j=1,126,3
        			sum_u(i)=sum_u(i)+mat_u(i,j)
        		end do
        	end do
        	
        	print*, ""
        	print*, "sum ux"
        	print*, ""
        	
        	do i=1, 14
        		print*, sum_u(i)/40
        	end do
        	
        	print*, ""
        	print*, "work"
        	print*, ""
        	sum_w=0
        	
        	do i=1, 14
        		sum_w=sum_p(i)*sum_u(i)*0.000025+sum_w
        		print*, sum_w
        	end do
        	
        	print*, "work = ", sum_w
      	
 	
    		close (400)
        	close (500)
        	close (600)
        	close (700)
        end subroutine
	
	subroutine delta(t, dt)
		implicit none
		
		real(kind=8), dimension(:) :: t, dt
		integer :: i, k
		
		k = size(t) - 1
		
		do i=1, k-1
			if (i==0) then 
				dt(i) = 0
			else
				dt(i)=t(i)-t(i-1)
			end if
		end do
	
	end subroutine
	
end module
program hello
	use m
	
	print*, "start"
	
	!call cop_carnot
	call work
	
	print*, "end"
	
end program
