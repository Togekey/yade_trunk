!deepak, 27/12/2014
!compile gfortran fibgen.f90 -o fibgen
!run ./fibgen
program fibgen 
! use random
implicit none 


! type defnitions 

! --- coords----


 type coords
    real::x(3) 
 end type coords 

 
! ---- fibre----- 
 
 type fibre 
  
    integer::myid 
    type(coords),dimension(:), allocatable::fseg 
    real :: rad !unused. 
    real :: ori 
    
    
  end type fibre 
  
 ! ----- create array of fibres----- 
 
 
 type(fibre),dimension(:),allocatable::nfib_array  
  
 
 ! --- some declarations 
  
  integer::seed=254135
  integer::i,j,k	! loop counters 
  integer::nofib 	! no. of fibres 
  integer::nseg 	! no. of nodes per fibres
  real::lfib 		! length of fibre 
  real,parameter::pi = 3.14159265
  real::box_len(3),el_len,dl 
  integer::ctr 
  
 !------------start----------- 
  
  
  nofib = 152			! always add 2 to to the total number required
  nseg  = 21 			! no. of nodes per fibre
  lfib = 0.02			! length of the fibre
  box_len(1) = 0.20             ! box dimensions
  box_len(2) = 0.20
  box_len(3) = 0.20
  ctr = int(nseg/2)
  dl = lfib/nseg 
  
  
  allocate(nfib_array(nofib))
  
  do i = 1,nofib 
    
    allocate(nfib_array(i)%fseg(nseg)) 
    
  end do 
  
  ! initialize the structures .. 
  
  
  do i = 1, nofib 
  
    nfib_array(i)%myid = i
    nfib_array(i)%ori = 0.0
    
    do j = 1, nseg 
      
      nfib_array(i)%fseg(j)%x(1)=0.0
      nfib_array(i)%fseg(j)%x(2)=0.0
      nfib_array(i)%fseg(j)%x(3)=0.0
    
    end do 
   
  end do 
  
  call srand(seed)
  
! set fibre coords  

  do i = 1, nofib 
    
    nfib_array(i)%fseg(ctr)%x(1) = abs(rand())*box_len(1)
    nfib_array(i)%fseg(ctr)%x(2) = abs(rand())*box_len(2)
    nfib_array(i)%fseg(ctr)%x(3) = abs(rand())*box_len(3)
    
  end do 
  
 
 ! generate orientation angles, 0 to pi 
 
 
 
  do i = 1, nofib 
    
    nfib_array(i)%ori  = (rand())*(pi)
 
  end do 
  
 ! generate segment coords 
 
 
 do i = 1, nofib 
  
 !change ori in x and y, z const 
 
  if (nfib_array(i)%myid.lt.int(nofib/3)) then 
  
    do j = 1,nseg 
      
      if(j.lt.ctr)then 
	
	el_len = (ctr-j)*dl
	nfib_array(i)%fseg(j)%x(1) = nfib_array(i)%fseg(ctr)%x(1) + (el_len*cos(nfib_array(i)%ori)) 
	nfib_array(i)%fseg(j)%x(2) = nfib_array(i)%fseg(ctr)%x(2) + (el_len*sin(nfib_array(i)%ori)) 
	nfib_array(i)%fseg(j)%x(3) = nfib_array(i)%fseg(ctr)%x(3)
  
      end if
      
     if(j.eq.ctr)then 
       nfib_array(i)%fseg(j)%x(1) = nfib_array(i)%fseg(ctr)%x(1)
       nfib_array(i)%fseg(j)%x(2) = nfib_array(i)%fseg(ctr)%x(2)
       nfib_array(i)%fseg(j)%x(3) = nfib_array(i)%fseg(ctr)%x(3)
     end if 
     
    if(j.gt.ctr)then 
      el_len = (j-ctr)*dl
      nfib_array(i)%fseg(j)%x(1) = nfib_array(i)%fseg(ctr)%x(1)+(-el_len*cos(nfib_array(i)%ori))
      nfib_array(i)%fseg(j)%x(2) = nfib_array(i)%fseg(ctr)%x(2)+(-el_len*sin(nfib_array(i)%ori)) 
      nfib_array(i)%fseg(j)%x(3) = nfib_array(i)%fseg(ctr)%x(3)
     end if 
 
    end do 
 
 
 end if 
 
  if (nfib_array(i)%myid.gt.int(nofib/3).and.nfib_array(i)%myid.lt.int(2*nofib/3)) then  ! change ori in x and z, y const
  
  
  
  do j = 1,nseg 
      
      if(j.lt.ctr)then 
	
	el_len = (ctr-j)*dl
	nfib_array(i)%fseg(j)%x(1) = nfib_array(i)%fseg(ctr)%x(1) + (el_len*cos(nfib_array(i)%ori)) 
	nfib_array(i)%fseg(j)%x(3) = nfib_array(i)%fseg(ctr)%x(3) + (el_len*sin(nfib_array(i)%ori))
	nfib_array(i)%fseg(j)%x(2) = nfib_array(i)%fseg(ctr)%x(2)
  
      end if
      
  if(j.eq.ctr)then 
       nfib_array(i)%fseg(j)%x(1) = nfib_array(i)%fseg(ctr)%x(1)
       nfib_array(i)%fseg(j)%x(2) = nfib_array(i)%fseg(ctr)%x(2)
       nfib_array(i)%fseg(j)%x(3) = nfib_array(i)%fseg(ctr)%x(3)
     end if 
        
    if(j.gt.ctr)then 
      el_len = (j-ctr)*dl
      nfib_array(i)%fseg(j)%x(1) = nfib_array(i)%fseg(ctr)%x(1)+(-el_len*cos(nfib_array(i)%ori))
      nfib_array(i)%fseg(j)%x(3) = nfib_array(i)%fseg(ctr)%x(3)+(-el_len*sin(nfib_array(i)%ori)) 
      nfib_array(i)%fseg(j)%x(2) = nfib_array(i)%fseg(ctr)%x(2)
     end if 
 
    end do 
  
  
  
  end if 
  

 
 ! change ori in y and z, x const
 
 
  if (nfib_array(i)%myid.gt.int(2*nofib/3)) then 
 
 do j = 1,nseg 
      
      if(j.lt.ctr)then 
	
	el_len = (ctr-j)*dl
	nfib_array(i)%fseg(j)%x(3) = nfib_array(i)%fseg(ctr)%x(3) + (el_len*cos(nfib_array(i)%ori)) 
	nfib_array(i)%fseg(j)%x(2) = nfib_array(i)%fseg(ctr)%x(2) + (el_len*sin(nfib_array(i)%ori)) 
	nfib_array(i)%fseg(j)%x(1) = nfib_array(i)%fseg(ctr)%x(1)
  
      end if
      
  if(j.eq.ctr)then 
       nfib_array(i)%fseg(j)%x(1) = nfib_array(i)%fseg(ctr)%x(1)
       nfib_array(i)%fseg(j)%x(2) = nfib_array(i)%fseg(ctr)%x(2)
       nfib_array(i)%fseg(j)%x(3) = nfib_array(i)%fseg(ctr)%x(3)
     end if 
        
    if(j.gt.ctr)then 
      el_len = (j-ctr)*dl
      nfib_array(i)%fseg(j)%x(3) = nfib_array(i)%fseg(ctr)%x(3)+(-el_len*cos(nfib_array(i)%ori))
      nfib_array(i)%fseg(j)%x(2) = nfib_array(i)%fseg(ctr)%x(2)+(-el_len*sin(nfib_array(i)%ori))
      nfib_array(i)%fseg(j)%x(1) = nfib_array(i)%fseg(ctr)%x(1)
     end if 
 
    end do 
 
  end if 
 
 
 
 end do 
 
 ! write fibre coords + connxn to file 
 
 open(unit=10,file="fibnodes.dat", action="write", status="replace")
 
 
 do i=1, nofib 
 
    if(i==int(nofib/3).or.i==int(2*nofib/3))then 
      continue
    else 
      write(10,*) nfib_array(i)%fseg
     end if
      
  
 end do 

 
 
 close(10)
 
 
end program fibgen
  
