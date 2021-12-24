!======================================================================!
!
       program randconf
!
       use datatypes
       use geometry
!
       implicit none
!
       include 'parameters.h'
       include 'info.h'
!
       character(len=leninp)                  ::  inp     !  Input file name
       character(len=lenout)                  ::  outp    !  Output file name
       character(len=lenarg)                  ::  molaux  !  Auxliary string
       real(kind=8),dimension(ndim)           ::  box     !  Simulation box
       real(kind=8)                           ::  rho     !  Density
       real(kind=8)                           ::  r0      !  Minimum separation
       real(kind=8)                           ::  rmin    !  Corrected miminum separation
       integer,dimension(:),allocatable       ::  nmol    !  Number of molecules
       integer,dimension(ndim)                ::  nbox    !  Unit cells per direction
       integer                                ::  ntrial  !  Maximum number of trials
       integer                                ::  nfile   !  Number of input files
       logical                                ::  debug   !  Debug mode
       character(len=lenarg)                  ::  straux  !  Auxiliary string
       integer                                ::  io      !  Input/Output status
       integer                                ::  i,j,k   !  Indexes
!
! Reading command line options
!
       call command_line(inp,outp,nfile,molaux,box,nbox,               &
                         rho,r0,ntrial,debug)
!
       allocate(mol(nfile),nmol(nfile))
!
!  General settings and fatal errors check
!
       do i = 1, nfile
         molaux = adjustl(molaux)
         j      = index(molaux,' ')
         straux = molaux(:j-1)
         read(straux,*,iostat=io) nmol(i)
         molaux = molaux(j+1:)
       end do
!
       do i = 1, nfile
         inp          = adjustl(inp)
         j            = index(inp,' ')
         mol(i)%fname = inp(:j-1)
         inp          = inp(j+1:)
       end do
!
!~        i = index(outp,'.')
!~        if ( i .gt. 0 ) outp = outp(:i-1)
!
!
! Processing Gromacs input files
!
       call read_gro(nfile)
!
       do i = 1, nfile
         call translate(mol(i)%nat,mol(i)%coord,                       &
                        -1.0d0,CofM_vector(mol(i)%nat,mol(i)%coord,    &
                                           mol(i)%mass),               &
                        mol(i)%coord)
       end do
!
! Memory allocation
!
       sys%nat = 0
       do i = 1, nfile
         sys%nat = sys%nat + mol(i)%nat*nmol(i)
       end do
!
       allocate(sys%renum(sys%nat)   ,  &
                sys%rename(sys%nat)  ,  &
                sys%atname(sys%nat)  ,  &
                sys%atnum( sys%nat)  ,  &
                sys%coord(3,sys%nat) )
!
! Generating the random initial configuration
!
       call randcoord(nfile,nmol,box,r0,ntrial,outp,debug)
!
!~        do i = 1, nfile
!~            write(*,*) trim(mol(i)%title)
!~            write(*,*) mol(i)%nat
!~          do k = 1, mol(i)%nat
!~            write(*,'(I5,2A5,I5,3F8.3)') mol(i)%renum(k),   &
!~                                         mol(i)%rename(k),  &
!~                                         mol(i)%atname(k),  &
!~                                         mol(i)%atnum(k),   &
!~                                         mol(i)%coord(:,k)
!~          end do
!~            write(*,'(3F10.5)') mol(i)%latvec
!~        end do
!
! Deallocate memory
!
       deallocate (mol,nmol)
       deallocate(sys%renum,sys%rename,sys%atname,sys%atnum,sys%coord)
!
       end program randconf
!
!======================================================================!
!
       subroutine command_line(inp,outp,nfile,molaux,box,nbox,rho,     &
                               r0,ntrial,debug)
!
       use flags
!
       implicit none
!
       include 'parameters.h'
       include 'info.h'
!
! Input/output variables
!
       character(len=leninp),intent(out)         ::  inp     !  Input file name
       character(len=lenout),intent(out)         ::  outp    !  Output file name
       character(len=lenarg),intent(out)         ::  molaux  !  Auxliary string
       real(kind=8),dimension(ndim),intent(out)  ::  box     !  Simulation box
       real(kind=8),intent(out)                  ::  rho     !  Density
       real(kind=8),intent(out)                  ::  r0      !  Minimum separation
       integer,dimension(ndim),intent(out)       ::  nbox    !  Unit cells per direction
       integer,intent(out)                       ::  nfile   !  Number of input files
       integer,intent(out)                       ::  ntrial  !  Maximum number of trials
       logical,intent(out)                       ::  debug   !  Debug mode
!
! Local variables
!
       character(len=8),parameter                ::  version = 'v0.0.0'
       character(len=lencmd)                     ::  cmd     !  Command executed
       character(len=lenarg)                     ::  code    !  Executable name
       character(len=lenarg)                     ::  arg     !  Argument read
       character(len=lenarg)                     ::  next    !  Next argument to be read
       integer                                   ::  naux    !  Auxliary variable
       integer                                   ::  io      !  Status
       integer                                   ::  i       !  Index
!
! Setting defaults
!
       inp    = 'conf.gro'
       outp   = 'sys.gro'
       nfile  = 1
       molaux = '512'
       box    = (/ 5.0d0, 5.0d0, 5.0d0 /)
       ntrial = 100
       r0     = 0.5d0
       debug  = .FALSE.
!
! Reading command line
!
       call get_command_argument(0,code)
       call get_command(cmd)
! Checking if any argument has been introduced
       if ( command_argument_count().eq.0) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')    'ERROR:  No argument introduced on'//    &
                              ' command-line'
         write(*,*)
         write(*,'(3X,A)')    'Please, to request help execute'
             write(*,*)
         write(*,'(4X,2(A))')  code(1:len_trim(code)), ' -h'
         write(*,'(2X,68("="))')
         write(*,*)
         call exit(0)
       end if
! Reading command line options
       i = 1
       do
         call get_command_argument(i,arg)
         if ( len_trim(arg) == 0 ) exit
         i = i+1
         select case ( arg )
           case ('-f','-file','-files','--file','--files')
             flginp = .TRUE.
             call read_string(i,leninp,inp,naux,arg,cmd)
             if ( flgmol ) then
               if ( naux .ne. nfile ) then
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A)')      'ERROR:  Rank mismatch in '//  &
                                        'arguments'
                 write(*,*)
                 write(*,'(3X,A,X,I4)') 'Number of values specifie'//  &
                                        'd after -nmol   :',nfile
                 write(*,'(3X,A,X,I4)') 'Number of input files spe'//  &
                                        'cified after -f :',naux
                 write(*,'(2X,68("="))')
                 call exit(0)
               end if
             else
               nfile = naux
             end if
           case ('-o','-out','-outp','-output','--output')
             call get_command_argument(i,outp,status=io)
             call check_arg(outp,io,arg,cmd)
             i = i + 1
           case ('-b','-box','--box','--simulation-box')
             flgbox = .TRUE.
             call read_realvec(i,ndim,box)
           case ('-n','-nmol','--nmol')
             flgmol = .TRUE.
             call read_string(i,lenarg,molaux,naux,arg,cmd)
             if ( flginp ) then
               if ( naux .ne. nfile ) then
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A)')      'ERROR:  Rank mismatch in '//  &
                                        'arguments'
                 write(*,*)
                 write(*,'(3X,A,X,I4)') 'Number of input files spe'//  &
                                        'cified after -f :',nfile
                 write(*,'(3X,A,X,I4)') 'Number of values specifie'//  &
                                        'd after -nmol   :',naux
                 write(*,'(2X,68("="))')
                 call exit(0)
               end if
             else
               nfile = naux
             end if
           case ('-r','-r0','--r0','--minimum-distance')
             flgdist = .TRUE.
             call get_command_argument(i,next,status=io)
             read(next,*) r0
             i = i + 1
           case ('-ntrial','-ntrials','--ntrial','--ntrials')
             flgntrl = .TRUE.
             call get_command_argument(i,next,status=io)
             read(next,*) ntrial
             i = i + 1
           case ('-nbox','--nbox')
             flgnbox = .TRUE.
             call read_intvec(i,ndim,nbox)
           case ('-rho','--rho','--density')
             flgrho = .TRUE.
             call get_command_argument(i,next,status=io)
             read(next,*) rho
             i = i + 1
           case ('-d','--debug','--verbose')
             debug = .TRUE.
             i = i + 1
           case ('-v','--version')
             write(*,*)
             write(*,'(2X,2A)') 'Current version ',version
             write(*,*)
             call exit(0)
           case ('-h','-help','--help')
             call print_help()
             call exit(0)
           case default
             write(*,*)
             write(*,'(2X,68("="))')
             write(*,'(3X,A)')    'ERROR:  Unknown statements from'//  &
                                  ' command line'
             write(*,*)
             write(*,'(4X,A)')     trim(cmd)
             write(*,*)
             write(*,'(3X,2(A))') 'Unrecognised command-line option'// &
                                  '  :  ', arg
             write(*,'(3X,A)')    'Please, to request help execute'
             write(*,*)
             write(*,'(4X,2(A))')  code(1:len_trim(code)), ' -h'
             write(*,'(2X,68("="))')
             write(*,*)
             call exit(0)
         end select
       end do
!
       return
       end subroutine command_line
!
!======================================================================!
!
       subroutine print_help()
!
       implicit none
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'Command-line options:'
       write(*,*)
       write(*,'(5X,A)') '-h,--help      Print usage information'//    &
                         ' and exit'
       write(*,'(5X,A)') '-f,--file      Input file name(s)'
       write(*,'(5X,A)') '-o,--output    Output file name'
       write(*,'(5X,A)') '-n,--nmol      Number of molecules'
       write(*,'(5X,A)') '-b,--box       Simulation box dimensions'
       write(*,'(5X,A)') '-r,--r0        Minimum separation'
       write(*,'(5X,A)') '-ntrial        Maximum number of trials '//  &
                         'before decreasing the minimum separation' 
       write(*,'(5X,A)') '-rho           Density (g/cm**3)'
       write(*,'(5X,A)') '-nbox          Unit cells per direction '//  &
                         '(nx,ny,nz)'
       write(*,'(5X,A)') '-d,--debug     Debug mode'
       write(*,'(2X,68("="))')
       write(*,*)
!
       return
       end subroutine print_help
!
!======================================================================!
!
       subroutine check_arg(opt,io,arg,cmd)
!
       implicit none
!
       character(len=*),intent(in)  ::  opt
       character(len=*),intent(in)  ::  arg
       character(len=*),intent(in)  ::  cmd
       integer,intent(in)           ::  io
!
       if ( (io .ne. 0) .or. (opt(1:1) .eq. '-') ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')    'ERROR:  No argument introduced'//      &
                              ' for command-line'//   &
                              ' option'
         write(*,*)
         write(*,'(4X,A)')    trim(cmd)
         write(*,*)
         write(*,'(3X,2(A))') 'Argument missing for command'//        &
                              '-line option  :  ', arg
         write(*,'(2X,68("="))')
         write(*,*)
         call exit(0)
       end if
!
       return
       end subroutine check_arg
!
!======================================================================!
!
       subroutine read_string(i,lenstr,inp,nfile,arg,cmd)
!
       implicit none
!
       include 'info.h'
!
! Input/output variables
       integer,intent(inout)              ::  i       !  Argument index
       integer,intent(in)                 ::  lenstr  !  Input length
       character(len=lenstr),intent(out)  ::  inp     !  Input file name
       integer,intent(out)                ::  nfile   !  Number of input files
       character(len=lenarg),intent(in)   ::  arg     !  Argument read
       character(len=lencmd),intent(in)   ::  cmd     !  Command executed
! Local variables
       character(len=lenarg)              ::  next    !  Next argument to be read
       integer                            ::  io      !  Status
!
       call get_command_argument(i,next,status=io)
       call check_arg(next,io,arg,cmd)
       nfile = 1
       inp   = next
       do
         call get_command_argument(i+nfile,next,status=io)
         if ( (io .ne. 0) .or. (next(1:1) .eq. '-') ) exit
         nfile = nfile + 1
         inp   = trim(inp)//' '//trim(next)
       end do
       i = i + nfile
!
       return
       end subroutine read_string
!
!======================================================================!
!
       subroutine read_realvec(i,n,box)
!
       implicit none
!
       include 'info.h'
!
! Input/output variables
       integer,intent(inout)                  ::  i       !  Argument index
       integer,intent(in)                     ::  n       !  Vector dimension
       real(kind=8),dimension(n),intent(out)  ::  box     !  Double precision vector
! Local variables
       character(len=lenarg)                  ::  next    !  Next argument to be read
       integer                                ::  io      !  Status
       integer                                ::  k       !  Index
!
       do k = 1, n
         call get_command_argument(i,next,status=io)
         read(next,*) box(k)
         i = i + 1
       end do
!
       return
       end subroutine read_realvec
!
!======================================================================!
!
       subroutine read_intvec(i,n,nbox)
!
       implicit none
!
       include 'info.h'
!
! Input/output variables
       integer,intent(inout)             ::  i       !  Argument index
       integer,intent(in)                ::  n       !  Vector dimension
       integer,dimension(n),intent(out)  ::  nbox    !  Integer vector
! Local variables
       character(len=lenarg)             ::  next    !  Next argument to be read
       integer                           ::  io      !  Status
       integer                           ::  k       !  Index
!
       do k = 1, n
         call get_command_argument(i,next,status=io)
         read(next,*) nbox(k)
         i = i + 1
       end do
!
       return
       end subroutine read_intvec
!
!======================================================================!
!
       subroutine read_gro(n)
!
       use datatypes
!
       implicit none
!
       include 'info.h'
!
! Input/output variables
       integer,intent(in)                       ::  n       !  Number of input files
! Local variables
       character(len=lenarg)                    ::  straux  !  Auxiliary string
       character(len=5)                         ::  aux     !
       integer                                  ::  io      !  Input/Output status
       integer                                  ::  i,j,k   !  Indexes
!
       do i = 1, n
         open(unit=uniinp,file=trim(mol(i)%fname),action='read',       &
              status='old',iostat=io)
!
         if ( io .ne. 0 ) then
           write(*,'(2X,68("="))')
           write(*,'(3X,A)')      'ERROR:  Missing input file'
           write(*,*)
           write(*,'(3X,3(A))')   'Input file ',trim(mol(i)%fname),    &
                                  ' not found in the current directory'
           write(*,'(2X,68("="))')
           call exit(0)
         end if
!
           read(uniinp,'(A)') mol(i)%title
           read(uniinp,*)     mol(i)%nat
!
         allocate(mol(i)%renum(mol(i)%nat)   ,  &
                  mol(i)%rename(mol(i)%nat)  ,  &
                  mol(i)%atname(mol(i)%nat)  ,  &
                  mol(i)%atnum( mol(i)%nat)  ,  &
                  mol(i)%mass(mol(i)%nat)    ,  &
                  mol(i)%coord(3,mol(i)%nat) )
!
         do k = 1, mol(i)%nat
           read(uniinp,'(I5,2A5,I5,3F8.3)') mol(i)%renum(k),   &
                                            mol(i)%rename(k),  &
                                            mol(i)%atname(k),  &
                                            mol(i)%atnum(k),   &
                                            mol(i)%coord(:,k)
         end do
! NO BORRRRRAR
!~            straux = sys(i)%rename(k)
!~            aux    = ''
!~            do
!~              select case ( straux(1:1) )
!~                case ( 'a':'z','A':'Z')
!~                  sys(i)%rename(k) = straux
!~                  read(aux,*) sys(i)%renum(k)
!~                  exit
!~                case default
!~                  aux    = trim(aux)//straux(1:1)
!~                  straux = straux(2:)
!~              end select
!~            end do
           read(uniinp,*) mol(i)%latvec
         close(uniinp)
       end do
!
       do i = 1, n
         do j = 1, mol(i)%nat
           aux    = adjustl(mol(i)%atname(j))
           straux = ''
           do
             select case ( aux(1:1) )
               case ( 'a':'z','A':'Z')
                 straux = trim(straux)//aux(1:1)
                 aux    = aux(2:)
               case default
                 exit
             end select
           end do
           select case ( straux )
             case ( 'H' )
               mol(i)%mass(j) = 1.007825d0
             case ( 'HE' )
               mol(i)%mass(j) = 4.002602d0  ! Not exact
             case ( 'LI' )
               mol(i)%mass(j) = 6.941d0     ! Not exact
             case ( 'BE' )
               mol(i)%mass(j) = 9.012182d0  ! Not exact
             case ( 'B' )
               mol(i)%mass(j) = 10.811d0    ! Not exact
             case ( 'C' )
               mol(i)%mass(j) = 12.0d0
             case ( 'N' )
               mol(i)%mass(j) = 14.003074d0
             case ( 'O' )
               mol(i)%mass(j) = 15.994915d0
             case ( 'F' )
               mol(i)%mass(j) = 18.998403d0 ! Not exact
             case ( 'NE' )
               mol(i)%mass(j) = 20.1797d0   ! Not exact
             case ( 'CL' )
               mol(i)%mass(j) = 35.453d0    ! Not exact
             case ( 'AR' )
               mol(i)%mass(j) = 39.948d0    ! Not exact
             case ( 'KR' )
               mol(i)%mass(j) = 83.798d0    ! Not exact
           end select
         end do
       end do
!
       return
       end subroutine read_gro
!
!======================================================================!
!
       subroutine randcoord(nfile,nmol,L,r0,ntrial,outp,debug)
!
       use datatypes
       use geometry
!
       implicit none
!
       include 'parameters.h'
       include 'info.h'
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(ndim),intent(in)   ::  L       !  Simulation box
       real(kind=8),intent(in)                   ::  r0      !  Minimum separation
       integer,dimension(nfile),intent(in)       ::  nmol    !  Number of molecules
       integer,intent(in)                        ::  nfile   !  Number of input files
       integer,intent(in)                        ::  ntrial  !  Maximum number of trials
       character(len=lenout),intent(in)          ::  outp    !  Output file name
       logical,intent(in)                        ::  debug   !  Debug mode
!
! Declaration of the local variables
!
       real(kind=8),dimension(ndim)              ::  r       !  Relative position
       real(kind=8)                              ::  dist    !  Minimum image distance
       real(kind=8)                              ::  rmin    !  Corrected miminum separation
       logical                                   ::  check   !  Checking variable
       integer                                   ::  ifile   !  Input file index
       integer                                   ::  jfile   !  Input file index
       integer                                   ::  imol    !  Molecule index
       integer                                   ::  jmol    !  Molecule index
       integer                                   ::  iat     !  Atom index
       integer                                   ::  jat     !  Atom index
       integer                                   ::  inext   !  System index
       integer                                   ::  jnext   !  System index
       integer                                   ::  itrial  !  Trial index
       integer                                   ::  irenum  !  Residue number
       integer                                   ::  iatnum  !  Atom number
!
! Generating the coordinates of the first atom
!
       call randmol(L,mol(1)%nat,mol(1)%coord,sys%coord(:,:mol(1)%nat))
!
       inext = 0
       do iat = 1, mol(1)%nat
         sys%renum(iat)  = mol(1)%renum(iat)
         sys%rename(iat) = mol(1)%rename(iat)
         sys%atname(iat) = mol(1)%atname(iat)
         sys%atnum(iat)  = mol(1)%atnum(iat)
         inext = inext + 1
       end do
!
       irenum = 1
!
       open(unit=uniout,file=trim(outp),action='write')

         write(uniout,*) 'Gromacs coordinates file'
         write(uniout,*) sys%nat
       do iatnum = 1, mol(1)%nat
         write(uniout,'(I5,2A5,I5,3F8.3)') irenum               ,  &
                                           mol(1)%rename(iatnum),  &
                                           mol(1)%atname(iatnum),  &
                                           iatnum               ,  &
                                           sys%coord(:,iatnum)
       end do
!
! Generating random coordinates
!
!~        inext = mol(1)%atnum(mol(1)%nat)
       do ifile = 1, nfile
!~ if (debug) write(*,'(A,I8)')       &
!~ 'Reference atom :',ifile
	     do imol = 1, nmol(ifile)
           if ( ( imol.eq.1 ) .and. (ifile.eq.1) ) cycle
           irenum = irenum + 1
!~ if (debug) write(*,'(A,2I8)')       &
!~ 'Reference atom :',ifile,imol
           if (debug) then
             write(*,*)
             write(*,*) 'Generating coordinates of molecule',ifile,imol
             write(*,'(1X,58("-"))')
           end if
! Initial settings
           itrial = 1
           rmin   = r0*r0
! While the distance between two atoms is lower than r0 generate new coordinates
          ! do while ( check )
999          continue
             if (debug) write(*,*) 'Generating new coordinates  : ',   &
                                    itrial
             if ( mod(itrial,ntrial+1) .eq. 0 ) then
               rmin   = 0.95d0*rmin	! FLAG
               itrial = 0
             end if
! Generating coordinates for the molecule imol of type ifile
             call randmol(L,mol(ifile)%nat,mol(ifile)%coord,           &
                          sys%coord(:,inext+1:inext+mol(ifile)%nat)  )
! Checking the distance of every new atom with respect to the previous particles
             do iat = 1, mol(ifile)%nat
!~ if (debug) write(*,'(2(A,3I8))')       &
!~ 'Reference atom :',ifile,imol,iat
               jnext = 0
               do jfile = 1, ifile 
!~ if (debug) write(*,'(A,3I8,A,I8)')       &
!~ 'Reference atom :',ifile,imol,iat,     &
!~ ' Target atom :',jfile
                 do jmol = 1, nmol(jfile)
                   if ( (ifile.eq.jfile) .and. (imol.le.jmol) ) cycle
!~ if (debug) write(*,'(A,3I8,A,2I8)')       &
!~ 'Reference atom :',ifile,imol,iat,     &
!~ ' Target atom :',jfile,jmol
                   do jat = 1, mol(jfile)%nat
!~ if (debug) write(*,'(2(A,3I8))')       &
!~ 'Reference atom :',ifile,imol,iat,     &
!~ ' Target atom :',jfile,jmol,jat
! Calculating the minimum image vector  !FLAG ERROR
                     r  = minimgvec(sys%coord(:,inext+iat), &
                                    sys%coord(:,jnext+jat),L)
! Calculating the minimum distance between two atoms
                     dist = sqrt(dot_product(r,r))
                     if (debug) write(*,'(2(A,3I8),2(A,F12.5))')    &
                                'Reference atom :',ifile,imol,iat,  &
                                ' Target atom :',jfile,jmol,jat,    &
                                ' Distance :',dist*10.0d0,          &
                                ' Reference :',rmin*10.0d0
! If there is any distance lower than r0 exit the loop and generate new coordinates
                     if ( dist .lt. rmin ) then
                       itrial = itrial + 1				
                       GO TO 999
                     end if
                   end do  ! end jat
                   jnext = jnext + mol(jfile)%nat
                 end do    ! end jmol
                 continue
!~                  jnext = jnext + mol(jfile)%nat
               end do      ! end jfile
             end do        ! end iat
        !   end do          ! end while
           do iatnum = 1, mol(ifile)%nat
             write(uniout,'(I5,2A5,I5,3F8.3)') irenum               ,  &
                                               mol(ifile)%rename(iatnum),  &
                                               mol(ifile)%atname(iatnum),  &
                                               inext + iatnum           ,  &
                                               sys%coord(:,inext+iatnum)
           end do
           inext = inext + mol(ifile)%nat
         end do            ! end imol
       end do              ! end ifile
         write(uniout,'(3F10.5)') L
!
       close(uniout)
!
       return
       end subroutine randcoord
!
!======================================================================!
!
       subroutine randmol(L,nat,inpmol,outmol)
!
       use datatypes
       use geometry
!
       implicit none
!
       include 'parameters.h'
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(ndim,nat),intent(in)   ::  inpmol  !  Input coordinates
       real(kind=8),dimension(ndim,nat),intent(out)  ::  outmol  !  Output coordinates
       real(kind=8),dimension(ndim),intent(in)       ::  L       !  Simulation box
!~        integer,intent(in)                            ::  ndim    !  Dimension of the Real Space
       integer,intent(in)                            ::  nat     !  Number of atoms
!
! Declaration of the local variables
!
       real(kind=8),dimension(ndim,nat)              ::  auxmol  !  Auxliary coordinates
       real(kind=8),dimension(ndim)                  ::  r       !  Translation vector
       real(kind=8)                                  ::  angle   !  Rotation angle
       integer                                       ::  iat     !  Atom index
integer :: k
!
!~            write(*,*) mol(1)%nat
!~  write(*,*) 'INITIAL MATRIX'
!~         do k = 1, mol(1)%nat
!~            write(*,'(A,3F8.4)') 'C',inpmol(:,k)*10.0d0
!~          end do
       call random_number(angle)
       angle = angle*2.0d0*pi
       call Rx(nat,inpmol,angle,outmol)
!~            write(*,*) mol(1)%nat
!~ write(*,*) 'Rx ROTATION'
!~          do k = 1, mol(1)%nat
!~            write(*,'(A,3F8.4)') 'C',outmol(:,k)*10.0d0
!~          end do
!
       call random_number(angle)
       angle = angle*2.0d0*pi
       call Ry(nat,outmol,angle,auxmol)
!~            write(*,*) mol(1)%nat
!~ write(*,*) 'Ry ROTATION'
!~          do k = 1, mol(1)%nat
!~            write(*,'(A,3F8.4)') 'C',auxmol(:,k)*10.0d0
!~          end do
!
       call random_number(angle)
       angle = angle*2.0d0*pi
       call Rz(nat,auxmol,angle,outmol)
!~            write(*,*) mol(1)%nat
!~ write(*,*) 'Rz ROTATION'
!~          do k = 1, mol(1)%nat
!~            write(*,'(A,3F8.4)') 'C',outmol(:,k)*10.0d0
!~          end do
!
       call random_number(r)
       r(:) = r(:)*L(:)
       call translate(nat,outmol,1.0d0,r,outmol)
!~            write(*,*) mol(1)%nat
!~ write(*,*) 'TRANSLATION'
!~          do k = 1, mol(1)%nat
!~            write(*,'(A,3F8.4)') 'C',outmol(:,k)*10.0d0
!~          end do
!~ write(*,*)
!
       return
       end subroutine randmol
!
!======================================================================!
