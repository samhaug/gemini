      subroutine stripmn(cmstripped,cmnotstr,clength)
c
c Cuts away the pathname from the model file name read in in 
c the main program.
c

      integer i,clength
      character cmstripped*(*),cmnotstr*(*)

      do i=clength,1,-1
        if (cmnotstr(i:i).eq.'/') goto 10
      enddo

 10   cmstripped=cmnotstr((i+1):(i+10))

      end
