       subroutine calsoft_read_codes

       use calibration_data_module
       use plant_data_module
       use input_file_module
       use soil_module
       use plant_module
       use hydrograph_module
       use hru_lte_module
       use sd_channel_module
       use organic_mineral_mass_module
       use mgt_operations_module
       use conditional_module
       use utils, only: open_file

       implicit none

       character (len=80) :: titldum = ""!           |title of file
       character (len=80) :: header = "" !           |header of file
       integer :: eof = 0              !           |end of file

       eof = 0

       if (open_file(107, in_chg%codes_sft)) then
         do
           read (107,*,iostat=eof) titldum
           if (eof < 0) exit
           read (107,*,iostat=eof) header
           if (eof < 0) exit
           read (107,*,iostat=eof) cal_codes
           if (eof < 0) exit
           exit
         enddo
         close(107)

         if (cal_codes%hyd_hru /= "n" .or. cal_codes%hyd_hrul == "y".or.    &
             cal_codes%plt == "y" .or. cal_codes%sed == "y" .or.            &
             cal_codes%nut == "y" .or. cal_codes%chsed == "y" .or.          &
             cal_codes%chnut == "y" .or. cal_codes%res == "y") cal_soft = "y"
       end if

       return
      end subroutine calsoft_read_codes
