pro choose_color, index=index, red=red, green=green, blue=blue, $
                      rainbow=rainbow, ryb=ryb, ash_color=ash_color, yellow=yellow, mag_color=mag_color, $
                      quiet=quiet, mono_mag=mono_mag

       library_dir = '~/astro/idl/lib/'
       vctr=bytarr(256) & vctg=vctr & vctb=vctr

       if keyword_set(index) then begin
           loadct_mld, index
           if not keyword_set(quiet) then print, 'loading color spectrum ', index, ' with loadct_mld'
           tvlct,vctr,vctg,vctb, /get

       endif else begin
           if not keyword_set(rainbow) and not keyword_set(ryb) and not keyword_set(ash_color) and not keyword_set(yellow) then ash_color=1b
           if keyword_set(rainbow) then begin
               ;rainbow color spectrum, red to blue (Temp shell slices)
               openr,color_lun,library_dir+'ct41.bin', /get_lun
               if not keyword_set(quiet) then print, 'loading rainbow color spectrum'
           endif
           if keyword_set(ryb) then begin
               ;my color spectrum, red to blue (Temp shell slices)
               ;openr,color_lun,library_dir+'ct60.bin', /get_lun
               openr,color_lun,library_dir+'ct61.bin', /get_lun
               if not keyword_set(quiet) then print, 'loading red-yellow-blue color spectrum'
           endif

           if keyword_set(yellow) then begin
               openr,color_lun,library_dir+'ct49.bin', /get_lun
               if not keyword_set(quiet) then print, 'loading Thermal Wind spectrum (yellow)'
           endif

           if keyword_set(ash_color) then begin
               ;orange to blue spectrum (Vr shell slices)
               openr,color_lun,'ct42new.bin', /get_lun
               if not keyword_set(quiet) then print, 'loading ASH Vr color spectrum'
           endif
     
           readu,color_lun,vctr,vctg,vctb

           if keyword_set(mono_mag) then begin
              openr,color_lun,library_dir+'mono_mag_ct.bin', /get_lun
              readu,color_lun,vctr,vctg,vctb
              vctr(255) = 255
              vctg(255) = 255
              vctb(255) = 255
              vctr(0) = 0
              vctg(0) = 0
              vctb(0) = 0
           endif

           if keyword_set(mag_color) then begin
              openr,color_lun,library_dir+'mag_color.bin', /get_lun
              if not keyword_set(quiet) then print, 'loading magnetic color spectrum'
              readu,color_lun,vctr,vctg,vctb
              vctr(255) = 255
              vctg(255) = 255
              vctb(255) = 255
              vctr(0) = 0
              vctg(0) = 0
              vctb(0) = 0
           endif

           if keyword_set(ash_color) then begin
               vctb[1] = 1b
           endif

           tvlct,vctr,vctg,vctb

       endelse

       red = vctr
       green = vctg
       blue = vctb
       close, color_lun       
       free_lun, color_lun
end 
