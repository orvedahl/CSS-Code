;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
;+
;  pro display_shell_slice
;  Benjamin Brown
;  Oct 22nd, 2004
;
;  Summary: This is a routine I've been meaning to spin off for quite
;  some time.  It displays a single shell slice in either molweide or
;  orthographic projection.  This is based on pro shell_display,
;  currently a subroutine of shell_slices_to_jpg.pro.
;
;  I assume in this that the shell slice has already been read into
;  memeory by some previous routine.  This routine is merely a display
;  routine. 
;
;
;       usage: display_shell_slice, lat, lon, shell, $
;                         title=title, textu=textu, Mollweide=Mollweide, $ 
;                         publication=publication, maxv=maxv, minv=minv
;
;-
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@scim.pro

pro display_shell_slice, lat, lon, shell, $
                         Mollweide=Mollweide, Cylindrical=Cylindrical, Orthographic=Orthographic,$ 
                         title=title, textu=textu, subtitle=subtitle, $
                         publication=publication, science=science, maxv=maxv, minv=minv, $
                         rainbow=rainbow, quiet=quiet, color_bar=color_bar, $
                         psf=psf, latlon_window=latlon_window, grid=grid, format=format, $
                         white_bg=white_bg, center_on_zero=center_on_zero, $
                         no_interpolation=no_interpolation, $
                         radial_scale_factor=radial_scale_factor, $
                         transparent=transparent, res_scale=res_scale,highlightbox=highlightbox,nolines=nolines

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Constants

  dimensions_of_shell = size(shell)
  N_Phi = dimensions_of_shell[1]-1
  N_Theta = dimensions_of_shell[2]
  original_p_multi = !p.multi


  missing_data = !VALUES.F_NaN

  ; colors should be set by the program calling display shell slice

  ; determine contour levels / image scale
  ; do this automatically; used to query user as to whether they agreed
  ; with the choices.  Streamline by using auto-determined levels.
  nlevs=10                        ;  number of contour levels
  levs=get_shell_slice_levels(shell,lat,nlevs,mn=mnv,mx=mxv,pct=0.04)

  if not keyword_set(quiet) then begin
      if not keyword_set(maxv) AND not keyword_set(minv) then begin
          print, 'using auto-determined levels:  min = ', mnv, ' max = ', mxv
          print, '        values for the field:  min = ', min(shell), ' max = ', max(shell)
      endif else begin
          print, ' using pre-determined levels:  min = ', minv,' max = ', maxv
          print, '      auto-determined levels:  min = ', mnv, ' max = ', mxv
          print, '        values for the field:  min = ', min(shell), ' max = ', max(shell)
      endelse
  endif

  ; however, if the user passed min and max values, use those.
  if n_elements(maxv) gt 0 then mxv=maxv
  if n_elements(minv) gt 0 then mnv=minv
  print, 'scaling values are: ', mnv, mxv

  ; begin display procedures
  if keyword_set(Mollweide) then begin      
      if n_elements(Mollweide) lt 2 then begin
          lcent=0.0  &  bcent=0.0 ;; 35.0
      endif else begin
          lcent=Mollweide[0]  &  bcent=Mollweide[1]
      endelse

    ; **********************************************************
    ;          plot on Mollweide projection    
    ; **********************************************************

      if keyword_set(publication) then begin

          if keyword_set(color_bar) then begin
              ymargin = [1,1] ;For horizontal bar use ymargin=[9,1] and xmargin=[1,1]
              xmargin = [5,1]
          endif else begin
              ymargin = [1,1]
              xmargin = [1,1]
          endelse
          map_set,bcent,lcent,/mollweide,/noer,/nobor,/iso,chars=1.0,$
            ymargin=ymargin, xmargin=xmargin, title=title, limit=latlon_window

          print, 'sized for publication'
          
      endif else begin
          map_set,bcent,lcent,/mollweide,/noer,/nobor,/iso,chars=1.8,$
            ymargin=[6,10],title=title, limit=latlon_window
      endelse

      clip=!p.clip              ; get pictures boundary x0,y0,x1,y1
      winsiz=double([!d.x_size,!d.y_size])
      original_pos=(clip)/[winsiz,winsiz]


      ;; here's a hell of a kludge to make sure we get shells to scale
      ;; properly with radius.  Map_set updates !p.clip, but not
      ;; !p.pos, so if we want to remain flexible with respect to
      ;; x and y margins and !p.multi, we need to cleverly use !p.clip
      ;; to infer !p.position, set this in a position call, and then
      ;; undo our damage later on.  BPB 2/12/2007 
      if keyword_set(radial_scale_factor) then begin
          map_clip = !p.clip
          map_x_size = map_clip[2]-map_clip[0]
          map_x_size_prime = radial_scale_factor*map_x_size
          map_y_size = map_clip[3]-map_clip[1]
          map_y_size_prime = radial_scale_factor*map_y_size
          map_clip[0] = map_clip[0] + 0.5*(map_x_size-map_x_size_prime)
          map_clip[1] = map_clip[1] + 0.5*(map_y_size-map_y_size_prime)
          map_clip[2] = map_clip[2] - 0.5*(map_x_size-map_x_size_prime)
          map_clip[3] = map_clip[3] - 0.5*(map_y_size-map_y_size_prime)
          
          map_pos=map_clip/[winsiz,winsiz]
          map_set,bcent,lcent,/mollweide,/noer,/nobor,/iso,chars=1.0,$
            limit=latlon_window, pos=map_pos
      endif


      if keyword_set(publication) then color_pos = [0.40,0.06,0.60,0.09] $;[0.30,0.06,0.70,0.09]
                                  else color_pos = [0.30,0.05,0.70,0.08] 

   
    endif else if keyword_set(cylindrical) then begin

    ; **********************************************************
    ;           Instead, plot on a lat-long map
    ; **********************************************************        
        if n_elements(Cylindrical) lt 2 then begin
            lcent=0.0  &  bcent=0.0
        endif else begin
            lcent=Cylindrical[0]  &  bcent=Cylindrical[1]
        endelse

        ;  set up coordinate system
        map_set,bcent,lcent,/cylindrical,/noer,/nobor,/iso,chars=1.0,$
          ymargin=[1,1], xmargin=[1, 1], limit=latlon_window ;, subtitle=subtitle
        xyouts, 0.75, 0.85, title,align=2,/norm
        
                                ;  remap and display image on sphere
        
        if keyword_set(publication) then color_pos = [0.40,0.15,0.60,0.18] $
                                    else color_pos = [0.30,0.05,0.70,0.08]
    endif else if keyword_set(Orthographic) then begin
    ; **********************************************************
    ;          plot on Orthographic (Spherical) projection    
    ; **********************************************************

        if n_elements(Orthographic) lt 2 then begin
            lcent=0.0  &  bcent=10.0 
        endif else begin
            lcent=Orthographic[0]  &  bcent=Orthographic[1]
        endelse

        ymargin = [1,2]
        xmargin = [1,1]

        if keyword_set(publication) then begin
            map_set,bcent,lcent,/orthographic,/noer,/nobor,/iso,chars=1.0,$
              ymargin=ymargin, xmargin=xmargin, title=title, limit=latlon_window ;, subtitle=subtitle
            print, 'sized for publication'
        endif else begin
            map_set,bcent,lcent,/orthographic,/noer,/nobor,/iso,chars=1.8,$
              ymargin=[6,10],title=title, limit=latlon_window ;, subtitle=subtitle
        endelse

        clip=!p.clip            ; get pictures boundary x0,y0,x1,y1
        winsiz=double([!d.x_size,!d.y_size])
        original_pos=(clip)/[winsiz,winsiz]


      ;; here's a hell of a kludge to make sure we get shells to scale
      ;; properly with radius.  Map_set updates !p.clip, but not
      ;; !p.pos, so if we want to remain flexible with respect to
      ;; x and y margins and !p.multi, we need to cleverly use !p.clip
      ;; to infer !p.position, set this in a position call, and then
      ;; undo our damage later on.  BPB 2/12/2007 
        if keyword_set(radial_scale_factor) then begin
            map_clip = !p.clip
            map_x_size = map_clip[2]-map_clip[0]
            map_x_size_prime = radial_scale_factor*map_x_size
            map_y_size = map_clip[3]-map_clip[1]
            map_y_size_prime = radial_scale_factor*map_y_size
            map_clip[0] = map_clip[0] + 0.5*(map_x_size-map_x_size_prime)
            map_clip[1] = map_clip[1] + 0.5*(map_y_size-map_y_size_prime)
            map_clip[2] = map_clip[2] - 0.5*(map_x_size-map_x_size_prime)
            map_clip[3] = map_clip[3] - 0.5*(map_y_size-map_y_size_prime)
            
            ymargin = [1,2]
            xmargin = [1,1]
            map_pos=map_clip/[winsiz,winsiz]
            map_set,bcent,lcent,/orthographic,/noer,/nobor,/iso,chars=1.0,$
              ymargin=ymargin, xmargin=xmargin, limit=latlon_window, pos=map_pos
            
        endif

        if keyword_set(publication) then color_pos = [0.40,0.05,0.60,0.08] $
        else color_pos = [0.30,0.05,0.70,0.08]

   
  endif

    ; **********************************************************
    ;           Regardless, do the following...
    ; **********************************************************        

      if n_elements(transparent) ne 0 then colorbar_text_color = byte(1)

        ;  remap and display image on sphere
      if keyword_set(psf) then $ 
        backval = byte(255) $                         ; white background for PS files
      else backval=byte(0)                            ; standard black background

      if keyword_set(white_bg) then begin
          backval=byte(255)                           ; make a white background on request
          colorbar_text_color = byte(0)               ; with black text for the colorbar
      endif

      scaled_shell = shell[0:n_phi-1,*]
      if keyword_set(no_interpolation) then begin
          remaped_shell=map_image(scaled_shell,x0,y0,xsiz,ysiz,latmin=lat(0),latmax=lat(n_theta-1),$
                          lonmin=lon(0),lonmax=lon(n_phi-1),missing=missing_data,$
                          /compress,scale=res_scale)
      endif else begin
          remaped_shell=map_image(scaled_shell,x0,y0,xsiz,ysiz,latmin=lat(0),latmax=lat(n_theta-1),$
                          lonmin=lon(0),lonmax=lon(n_phi-1),/bilinear,missing=missing_data,$
                          /compress,scale=res_scale)
      endelse


      ; rescale input image for display, see scim.pro
      scim,remaped_shell,scale=[mnv,mxv],outim=scaled_remap,/nowin,/quiet,/interp, zero_value=byte(255) ;255b ; original      

      print, size(scaled_remap)

      tv,scaled_remap,x0,y0,xsiz=xsiz,ysiz=ysiz, /device

      If (not keyword_set(nolines)) Then Begin
         ;  put in equator and redraw horizon
         oplot,findgen(361)-lcent-90.,replicate(0,361),linestyle=5,thick=2,color=colorbar_text_color
      EndIf

      if keyword_set(grid) then map_grid, londel=45.0, latdel=30.0, color=colorbar_text_color, glinethick=1.0, glinestyle=2
      if keyword_set(highlightbox) then begin
         map_grid, lats=[highlightbox(0),highlightbox(2)], lons=[highlightbox(1),highlightbox(3)], $
                   color=colorbar_text_color, glinethick=3, glinestyle=0
      endif

      If (not keyword_set(nolines)) Then Begin
         ;nd = n_elements(lat)/6
         ;latgrid = [lat(0),lat(nd),lat(2*nd),lat(3*nd),lat(4*nd),lat(5*nd),lat(6*nd)]
         ;nd = n_elements(lon)/6
         ;longrid = [lon(0),lon(nd),lon(2*nd),lon(3*nd),lon(4*nd),lon(5*nd),lon(6*nd)]
         ;map_grid,lats=latgrid,lons=longrid,thick=2,/horizon
         map_grid,thick=2,/horizon
      EndIf

      if keyword_set(radial_scale_factor) then begin
          map_set,bcent,lcent,mollweide=mollweide, orthographic=orthographic,/noer,/nobor,/horizon,/iso,chars=1.0,$
            ymargin=ymargin, xmargin=xmargin, title=title, limit=latlon_window, pos=original_pos
      endif

  print, 'displaying shell slice, about to make colorbar'
    if keyword_set(color_bar) then begin
        if keyword_set(publication) then begin
            ;; output a publication grade colourbar.  Mark the top and bottom
            ;; with rounded values of the min and max.  Call it good.
            
            pos_colorbar_side = [0.02,0.3,0.04,0.7]
            if max([abs(mnv),abs(mxv)]) gt 1 and max([abs(mnv),abs(mxv)]) lt 1d6 then begin
                colorbar_format='(I8)'
                color_min = discrete_round(mnv,1)
                color_max = discrete_round(mxv,1)
            endif else begin
                colorbar_format='(G10.3)'
                color_min = mnv
                color_max = mxv
            endelse
               colourbar2,pos_colorbar_side,mn=color_min, mx=color_max, format=colorbar_format, $
               zero=colourbar_zero,textu=textu, /right_zero_tick, text_color=colorbar_text_color, text_char_size=2.0
            ;  colourbar, color_pos, mn=mnv,mx=mxv,textu=textu, text_color=colorbar_text_color ;, format=format 
        endif else begin
            ;; do a simple horizontal colourbar
            colourbar, color_pos, mn=mnv,mx=mxv,textu=textu, text_color=colorbar_text_color ;, format=format       
        endelse
    endif

    if keyword_set(science) then begin
    ; display PDF of shell quantities

        ; save plotting variables
         Map_P = !P
         Map_X = !X
         Map_Y = !Y

         hist_pos = [0.05, 0.05, 0.2, 0.15]
         shell_hist = spherical_histogram(shell, lat, nbins=100, locations=hist_bins, /latitude)
         non_zero_values =  where(shell_hist gt 0)

         plot, hist_bins[non_zero_values], shell_hist[non_zero_values], $
           pos=hist_pos, /normal, /noerase, /ylog, ys=1, xs=3, charsize=0.5

         ; overlay the color table min and max cuts
         if mnv gt min(hist_bins) then plots, [mnv, mnv], [min(shell_hist[non_zero_values]), max(shell_hist)], linestyle=1
         if mxv lt max(hist_bins) then plots, [mxv, mxv], [min(shell_hist[non_zero_values]), max(shell_hist)], linestyle=1

         ; print out the actual min and max values in the shell
         x_out = hist_pos[2]+0.01
         y_out = 0.5*(hist_pos[1]+hist_pos[3])
         xyouts, x_out, y_out - 0.0125, align=0, 'shell min: ' + string(min(shell), format='(g10.4)'), /norm, charsize=0.75
         xyouts, x_out, y_out + 0.0125, align=0, 'shell max: ' + string(max(shell), format='(g10.4)'), /norm, charsize=0.75
         ; restore original plot coordinates
         !P = Map_P
         !X = Map_X
         !Y = Map_Y
     endif
    !p.multi = original_p_multi
    !p.multi[0] = !p.multi[0]+1
end
