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
;                         Mollweide=Mollweide, $ 
;                         publication=publication, maxv=maxv, minv=minv
;
;-
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

pro display_css_shell_slice, lat, lon, shell, $
                         Mollweide=Mollweide, Cylindrical=Cylindrical, Orthographic=Orthographic,$ 
                         publication=publication, science=science, maxv=maxv, minv=minv, $
                         psf=psf, latlon_window=latlon_window, grid=grid, format=format, $
                         white_bg=white_bg, $
                         no_interpolation=no_interpolation, $
                         radial_scale_factor=radial_scale_factor, $
                         res_scale=res_scale, nolines=nolines, pos=pos

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
              ymargin=ymargin, xmargin=xmargin, limit=latlon_window
            print, 'sized for publication'
        endif else begin
            map_set,bcent,lcent,/orthographic,/noer,/nobor,/iso,chars=1.8,$
              ymargin=[6,10], limit=latlon_window
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

    ;  remap and display image on sphere

      if keyword_set(psf) then $ 
        backval = byte(255) $                         ; white background for PS files
      else backval=byte(0)                            ; standard black background

      if keyword_set(white_bg) then begin
          backval=byte(255)                           ; make a white background on request
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

      ;Map lines
      ;Lats
      dx = xsiz*cos(max(lat)*!dpi/180d0)
      ddx = (xsiz-dx)/2d0
      plots, /dev, x0 + [ddx,xsiz-ddx], y0+[0,0], thick=2
      plots, /dev, x0 + [ddx,xsiz-ddx], y0+ysiz+[0,0], thick=2

      ;Meridians
      x = x0 + xsiz/2d0 - 0.5d0*xsiz*cos(lat*!dpi/180d0)
      y = ysiz*(1d0-(lat-min(lat))/(max(lat)-min(lat))) + y0
      plots, /dev, x, y, thick=2

      x = x0 + xsiz/2d0 + 0.5d0*xsiz*cos(lat*!dpi/180d0)
      plots, /dev, x, y, thick=2

      print, 'minmax(x)=', minmax(x)

      ;Equator
      plots, /dev, x0 + [0,xsiz], y0+ysiz/2+[0,0], thick=2, linestyle=2

end
