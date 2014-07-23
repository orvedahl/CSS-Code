pro image_contour, a, WINDOW_SCALE = window_scale, ASPECT = aspect, $
                   INTERP = interp,title=title,xtitle=xtitle, $
                   ytitle=ytitle, xx=xx,yy=yy,levels=levels,c_labels=c_labels,$
                   charsize=charsize, data_overlay=data_overlay,$
                   noaxis=noaxis, overplot=overplot, c_thick=c_thick,$
                   mini=mini,maxi=maxi,blackline=blackline, $
                   position=position, nocontour=nocontour

;+
; NAME:
;	IMAGE_CONT
;
; PURPOSE:
;	Overlay an image and a contour plot.
;
; CATEGORY:
;	General graphics.
;
; CALLING SEQUENCE:
;	IMAGE_CONT, A
;
; INPUTS:
;	A:	The two-dimensional array to display as an image.
;       data_overlay:     the two-dimensional array of the same dimensions as a to overlay 
;                            (default is A itself) as a contour plot.
;
; KEYWORD PARAMETERS:
; WINDOW_SCALE:	Set this keyword to scale the window size to the image size.
;		Otherwise, the image size is scaled to the window size.
;		This keyword is ignored when outputting to devices with 
;		scalable pixels (e.g., PostScript).
;
;	ASPECT:	Set this keyword to retain the image's aspect ratio.
;		Square pixels are assumed.  If WINDOW_SCALE is set, the 
;		aspect ratio is automatically retained.
;
;	INTERP:	If this keyword is set, bilinear interpolation is used if 
;		the image is resized.
;
; OUTPUTS:
;	No explicit outputs.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	The currently selected display is affected.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	If the device has scalable pixels, then the image is written over
;	the plot window.
;
; MODIFICATION HISTORY:
;	DMS, May, 1988.
;       Add overplot keyword, 9/17/99
;
;       B.Brown    Dec  2, 2005
;                  Cleaning up, standardizing image scaling.
;                  Jun 12, 2006 
;                  Adding a position keyword
;-

;on_error,2                      ;Return to caller if an error occurs
sz = size(a)			;Size of image
if sz(0) lt 2 then message, 'Parameter not 2D'

	;set window used by contour

contour,[[0,0],[1,1]],/nodata, xstyle=4, ystyle = 4, pos=position

px = !x.window * !d.x_vsize	;Get size of window in device units
py = !y.window * !d.y_vsize
swx = px(1)-px(0)		;Size in x in device units
swy = py(1)-py(0)		;Size in Y
six = float(sz(1))		;Image sizes
siy = float(sz(2))
aspi = six / siy		;Image aspect ratio
aspw = swx / swy		;Window aspect ratio
f = aspi / aspw			;Ratio of aspect ratios
print, 'f=', f
if not keyword_set(title) then title=' '
if not keyword_set(xtitle) then xtitle=' '
if not keyword_set(ytitle) then ytitle=' '
if not keyword_set(c_labels) then c_labels=[0]
if not keyword_set(c_thick) then c_thick=[1]  ; modif add c_thick
if not keyword_set(charsize) then charsize=1.0
if not keyword_set(data_overlay) then data_overlay=a

if not keyword_set(blackline) then blackline=0b 

mx = !d.n_colors-1		;Brightest color

if keyword_set(blackline) then begin
    contour_colors=0
endif else begin
    contour_colors = mx
endelse

if (!d.flags and 1) ne 0 then begin ;Scalable pixels? i.e. ps or psf
    ;; **************************************************************************
    ;; *              Postscript display follows this branch                    *     
    ;; **************************************************************************

    background_color = !d.n_colors-1 ; = mx ; this is the brightest color in the scaling

    if not keyword_set(mini) then mini=min(a)
    if not keyword_set(maxi) then maxi=max(a)

    background_pixels=!VALUES.F_NAN ;min(a)  ; by default, image_polar places missing pixels at a
                              ; value corresponding to min(a), far below the 
                              ; floor of the actual data.
        
    invert_index = where(a eq background_pixels, invert_count)

    scim,a,scale=[mini,maxi],outim=a_image,/nowin,/quiet,/interp,zero_value=background_color ; original

    print, 'resetting background colors at ', invert_count, ' pixels'

    a_image[invert_index] = background_color
 

    if keyword_set(aspect) then begin ;Retain aspect ratio?
        ;Adjust window size to maintain x/y aspect ratio
        if f ge 1.0 then swy = swy / f else swx = swx * f
    endif

    tv,a_image,px(0),py(0),xsize=swx,ysize=swy,/device 
        
    if (keyword_set(noaxis)) then begin
        xstyle=5 & ystyle=5
    endif else begin
        xstyle=1 & ystyle=1
    endelse

    ;contour,a,/noerase,xst=xstyle,yst=ystyle,$ ;Do the contour 
    ;  pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
    ;  c_color = contour_colors,levels=levels,title=title,xtitle=xtitle, $
    ;  ytitle=ytitle,c_labels=c_labels,c_thick=c_thick,charsize=charsize,/overplot
    
endif else begin           
    ;; **************************************************************************
    ;; *                  Window display follows this branch                    *     
    ;; **************************************************************************

    background_color = 0

    if not keyword_set(mini) then mini=min(a)
    if not keyword_set(maxi) then maxi=max(a)

    background_pixels=!VALUES.F_NAN ;min(a)  ; by default, image_polar places missing pixels at a
                              ; value corresponding to min(a), far below the 
                              ; floor of the actual data.
        
    invert_index = where(a eq background_pixels, invert_count)

    scim,a,scale=[mini,maxi],outim=a_image,/nowin,/quiet,/interp,zero_value=background_color ; original
    print, 'resetting background colors at ', invert_count, ' pixels'

    a_image[invert_index] = background_color

    if keyword_set(window_scale) then begin 
        ;Scale the window to the image size
        tv,a_image,px(0),py(0)	;Output image
        swx = six		;Set window size from image
        swy = siy
    endif else begin
        ;Scale the image to the current window size
        ;if keyword_set(aspect) then begin
           ;Maintain an aspect ratio?
            if f ge 1.0 then swy = swy / f else swx = swx * f
        ;endif  ;aspect
        tv,poly_2d(a_image,$   ;Resample and expand to match the window size.  Do this in scim?
                   [[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],$
                   keyword_set(interp),swx,swy), $
          px(0),py(0) 
    endelse			;window_scale
endelse                         ;scalable pixels

if (keyword_set(levels)) then begin
    if (keyword_set(noaxis)) then begin
       xstyle=5 & ystyle=5
    endif else begin
       xstyle=1 & ystyle=1
    endelse
    if keyword_set(overplot) then begin
        If (keyword_set(nocontour)) Then c_thick=0
        if keyword_set(xx) then begin
            contour,data_overlay,xx,yy,/noerase,xst=xstyle,yst=ystyle,$ ;Do the contour 
              pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
              c_color = contour_colors,levels=levels,title=title,xtitle=xtitle, $
              ytitle=ytitle,c_labels=c_labels,charsize=charsize,c_thick=c_thick,/overplot
        endif else begin
            contour,data_overlay,/noerase,xst=xstyle,yst=ystyle,$ ;Do the contour 
              pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
              c_color = contour_colors,levels=levels,title=title,xtitle=xtitle, $
              ytitle=ytitle,c_labels=c_labels,charsize=charsize,c_thick=c_thick,/overplot
        endelse
    endif else begin
        If (keyword_set(nocontour)) Then c_thick=0
        if keyword_set(xx) then begin
            contour,data_overlay,xx,yy,/noerase,xst=xstyle,yst=ystyle,$ ;Do the contour 
              pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
              c_color = contour_colors,levels=levels,title=title,xtitle=xtitle, $
              ytitle=ytitle,c_labels=c_labels,c_thick=c_thick,charsize=charsize
        endif else begin
            contour,data_overlay,/noerase,xst=xstyle,yst=ystyle,$ ;Do the contour 
              pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
              c_color = contour_colors,levels=levels,title=title,xtitle=xtitle, $
              ytitle=ytitle,c_labels=c_labels,c_thick=c_thick,charsize=charsize
        endelse
    endelse                    ; overplot;
endif

return
end






