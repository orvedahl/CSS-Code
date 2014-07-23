; idem colourbar but vertical
;******************************************************************

@map.pro
   pro colourbar2,pos,mx=mx,mn=mn,zero=zero,plttype=plttype,textu=textu,$
       textmn=textmn,textmx=textmx, text_char_size=text_char_size, format=format, charthick=charthick
;------------------------------------------------------------------
; pos=[x0,y0,x1,y1]
;
   if not keyword_set(format) then format='(g9.3)'
   if keyword_set(text_char_size) then begin
       text_char_size = text_char_size
   endif else begin
       text_char_size = 1.4
   endelse

   if n_elements(plttype) eq 0 then plttype = ' '

   xp = abs(pos(2)-pos(0))
   yp = abs(pos(3)-pos(1))

   xx0=convert_coord(pos(0),pos(1),/normal,/to_device)
   xx1=convert_coord(pos(2),pos(3),/normal,/to_device)
;   ncols = !d.n_colors-1 ; original version
   ncols=255

   x0 = xx0(0)
   y0 = xx0(1)
   x1 = xx1(0)
   y1 = xx1(1)

;   bar=(indgen(ncols)+1)#replicate(1,fix(abs(y1-y0))) ; previous
   bar=replicate(fix(abs(y1-y0)),1)#(indgen(ncols)+1)
   bar=map(bar,abs(x1-x0)+2,abs(y1-y0))+2

   tv,bytscl(bar,top=ncols),pos(0),pos(1),xs=xp,ys=yp,/norm

   thick=2
   p0 = pos(0)*!d.x_size
   p1 = pos(1)*!d.y_size
   p2 = pos(2)*!d.x_size
   p3 = pos(3)*!d.y_size
   y_char = !D.Y_CH_SIZE
   plots,[p0+2,p2+1,p2+1,p0+2,p0+2],[p1,p1,p3-1,p3-1,p1], $
         thick=thick,/dev;,color=2
;
; Scaling fact from PS to regular window device ...
  if plttype eq 'ps' then sf = 1/32. else sf = 1.
;
; A letter width is about 7 pixels ...

  if n_elements(textu) eq 0 then textu=' '

;
;  Output min and max
;

  if (n_elements(textmx) ne 0) then begin
      xyouts, p0+(p2-p0)/2.,p3+10/sf,textmx, /dev, charsize=text_char_size, align=0.5, charthick=charthick
  endif else begin 
      if n_elements(mx) eq 0 then xyouts,p2,p3+5/sf,'!6max',/dev else $
        xyouts,p2,p3+y_char/sf/2., strtrim(string(mx,form=format),2), $
        /dev, charsize=text_char_size,charthick=1.4, align=1.0 ;,color=2
        xyouts,p2,p3+y_char/sf/2., ' '+textu, $
        /dev, charsize=text_char_size,charthick=1.4, align=0.0 ;,color=2

        ; original outputs number and units as a single string
        ;xyouts,p0+(p2-p0)/2,p3+6/sf, strtrim(string(mx,form=format),2)+textu, $
        ;   /dev, charsize=text_char_size,charthick=1.4, align=0.5 ;,color=2
  endelse 


  if (n_elements(textmn) ne 0) then begin
     xyouts, p0+(p2-p0)/2.,p1-15/sf,textmn, /dev, charsize=text_char_size, align=0.5, charthick=charthick
  endif else begin 
      if n_elements(mn) eq 0 then begin
          xyouts,p0+10/sf,p1-10/sf,'!6min',/dev 
      endif else begin
        xyouts,p2, p1-1.5*y_char/sf, strtrim(string(mn,form=format),2), $
          /dev, charsize=text_char_size,charthick=1.4, align=1.0 ;,color=2
        xyouts,p2, p1-1.5*y_char/sf, ' '+textu, $
          /dev, charsize=text_char_size,charthick=1.4, align=0.0 ;,color=2

        ; original outputs number and units as a single string
        ;xyouts,p0+(p2-p0)/2, p1-17/sf, strtrim(string(mn,form=format),2)+textu, $
        ;/dev, charsize=text_char_size,charthick=1.4, align=0.5 ;,color=2
    endelse
  endelse




; evaluate the position of the central tick (i.e colors are not necessary symetric)
; do not output anything if min and max values have not been passed;
; instead make only a 'qualitative' colourbar.
  if keyword_set(mn) and keyword_set(mx) then begin
      if not keyword_set(zero) then begin
          if mn lt 0 then begin
              zero=0.
          endif else begin
              zero=(mx+mn)/2.
          endelse
      endif
      posc=(p3-p1)/(mx-mn)*zero+p1-(p3-p1)/(mx-mn)*mn
      ; make the tick mark
      plots,[p2,p2+5/sf],[posc,posc],thick=thick,/dev ;,color=2
  endif 
end




