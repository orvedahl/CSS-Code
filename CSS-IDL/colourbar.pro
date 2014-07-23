;******************************************************************
   pro colourbar,pos,mx=mx,mn=mn,zero=zero,plttype=plttype,textu=textu,$
                 textmn=textmn,textmx=textmx, text_char_size=text_char_size, format=format, $
                 font=font, text_color=text_color, text_thick=text_thick, pmvals=pmvals

;------------------------------------------------------------------
; pos=[x0,y0,x1,y1]
;
 
   ;original_font = !p.font
   ;if keyword_set(font) then !p.font = font
   ;!p.font=3
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

   bar=(indgen(ncols)+1)#replicate(1,fix(abs(y1-y0)))

   bar=map(bar,abs(x1-x0)+2,abs(y1-y0))+2

   tv,bytscl(bar,top=ncols),pos(0),pos(1),xs=xp,ys=yp,/norm

   thick=2

   p0 = pos(0)*!d.x_size
   p1 = pos(1)*!d.y_size
   p2 = pos(2)*!d.x_size
   p3 = pos(3)*!d.y_size
   y_char = !D.Y_CH_SIZE*text_char_size
   x_char = !D.X_CH_SIZE*text_char_size

   plots,[p0+2,p2+1,p2+1,p0+2,p0+2],[p1,p1,p3-1,p3-1,p1], $
         thick=thick,/dev, color=text_color;,color=2

;
; Scaling fact from PS to regular window device ...
 if plttype eq 'ps' then sf = 1/32. else sf = 1.
;
; A letter width is about 7 pixels ...
sf=1
  if n_elements(textu) eq 0 then textu=' '

if keyword_set(pmvals) then begin
   val = max([abs(mx),abs(mn)])
   if (val eq abs(mx)) then begin
      valtxt = textmx
   endif else begin
      valtxt = strsplit(textmn,'-',/EXTRACT)      
   endelse
   xyouts,0.5*(p2+p0),p1-1.25*y_char/sf,textoidl('\pm')+valtxt,/dev, charsize=text_char_size, color=text_color, align=0.5, charthick=text_thick
endif else begin
  ; max first
  if (n_elements(textmx) ne 0) then begin
      xyouts,p2+1.5*x_char,0.5*(p1+p3)-0.5*y_char,textmx,/dev, charsize=text_char_size, color=text_color, align=0.0, charthick=text_thick
      ;xyouts, p2,p1-1.25*y_char/sf, textmx, /dev, charsize=text_char_size, color=text_color, align=0.5, charthick=text_thick
  endif else begin  
      if n_elements(mx) eq 0 then begin
          xyouts,p2,p1-1.25*y_char/sf,'!6max',/dev, charsize=text_char_size, color=text_color, align=0.5, charthick=2.0
      endif else begin
              xyouts,p2,p1-1.25*y_char/sf,strtrim(string(mx,form=format),2),/dev,$
                charsize=text_char_size, color=text_color, align=0.5, charthick=2.0
;              xyouts,p2,p1-1.25*y_char,' '+textu,/dev,$
;                charsize=text_char_size, color=text_color, align=0.0, charthick=2.0
      endelse
  endelse


  ; output min marker
  if (n_elements(textmn) ne 0) then begin
      xyouts,p0-1.5*x_char,0.5*(p1+p3)-0.5*y_char,textmn,/dev, charsize=text_char_size, color=text_color, align=1.0, charthick=text_thick
      ;xyouts, p0,p1-1.25*y_char/sf, textmn, /dev, charsize=text_char_size, color=text_color, align=0.5, charthick=text_thick
  endif else begin  
      if n_elements(mn) eq 0 then begin 
          xyouts,p0,p1-1.25*y_char/sf,'!6min',/dev, charsize=text_char_size, color=text_color, align=0.5, charthick=2.0
      endif else begin
              xyouts,p0,p1-1.25*y_char/sf,strtrim(string(mn,form=format),2),/dev, $
                charsize=text_char_size, color=text_color, align=0.5, charthick=2.0
;              xyouts,p0,p1-1.25*y_char,' '+textu,/dev, $
;                charsize=text_char_size, color=text_color, align=0.0, charthick=2.0
      endelse
  endelse
endelse
; evaluate the position of the central tick (i.e colors are not necessary symetric)
; do not output anything if min and max values have not been passed;
; instead make only a 'qualitative' colourbar.
if not keyword_set(pmvals) then begin
  if keyword_set(mn) and keyword_set(mx) then begin
      if not keyword_set(zero) then begin
          if mn lt 0 then begin
              zero=0.
          endif else begin
              zero=(mx+mn)/2.
          endelse
      endif
  
      posc=(p2-p0)/(mx-mn)*zero+p0-(p2-p0)/(mx-mn)*mn
      
      plots,[posc,posc],[p1,p1-y_char/2.0],thick=thick,/dev, color=text_color ;,color=2
  endif
endif
  ;!p.font = original_font
end




