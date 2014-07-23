;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
;
;  scim.pro - (Re)scales input image for display
;
;  usage:  scim,im,mag=mag,win=win,scale=scale,top=top,bot=bot,outim=outim,$
;            /interp,/nowin,/quiet
;          where im = image to be scaled
;                mag = magnification factor (default=1)
;                win = window number to open (default=0)
;                scale = image scaling: if not set, then image scale is
;                        [min(im),max(im)]; if set to a scalar, then image
;                        scale is [-scale,scale]; if set to a vector, then
;                        image scale is [scale(0),scale(1)]
;                top = maximum value of scaled result (default=255)
;                bot = minimum value of scaled result (default=0)
;                outim = variable into which to put output image
;                interp = for enlargement, use interpolation rather than
;                         pixel replication if magnification factor is
;                         integral (for non-iontegral magnification
;                         factors, interpolation is always used)
;                nowin = if set, no window is opened (usually used in
;                        conjunction with outim keyword)
;                quiet = it set, output printed to screen is suppressed
;
;  M.DeRosa -  4 Nov 1998 - added this intro
;              4 Nov 1998 - added capability for non-integral magnifications
;             20 Aug 1999 - nowin keyword supersedes opw and pixmap
;                           keywords
;             24 Aug 1999 - changed default enlargement to pixel
;                           replication rather than interpolation for
;                           integral magnification factors
;             11 Dec 1999 - added quiet keyword
;
;  B.Brown
;             30 Nov 2005 - changed default floor and ceiling for
;                           scaling; this lets us put special colors
;                           at the first and last index (black and
;                           white typically) of the color table 
;
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

pro scim,im,mag=mag,win=win,interp=interp,scale=scale,top=top,bot=bot,$
         outim=outim,nowin=nowin,quiet=quiet, auto_levels=auto_levels, $
         missing_data=missing_data,zero_value=zero_value

;  usage message
if n_params() eq 0 then begin
  print,'  scim,im,mag=mag,win=win,scale=scale,top=top,bot=bot,outim=outim,'+$
    '/interp,/nowin,/quiet'
  return
endif

;  get image parameters
szim=size(im)
if szim(0) ne 2 then begin
  print,'  scim.pro:  input image must have exactly two dimensions'
  return
endif
nax=szim(1:2)
byt=szim(3) eq 1
if byt then immin=fix(min(im)) else immin=min(im)
if byt then immax=fix(max(im)) else immax=max(im)

;  check keywords and set parameters
if keyword_set(interp) then pixrep=0 else pixrep=1
if keyword_set(nowin) then nowin=1 else nowin=0
if n_elements(mag) eq 0 then mag=1.0 else mag=float(mag)
if n_elements(win) eq 0 then win=0 else win=fix(win)
if n_elements(top) eq 0 then top=254b else top=byte(top(0))
if n_elements(bot) eq 0 then bot=1b else bot=byte(bot(0))
case n_elements(scale) of
  0:  begin  ;  scale keyword is not set
    mn=immin
    mx=immax
    end
  1:  begin  ;  scale keyword is set to a scalar
    mn=-scale(0)
    mx=scale(0)
    end
  else:  begin  ;  scale keyword is set to a vector
    mn=scale(0)
    mx=scale(1)
    end
endcase

if keyword_set(auto_levels) then begin
  ; determine contour levels / image scale
  ; do this automatically.

  nlevs=10                        ;  number of contour levels
  levs=get_contour_levels(im,nlevs,mn=mn,mx=mx,pct=0.04)

  scale = [mn, mx]
endif    

;  random output
if not keyword_set(quiet) then begin
  print,'  scim:  original image min/max is '+strcompress(immin,/r)+'/'+$
    strcompress(immax,/r)
  print,'  scim:  display scale is '+strcompress(mn,/r)+' to '+$
    strcompress(mx,/r)
endif

;  now display image
naxm=nax*mag
if mag mod 1 ne 0.0 then outim=congrid(im,naxm(0),naxm(1),/interp) $
  else outim=rebin(im,nax(0)*mag,nax(1)*mag,sample=pixrep)

outim=bytscl(outim,min=mn,max=mx,top=top-bot,/NaN)+bot


if keyword_set(missing_data) then begin
    if finite(missing_data) then missing_data = where(im eq missing_data, missing_data_count) $
                            else good_data = where(finite(im), complement=missing_data, ncomplement=missing_data_count)
endif else begin
    good_data = where(finite(im), complement=missing_data, ncomplement=missing_data_count)
endelse
if missing_data_count ge 1 then outim[missing_data] = zero_value    

if not nowin then begin
  window,win,xs=fix(naxm(0)),ys=fix(naxm(1))
  tv,outim
endif

end
