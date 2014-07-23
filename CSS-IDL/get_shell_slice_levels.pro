;+^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
;
;  get_contour_levels.pro - This function returns contour levels given
;                           an input field using a histogram trim
;                           to determine upper and lower bounds
;
;  usage:  levels = get_contour_levels(field,nlev,spacing=spacing,$
;                     mn=mn,mx=mx,pct=pct)
;       where levels = array of contour levels
;             field = input field
;             nlev = number of contour levels desired
;             spacing = on output, contains spacing between contour levels
;             mn = on output, contains (pct) histogram field value
;             mx = in output, contains (1-pct) histogram field value
;             pct = percentage trim (default=3)
;             avg_centered = if centering min and max, use average
;                            separation from zero level
;
;  M.DeRosa - 22 Oct 1999 - created
;  M.DeRosa -  2 Nov 1999 - added pct keyword
;  B.Brown  - 23 May 2006 - added avg_centered keyword
;                           and saturate_max, saturate_min
;
;-^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@sign.pro
function get_shell_slice_levels,field,lats,nlev,spacing=spacing,mn=mn,mx=mx,pct=pct, $
                                avg_centered=avg_centered, saturate_max=saturate_max, saturate_min=saturate_min, $
                                uneven_range=uneven_range


if not keyword_set(saturate_max) and not keyword_set(saturate_min) and $
   not keyword_set(avg_centered) and not keyword_set(uneven_range) then saturate_min=1b

;  preliminaries
nlev=fix(nlev(0))
if keyword_set(pct) then pct=pct(0) else pct=0.03

;  sort field and get statistics

cumulative_hist = spherical_histogram(field, lats, locations=bins, /latitude, nbins=100, /cumulative)
mn = bins(value_locate(cumulative_hist,   pct))
mx = bins(value_locate(cumulative_hist, 1-pct))
print, mn, mx
meanf=mean(field) ; this is still incorrect for our spherical geometry...

;  if range spans both positive and negative values, then set contours so that
;  there are an equal number of positive and negative contours
if (sign(mn) lt 0) and (sign(mx) gt 0) then begin
    if keyword_set(avg_centered) then begin
        mx=mean(abs([mn,mx]))
        mn=-mx
        cent=0.0
    endif
    if keyword_set(saturate_min) then begin
        mx=max(abs([mn,mx]))
        mn=-mx
        cent=0.0
    endif
    if keyword_set(saturate_max) then begin
        mx=min(abs([mn,mx]))
        mn=-mx
        cent=0.0
    endif
    if keyword_set(uneven_range) then cent=meanf
endif else cent=meanf

;  now determine spacing and calculate contour levels
spacing=(mx-mn)/(nlev-1.)
spacing=(fix(spacing/10.^fix(alog10(spacing)))+1)*10.^fix(alog10(spacing))
levels=linrange(nlev,cent-spacing*(0.5*(nlev-1)),cent+spacing*(0.5*(nlev-1)))
;print,'  number of contour levels:  '+strtrim(nlev,2)
;print,'  contour spacing:  ',spacing
;print,'  min,max,avg: ',[minmax(field),meanf]

return,levels
end
