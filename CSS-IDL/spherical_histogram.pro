;+
;     spherical_histogram.pro
;     Benjamin Brown
;     May 26, 2006
;     
;     This IDL function calculates a latitude (or co-latitude)
;     weighted histogram, for use particularly with shell slice
;     outputs from ASH.
;
;          function spherical_histogram, shell, theta, latitude=latitude, locations=locations,$
;                                        binsize=binsize, nbins=nbins
;     
;     Feb 14, 2007 Just caught a subtle error in implementation of the angular weighting function!
;
;-
function spherical_histogram, shell, theta, latitude=latitude, locations=locations,$
                              binsize=binsize, nbins=nbins, cumulative=cumulative

     if not keyword_set(shell) then begin
         print, 'incorrect calling sequence in spherical_histogram'
         stop
     endif
     if not keyword_set(theta) then begin
         print, 'when using spherical_histogram you must pass in the angular distribution of the data array'
         stop
     endif

     shell_max = max(shell)
     shell_min = min(shell)


     if keyword_set(latitude) then angular_weight = cos(theta*!pi/180.) $
     else angular_weight = sin(theta)
     
     if not keyword_set(binsize) and not keyword_set(nbins) then nbins=100
     if keyword_set(nbins) then begin
         binsize = (shell_max-shell_min)/(nbins)
     endif

     n_lats = n_elements(theta)
     working_hist = fltarr(nbins)
     hist = fltarr(nbins)
     
     for i=0, n_lats-1 do begin
         ; IDL already lets you do iterative histogram building!
         ; but will that work with our sin(theta) weighting?
         ;print, working_histogram
         working_hist = histogram(shell[*,i], binsize=binsize, min=shell_min, max=shell_max)*angular_weight[i]
         hist = hist + working_hist

     endfor
     ;hist = hist/(nbins*n_lats)
     

     bins = binsize*findgen(n_elements(hist)) + shell_min
     
     ;hist = hist/int_tabulated(bins, hist)
     locations = bins
     if keyword_set(cumulative) then begin
         cumulative_hist = total(hist, /cumulative)
         ; do a clever trapazoidal rule, by Geoff Vasil's suggestion:
         ;
         ; f[i] = 1/2 f[0] + total(f[1:i-1]) + 1/2 f[i]
         ;
         cumulative_hist = cumulative_hist - 0.5*(hist[0] + hist)
         ; now normalize
         normalized_cumulative_hist = cumulative_hist/cumulative_hist(n_elements(cumulative_hist)-1)
         return, normalized_cumulative_hist
     endif
     ; else return the histogram itself
     return,hist
end
