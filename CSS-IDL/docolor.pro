pro docolor, pseudocolor=pseudocolor

if (!d.name eq 'X') or (!d.name eq 'MAC') then begin
  if keyword_set(pseudocolor) then begin
    device, pseudo_color=8, retain=2, /install
  endif else begin
    device, true_color=24, retain=2, /install
  endelse
endif

device, decomposed=0

ctr=bytarr(256) & ctg=ctr & ctb=ctr
close,3 & openr,3,'ct42.bin'
readu,3,ctr,ctg,ctb
tvlct,ctr,ctg,ctb
close,3
end
