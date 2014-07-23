pro view_shell, image, result, lats, lons, xsize = xsize, ysize = ysize, pos=pos
xsizet = xsize
ysizet = ysize
if (not keyword_set(pos)) then pos = [0.16,0.16,0.84,0.84]

backval = 255b
map_set,lats(2),lons(2),limit=[lats(0),lons(0),lats(1),lons(1)], /noerase, /noborder, /ortho, pos=pos
;result=map_image(image,xone,yone,xsizet,ysizet,latmin=lats(0),latmax=lats(1),lonmin=lons(0),lonmax=lons(1),/bilinear,missing = backval)
result=map_proj_image(image,/bilinear,dimensions=[xsizet*0.68,ysizet*0.68],missing=backval)
xsr = n_elements(result(*,0))
ysr = n_elements(result(0,*))

dx = (xsize-xsr)/2
dy = (ysize-ysr)/2

result2 = bytarr(xsize,ysize)

result2(*,*) = 255b
result2(dx:dx+xsr-1,dy:dy+ysr-1) = result(*,*)

xsizet = xsize
ysizet = ysize

tv, result2, xsize = xsizet, ysize = ysizet,/NORM
if (lats(2)+11.0d0 gt lats(1)) then begin
   lats(1)=lats(1)+0.05d0
endif
;if (lons(2)+11.0d0 gt lons(1)) then begin
;   lons(1)=lons(1)+2.0d0
;endif
map_set,lats(2),lons(2),limit=[lats(0),lons(0),lats(1),lons(1)], /noerase, /noborder, /grid, londel = 20.0, latdel = 10.0,GLINESTYLE=0,GLINETHICK = 1.0, /ortho,pos=pos


END
