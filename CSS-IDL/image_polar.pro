;+
;   image_polar.pro
;   ???
;   Sacha Allan Brun (?)
;
;   May 22, 2006  Ben Brown -- I've modified this a lot over the past
;                              year, trying to comment and clean up
;                              code.  It continues to evolve.
;
;   May 22, 2006  Ben Brown -- Major cleanup and rewrite of the border
;                              drawing section.  Old code is in
;                              image_polar_bak.pro
;
;   July 2006, Kyle Augustson -- A few minor changes (removed
;                                unnecessary contour mess at the boundaries)
;
;   May 2008, Kyle Augustson -- Major revisions for CSS, it is now
;                               maleable enough to plot differing
;                               latitude bounds
;
;-
@image_contour.pro

pro image_polar,a,r,costheta,spacing=spacing,amap=amap,title=title,$
                vmin=vmin,noborder=noborder,$
                levels=levels,c_labels=c_labels,c_thick=c_thick, $
                mini=mini,maxi=maxi,$
                blackline=blackline, $
                position=position,$
                special_level=special_level, partial_equator=partial_equator, $
                interface=interface, nointerp=nointerp, $
                publication=publication, bound=bound, nocontour=nocontour, window_scale=window_scale

   theta = !pi*.5 - acos(costheta)
   nth = n_elements(theta)
   nr = n_elements(r)

   rmax = max(r)
   rmin = min(r)

   if (not keyword_set(halfdomain)) Then halfdomain=0
   If (not keyword_set(bound)) Then bound=[0,-rmax,rmax,rmax]

   if not keyword_set(spacing) then spacing = 1d-2*rmax
   if not keyword_set(title) then title = ' '

   ;If (keyword_set(vmin)) Then Begin
   ;    vmin = 0d0
   ;    amin = 0d0
   ;Endif Else Begin
   ;    amin = min(a)
   ;    if (amin gt 0d0) then amin=max(a)+1d0
   ;Endelse
   amin = !VALUES.F_NAN

   if (not keyword_set(maxi)) then maxi = max(a)
   if (not keyword_set(mini)) then mini = min(a)

   margin = 0.0
   bound=[rmin*min(cos(theta)),rmax*min(sin(theta)),rmax,rmax*max(sin(theta))]
   print, bound

   If (not keyword_set(nointerp)) Then Begin
       ntot = n_elements(a)
       radii = dblarr(ntot)
       thetas = dblarr(ntot)
       plarr = dblarr(ntot)
       dr = (rmax-rmin)/double(nr)/100d0
       dth = (max(theta)-min(theta))/double(nth)/200d0
       randsr = (2d0*randomu(1312431L,ntot)-1d0)*dr
       randst = (2d0*randomu(2312432L,ntot)-1d0)*dth
       k=0L
       for it=0,nth-1 do begin
           for ir=0,nr-1 do begin
               plarr[k]=a[ir,it]
               radii[k]=r[ir]+randsr[k]
               thetas[k]=theta[it]+randst[k]
               k = k+1L
           endfor
       endfor

       amap = polar_surface(plarr,radii,thetas,missing=amin,spacing=spacing*[1,1],/quintic)
      ;stop
   endif else begin
       plarr = a
       amap = polar_surface(plarr,r,theta,missing=amin,bound=bound,spacing=spacing*[1,1],/quintic, /grid)
   endelse

   print, 'spacing is ', spacing
   sz=size(amap) & nx=sz(1) & ny=sz(2)
   x=amap & y=x
   for i=long(0),ny-1 do x(*,i) = (bound(2)-bound(0))*dindgen(nx)/double(nx-1)+bound(0)
   for i=long(0),nx-1 do y(i,*) = (bound(3)-bound(1))*dindgen(ny)/double(ny-1)+bound(1)
   rad = sqrt(x^2+y^2)
   ind=where(rad lt rmin,comp=cind)
   amap[ind] = amin 
 
   if (not keyword_set(levels)) then begin
      nlevels=10
      levels=(0.95*max(a)-0.95*min(a))*dindgen(nlevels)/(nlevels-1d0)+1d4*max(a)
   endif

   image_contour,amap,xx=findgen(nx),yy=findgen(ny),title=title,levels=levels,/noaxis,$
                 mini=mini,maxi=maxi, pos=position, nocontour=nocontour, window_scale=window_scale ; original

   if (not keyword_set(noborder)) then begin
      ; now we need to do everything in pixel space

      nym = abs(sin(theta[nth-1]))*double(ny)/(abs(sin(theta[0]))+abs(sin(theta[nth-1])))
      x = (rmax*max(cos(theta))-rmin*min(cos(theta)))*dindgen(nx)/double(nx-1) + rmin*min(abs(cos(theta)))

      ; draw the equator
      If (theta[nth-1] lt 0d0) Then Begin
         ix3 = max(where(x le rmin))
         plots,[nx*(x[ix3]-x[0])/(x[nx-1]-x[0]), nx],[nym, nym],linestyle=2,thick=2 ;, color=back_color
      Endif

      ;Draw latitudinal boundaries
      ;Upper
      s1 = tan(theta[0])/2d0
      ix1 = min(where(x ge rmin*cos(theta[0])))
      ix2 = min(where(x ge rmax*cos(theta[0])))
      ytop = double(ny)*s1*(double(nx)-dindgen(nx)-1d0)/double(nx-1)+1
      xp = indgen(nx)
      oplot, xp[ix1:ix2], ny-ytop[ix1:ix2], thick=2

      ;Lower
      s2 = tan(theta[nth-1])/2d0
      ix1 = min(where(x ge rmin*cos(theta[nth-1])))
      ix2 = min(where(x ge rmax*cos(theta[nth-1])))
      ybot = double(ny)*s2*(double(nx)-dindgen(nx)-1d0)/double(nx-1)
      oplot, xp[ix1:ix2], 2-ybot[ix1:ix2], thick=2

      ;Radial boundaries
      ;Lower
      If (theta[nth-1] lt 0d0) Then Begin
         nyb = ny-nym
         ix1 = 0
         ix3 = min(where(x ge rmin))
         y = sqrt(x[ix3]^2-x[ix1:ix3]^2)*double(nyb)/sqrt(x[ix3]^2-x[ix1]^2)-ytop[ix1]
         nny = n_elements(y)
         oplot, xp[ix1:ix3+1], nym+[y,0], thick=2 ;Portion above eq
         ix1 = min(where(x ge rmin*cos(theta[nth-1])))
         y = sqrt(x[ix3]^2-x[ix1:ix3]^2)*double(nym)/sqrt(x[ix3]^2-x[ix1]^2)+ybot[ix1]
         oplot, xp[ix1:ix3+1], nym-[y,0], thick=2 ;Portion below eq
      Endif Else Begin
         ix1 = min(where(x ge rmin*cos(theta[0])))
         ix2 = min(where(x ge rmin*cos(theta[nth-1])))
         y = sqrt(rmin^2-x^2)/sqrt(rmin^2-x[0]^2)*double(ny)
         oplot, xp[ix1:ix2], y[ix1:ix2], thick=2
      EndElse

      ;Upper
      If (theta[nth-1] lt 0d0) Then Begin
         ix3 = nx-1
         nyb = ny-nym
         ix1 = min(where(x ge rmax*cos(theta[0])))
         y = sqrt(rmax^2-x[ix1:ix3]^2)*double(nyb)/sqrt(rmax^2-x[ix1]^2)
         nny = n_elements(y)
         oplot, xp[ix1:ix3], nym+[y[0:nny-2],0d0], thick=2 ;Portion above eq
         ix1 = min(where(x ge rmax*cos(theta[nth-1])))
         y = sqrt(rmax^2-x[ix1:ix3]^2)*double(nym)/sqrt(rmax^2-x[ix1]^2)
         oplot, xp[ix1:ix3], nym-y, thick=2 ;Portion below eq
      Endif Else Begin
         ix1 = min(where(x ge rmax*cos(theta[0])))
         ix2 = min(where(x ge rmax*cos(theta[nth-1])))
         y = sqrt(rmax^2-x^2)/sqrt(rmax^2-x[0]^2)*double(ny)
         oplot, xp[ix1:ix2], y[ix1:ix2], thick=2
      EndElse

      if (keyword_set(interface)) then begin
         if (interface gt rmin) then begin
            interface_edge = interface/rmax*outer_edge
            y_interface=sqrt(interface_edge^2 - x^2)
            oplot,x,y_middle + y_interface,line=2,thick=2 ;, color=back_color
            oplot,x,y_middle - y_interface,line=2,thick=2 ;, color=back_color
         endif
      endif
   endif
;stop
end


