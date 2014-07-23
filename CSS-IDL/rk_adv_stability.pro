function argument, c

   x = real_part(c)
   y = imaginary(c)
   arg = atan(y,x)
   return, arg
end

pro rk_adv_stability, nx
docolor
choose_color, /rain

xrk3 = sqrt(3d0)*dindgen(nx)/double(nx-1)
xrk4 = 2d0*sqrt(2d0)*dindgen(nx)/double(nx-1)
   
GRK3 = dcomplexarr(nx)
GRK4 = dcomplexarr(nx)

img = Complex(0d0,1d0)

GRK3 = 1d0-img*xrk3-0.5d0*xrk3^2+img*xrk3^3/6d0
GRK4 = 1d0-img*xrk4-0.5d0*xrk4^2+img*xrk4^3/6d0+xrk4^4/24d0

p1 = 1d0-abs(GRK4)
p2 = 1d0-abs(GRK3)
p3 = abs(abs(-argument(GRK4)/xrk4)-1d0)
p4 = abs(abs(-argument(GRK3)/xrk3)-1d0)

plot, xrk4, p1, thick=4, yrange=[1d-7,1d0], /ylog
oplot, xrk4, p1, color=20, thick=4
oplot, xrk3, p2, color=80, thick=4
oplot, xrk4, p3, color=150, linestyle=2, thick=4
oplot, xrk3, p4, color=220, linestyle=2, thick=4
oplot, xrk4, 0d0*xrk4

stop

end
