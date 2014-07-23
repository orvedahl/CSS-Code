function ud6, dr, arr, order

If (order eq 1) Then Begin
    np = n_elements(arr)
    a0 = max([max(arr),1d0])
    array = arr/a0
    der = dblarr(np)
;Sixth order accurate first derivatives at the ends
    left = (-147.0d0*array(0)+360.0d0*array(1)-450.0d0*array(2)+400.0d0*array(3)-225.0d0*array(4)+72.0d0*array(5)-10.0d0*array(6))/60.0d0
    left1 = (-10.0d0*array(0)-77.0d0*array(1)+150.0d0*array(2)-100.0d0*array(3)+50.0d0*array(4)-15.0d0*array(5)+2.0d0*array(6))/60.0d0
    left2 = (2.0d0*array(0)-24.0d0*array(1)-35.0d0*array(2)+80.0d0*array(3)-30.0d0*array(4)+8.0d0*array(5)-1.0d0*array(6))/60.0d0
    
    right2 = (1.0d0*array(np-7)-8.0d0*array(np-6)+30.0d0*array(np-5)-80.0d0*array(np-4)+35.0d0*array(np-3)+24.0d0*array(np-2)-2.0d0*array(np-1))/60.0d0
    right1 = (-2.0d0*array(np-7)+15.0d0*array(np-6)-50.0d0*array(np-5)+100.0d0*array(np-4)-150.0d0*array(np-3)+77.0d0*array(np-2)+10.0d0*array(np-1))/60.0d0
    right = (10.0d0*array(np-7)-72.0d0*array(np-6)+225.0d0*array(np-5)-400.0d0*array(np-4)+450.0d0*array(np-3)-360.0d0*array(np-2)+147.0d0*array(np-1))/60.0d0
    
    der(0)=left
    der(1)=left1
    der(2)=left2
    
    der(np-1)=right
    der(np-2)=right1
    der(np-3)=right2
    
;Sixth order accurate first derivatives in the interior
    
    for i=3,np-4 do begin
        der(i) = 1.0d0*array(i+3)-9.0d0*array(i+2)+45.0d0*array(i+1)-45.0d0*array(i-1)+9.0d0*array(i-2)-1.0d0*array(i-3)
        der(i) = der(i)/60.0d0
    endfor

    der = der*a0/dr
Endif Else If (order eq 2) Then Begin
    np = n_elements(arr)
    a0 = max([max(arr),1d0])
    array = arr/a0
    der = dblarr(np)

;Sixth order accurate first derivatives in the interior
    left = (812d0*array(0)-3132d0*array(1)+5265d0*array(2)-5080d0*array(3)+2970d0*array(4)-972d0*array(5)+137d0*array(6))/180d0
    left1 = (137d0*array(0)-147d0*array(1)-255d0*array(2)+470d0*array(3)-285d0*array(4)+93d0*array(5)-13d0*array(6))/180d0
    left2 = (-13d0*array(0)+228d0*array(1)-420d0*array(2)+200d0*array(3)+15d0*array(4)-12d0*array(5)+2d0*array(6))/180d0
    for i=3,np-4 do begin
        der(i) = 2d0*array(i+3)-27d0*array(i+2)+270d0*array(i+1)-490d0*array(i)+270d0*array(i-1)-27d0*array(i-2)+2d0*array(i-3)
        der(i) = der(i)/180d0
    endfor

    right = (812d0*array(np-1)-3132d0*array(np-2)+5265d0*array(np-3)-5080d0*array(np-4)+2970d0*array(np-5)-972d0*array(np-6)+137d0*array(np-7))/180d0
    right1 = (137d0*array(np-1)-147d0*array(np-2)-255d0*array(np-3)+470d0*array(np-4)-285d0*array(np-5)+93d0*array(np-6)-13d0*array(np-7))/180d0
    right2 = (-13d0*array(np-1)+228d0*array(np-2)-420d0*array(np-3)+200d0*array(np-4)+15d0*array(np-5)-12d0*array(np-6)+2d0*array(np-7))/180d0

    der(0)=left
    der(1)=left1
    der(2)=left2
    
    der(np-1)=right
    der(np-2)=right1
    der(np-3)=right2

    der = der*a0/dr^2

EndIf

return, der

end
