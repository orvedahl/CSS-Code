PRO TextBox_Event, event
      ; This event handler responds to all events. Widget
      ; is always destoyed. The text is recorded if ACCEPT
      ; button is selected or user hits CR in text widget.

   Widget_Control, event.top, Get_UValue=info

   CASE event.ID OF
      info.cancelID: Widget_Control, event.id, /Destroy
      ELSE: BEGIN

            ; Get the text and store it in the pointer location.
         Widget_Control, info.textID, Get_Value=theText
         (*info.ptr).text = theText[0]
         (*info.ptr).cancel = 0 ; The user hit ACCEPT.
         Widget_Control, event.top, /Destroy
         ENDELSE
   ENDCASE
END

FUNCTION TextBox, Title=title, Label=label, Cancel=cancel, $
   Group_Leader=groupleader, XSize=xsize, Value=value

      ; Return to caller if there is an error. Set the cancel
      ; flag and destroy the group leader if it was created.

   Catch, theError
   IF theError NE 0 THEN BEGIN
      Catch, /Cancel
      ok = Dialog_Message(!Error_State.Msg)
      IF destoy_groupleader THEN Widget_Control, groupleader, /Destroy
      cancel = 1
      RETURN, ""
   ENDIF

      ; Check parameters and keywords.

   IF N_Elements(title) EQ 0 THEN title = 'Provide Input:'
   IF N_Elements(label) EQ 0 THEN label = ""
   IF N_Elements(value) EQ 0 THEN value = ""
   IF N_Elements(xsize) EQ 0 THEN xsize = 200

   IF N_Elements(groupleader) EQ 0 THEN BEGIN
      groupleader = Widget_Base(Map=0)
      Widget_Control, groupleader, /Realize
      destroy_groupleader = 1
   ENDIF ELSE destroy_groupleader = 0

      ; Create modal base widget.

   tlb = Widget_Base(Title=title, Column=1, /Modal, $
      /Base_Align_Center, Group_Leader=groupleader)

   labelbase = Widget_Base(tlb, Row=1)
   IF label NE "" THEN label = Widget_Label(labelbase, Value=label)
   textID = Widget_Text(labelbase, /Editable, Scr_XSize=xsize, Value=value)
   buttonBase = Widget_Base(tlb, Row=1)
   cancelID = Widget_Button(buttonBase, Value='Cancel')
   acceptID = Widget_Button(buttonBase, Value='Accept')

      ; Center the widgets on display.

   Device, Get_Screen_Size=screenSize
   IF screenSize[0] GT 2000 THEN screenSize[0] = screenSize[0]/2 ; Dual monitors.
   xCenter = screenSize(0) / 2
   yCenter = screenSize(1) / 2
   geom = Widget_Info(tlb, /Geometry)
   xHalfSize = geom.Scr_XSize / 2
   yHalfSize = geom.Scr_YSize / 2
   Widget_Control, tlb, XOffset = xCenter-xHalfSize, YOffset = yCenter-yHalfSize

      ; Realize the widget hierarchy.

   Widget_Control, tlb, /Realize

   ptr = Ptr_New({text:"", cancel:1})

   info = {ptr:ptr, textID:textID, cancelID:cancelID}
   Widget_Control, tlb, Set_UValue=info, /No_Copy

   XManager, 'TextBox', tlb

   theText = (*ptr).text
   cancel = (*ptr).cancel
   Ptr_Free, ptr
   IF destroy_groupleader THEN Widget_Control, groupleader, /Destroy

   RETURN, theText
END
