#!/usr/bin/env python
#
# Plot CSS Shell Slices as a GUI
#
#
# 7-3-2014 R. Orvedahl
#

#--------------------------------------------------------------------
import matplotlib
matplotlib.use('WXAgg')   # enable wxWidget backend stuff
matplotlib.use('PS')      # generate PS/EPS figures as default
import os
import sys
import getopt
import wx
import numpy
import string
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx
#--------------------------------------------------------------------
from read_shell_slice import *
import defaults
import constants
import derivatives
#--------------------------------------------------------------------


#####################################################################
# GUI class
#####################################################################
class window(wx.Frame):

    def __init__(self, parent, title, tex=True):

        self.parent = parent

        self.tex = tex

        if (self.tex):
            matplotlib.rc('text', usetex=True)
        else:
            matplotlib.rc('text', usetex=False)

        # return monitor dimensions
        xscreen, yscreen = wx.DisplaySize()

        # Application window size
        self.ys = defaults.fact_y*yscreen
        self.xs = defaults.fact_x*self.ys

        # (0,0) --> top left
        posx = defaults.pos_fact_x*xscreen
        posy = defaults.pos_fact_y*yscreen

        # initialize the title and size
        wx.Frame.__init__(self, parent, title=title, size=(self.xs, self.ys),
                    pos=(posx, posy),
                    style=wx.MINIMIZE_BOX|wx.MAXIMIZE_BOX|wx.RESIZE_BORDER|
                    wx.SYSTEM_MENU|wx.CAPTION|wx.CLOSE_BOX|wx.CLIP_CHILDREN)

        # setup Main window with a plot pane and a buttons pane
        self.mainpane = wx.Panel(self)
        self.plotpane = plotpanel(self.mainpane, self.xs, self.ys)
        self.buttons = buttonpanel(self.mainpane)

        # define a few attributes
        self.DefineAttributes()

        # make a menu
        self.SetMainMenuBar()

        # set sizers for the whole application window and one for the 
        # plot window and one for the button window
        self.mainsize = wx.BoxSizer(wx.HORIZONTAL)
        self.plotsize = wx.BoxSizer(wx.HORIZONTAL)
        self.buttonsize = wx.BoxSizer(wx.HORIZONTAL)

        # add the plot window and buttons window to the sizers
        self.plotsize.Add(self.plotpane, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        self.buttonsize.Add(self.buttons, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)

        # set 3 status bars
        self.statusbar = self.CreateStatusBar()
        self.statusbar.SetFieldsCount(3)
        self.statusbar.SetStatusWidths([-1,-1,-3])

        # this can be called anytime to change the status bar text with a 
        # little care as to which "parent" does the call
        self.statusbar.SetStatusText("Welcome", 0)
        self.statusbar.SetStatusText("Current Quantity: ", 1)
        self.statusbar.SetStatusText("Current File: ", 2)

        # add both sub-windows to main window
        self.mainsize.Add(self.plotsize, 3, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        self.mainsize.Add(self.buttons, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)

        # center the window on the screen
        #self.Centre()

        self.mainpane.SetSizer(self.mainsize)
        self.mainsize.Layout()

    def Radius(self, event):

        # the 100 is so it is not confused with ID for Quantity menu
        id = event.GetId() - 100

        # select radius
        if (id not in range(self.nradii)):
            self.ThrowError("Radius Not Available")

        else:
            self.current_rad = id
            self.plotpane.plot()

    def Quantity(self, event):

        id = event.GetId()

        if (id not in self.quantities.keys()):
            self.ThrowError("Quantity Undefined")
        else:
            quant = self.quantities[id]
            self.statusbar.SetStatusText("Quantity: " + quant, 1)

            self.current_qnt = id

            self.plot_title = "Shell Slice of "+quant
            self.plot_cmap_label = ""
            self.plot_vmin = None
            self.plot_vmax = None
            self.slider_min = 0.0
            self.slider_max = 100.0
            self.plot_enthalpy = False
            self.plot_vorticity = False
            self.plotpane.plot()

    def Annotate(self, event):

        self.plotpane.annotate_plot()
        
    def OnQuit(self, event):
        self.plotpane.Destroy()
        self.buttons.Destroy()
        self.mainpane.Destroy()
        self.Close()

    def ThrowError(self, message):
        msg = wx.MessageDialog(None, message, 'ERROR', wx.OK|wx.ICON_ERROR)
        msg.ShowModal()

    def ThrowWarning(self, message, cancel=False):
        if (cancel):
            msg = wx.MessageDialog(None, message, 'WARNING', 
                               wx.OK|wx.CANCEL|wx.ICON_EXCLAMATION)
            if (msg.ShowModal() == wx.ID_OK):
                return True
            else:
                return False
        else:
            msg = wx.MessageDialog(None, message, 'WARNING', 
                               wx.OK|wx.ICON_EXCLAMATION)
            msg.ShowModal()

    def DefineAttributes(self):

        self.data = None           # hold data to be displayed (many-D)
        self.data_fen = None       # hold enthalpy data
        self.data_vort = None      # hold vorticity data

        self.quantity_arr = None   # array of quantity indices
        self.nq = 0                # number of quantities
        self.quantities = {}       # dictionary of ID:Name pairs
        self.quantity_names = None # names of quantities

        self.radii = None          # hold radii data
        self.theta = None          # hold theta data
        self.phi = None            # hold phi data
        self.nradii = 0            # number of elements in radii data
        self.nth = 0               # number of elements in theta data
        self.np = 0                # number of elements in phi data

        self.current_qnt = 0       # hold index of current quantity
        self.current_rad = 0       # hold index of current radial slice
        self.current_rec = 0       # hold index of current record

        self.dirname = ""          # hold directory
        self.filename = ""         # hold filename
        self.fullfilename = ""     # hold full path to file

        self.slider_min = 0.0      # values of the sliders
        self.slider_max = 100.0    # values of the sliders

        # plot attributes
        self.plot_title = "Shell Slice"
        self.plot_xlabel = "phi (deg)"
        self.plot_ylabel = "theta (deg)"
        self.plot_xmin = None
        self.plot_xmax = None
        self.plot_ymin = None
        self.plot_ymax = None
        self.plot_cmap = defaults.cmap
        self.plot_rev_cmap = False
        self.plot_cmap_label = ""
        self.plot_vmin = None
        self.plot_vmax = None
        self.plot_manual_colorbar_range = False
        self.plot_manual_slider = False
        self.plot_enthalpy = False
        self.plot_vorticity = False
        self.plot_save_figure = False
        self.plot_save_output = "shell_slice_frame"

    def SetMainMenuBar(self):

        # setup a menubar (hard coded IDs that are different make it easier
        #                  to handle each event)
        self.mainmenuBar = wx.MenuBar()

        fileMenu = wx.Menu()
        fileMenu.Append(-10, '&Load File')
        fileMenu.Append(-20, '&Plot Data')
        fileMenu.Append(-30, '&Clear Plot')

        quit_item = wx.MenuItem(fileMenu, wx.ID_EXIT, '&Quit\tCtrl+w')
        fileMenu.AppendItem(quit_item)

        # empty on initialization
        quantity = wx.Menu()
        quantity.Append(-60, 'No File Loaded')

        # empty on initialization
        radii = wx.Menu()
        radii.Append(-70, 'No File Loaded')

        # for a submenu
        #fileMenu.AppendMenu(wx.ID_ANY, '&Select Quantity', quantity)

        self.mainmenuBar.Append(fileMenu, '&File')
        self.mainmenuBar.Append(quantity, '&Select Quantity')
        self.mainmenuBar.Append(radii, '&Radius')
        self.SetMenuBar(self.mainmenuBar)

        # make each object do something
        self.Bind(wx.EVT_MENU, self.buttons.Load, id=-10)
        self.Bind(wx.EVT_MENU, self.buttons.Plot, id=-20)
        self.Bind(wx.EVT_MENU, self.buttons.Clear_Plot, id=-30)
        self.Bind(wx.EVT_MENU, self.OnQuit, quit_item)
        self.Bind(wx.EVT_MENU, self.Quantity, id=-60)
        self.Bind(wx.EVT_MENU, self.Radius, id=-70)

    def UpdateMenu(self):

        # change Quantities Menu
        newquantityQ = wx.Menu()

        for i in range(self.nq):
            newquantityQ.Append(i, self.quantities[i])
            self.Bind(wx.EVT_MENU, self.Quantity, id=i)

        # replace old "Select Quantity" menu with updated quantities
        self.mainmenuBar.Replace(1, newquantityQ, '&Select Quantity')

        # change Radii Menu
        newquantityR = wx.Menu()

        for i in range(100,self.nradii+100):
            newquantityR.Append(i, str(self.radii[i-100]/constants.rsol))
            self.Bind(wx.EVT_MENU, self.Radius, id=i)

        # replace old "Radius" menu with updated radii
        self.mainmenuBar.Replace(2, newquantityR, '&Radius')

#####################################################################
# setup the button panel
#####################################################################
class buttonpanel(wx.Panel):

    def __init__(self, parent):

        wx.Panel.__init__(self, parent, wx.ID_ANY, style=wx.SUNKEN_BORDER)

        self.parent = parent

        self.mainparent = self.GetParent().GetParent()

        self.sizer = wx.BoxSizer(wx.VERTICAL)

        # setup a few buttons
        #--------------------------------------------------------------
        self.Lfile = wx.Button(self, -1, "Load File")
        self.Lfile.Bind(wx.EVT_BUTTON, self.Load)
        self.sizer.Add(self.Lfile, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.vortBut = wx.Button(self, -1, "Compute Vorticity")
        self.vortBut.Bind(wx.EVT_BUTTON, self.Vorticity)
        self.sizer.Add(self.vortBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.enthBut = wx.Button(self, -1, "Compute F_en")
        self.enthBut.Bind(wx.EVT_BUTTON, self.Compute_Enthalpy)
        self.sizer.Add(self.enthBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.animBut = wx.Button(self, -1, "Animate over Radii")
        self.animBut.Bind(wx.EVT_BUTTON, self.Animate)
        self.sizer.Add(self.animBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.annotateBut = wx.Button(self, -1, "Change Plot Attributes")
        self.annotateBut.Bind(wx.EVT_BUTTON, self.Annotate_Plot)
        self.sizer.Add(self.annotateBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.cmapBut = wx.Button(self, -1, "Change Colormap")
        self.cmapBut.Bind(wx.EVT_BUTTON, self.Change_Colormap)
        self.sizer.Add(self.cmapBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.quitBut = wx.Button(self, -1, "Quit")
        self.quitBut.Bind(wx.EVT_BUTTON, self.Quit)
        self.sizer.Add(self.quitBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------

        # setup some sliders
        #self.sizer2 = wx.BoxSizer(wx.VERTICAL)
        #--------------------------------------------------------------
        self.slider_rec_default = 0
        self.slider_r = wx.Slider(self, -1, self.slider_rec_default, 0, 100, 
                                wx.DefaultPosition,
                                wx.DefaultSize, wx.SL_LABELS|
                                wx.SL_HORIZONTAL|wx.SL_TOP, name="Record")
        self.slider_r.Bind(wx.EVT_SCROLL_CHANGED, self.SliderRec)
        self.sizer.Add(self.slider_r, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.slider_min_default = 0
        self.slider_m = wx.Slider(self, -1, self.slider_min_default, 0, 100, 
                                wx.DefaultPosition,
                                wx.DefaultSize, wx.SL_LABELS|
                                wx.SL_HORIZONTAL|wx.SL_TOP, name="Min Value")
        self.slider_m.Bind(wx.EVT_SCROLL_CHANGED, self.SliderMin)
        self.sizer.Add(self.slider_m, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.slider_max_default = 100
        self.slider_M = wx.Slider(self, -1, self.slider_max_default, 0, 100, 
                                wx.DefaultPosition,
                                wx.DefaultSize, wx.SL_LABELS|
                                wx.SL_HORIZONTAL|wx.SL_BOTTOM, name="Max Value")
        self.slider_M.Bind(wx.EVT_SCROLL_CHANGED, self.SliderMax)
        self.sizer.Add(self.slider_M, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------

        #self.sizer.Add(self.sizer2, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        self.SetSizer(self.sizer)
        self.Fit()

    def Load(self, event):
        self.mainparent.statusbar.SetStatusText("Loading Files ...", 0)

        self.dirname = defaults.dir
        dlg = wx.FileDialog(self, "Select File", self.dirname, 
                            "", "*", wx.FD_OPEN)
        if (dlg.ShowModal() == wx.ID_OK):
            self.mainparent.statusbar.SetStatusText("Please Wait ...", 0)
            self.mainparent.filename = dlg.GetFilename()
            self.mainparent.dirname = dlg.GetDirectory()
            self.mainparent.fullfilename = dlg.GetPath()
            self.mainparent.statusbar.SetStatusText("Current File: " + 
                                         self.mainparent.fullfilename, 2)
            # read file
            fil = str(self.mainparent.fullfilename)
            data, nrecs, theta, phi, radii, quantities, quant_names, ierr = \
               read_shell_slice(fil, rThetaPhiQuantities=True)

            if (ierr):
                self.mainparent.ThrowError("File Does Not Appear To Be Binary")

            else:

                nradii = len(radii)
                np = len(phi)
                nth = len(theta)
                nq = len(quant_names)

                self.mainparent.quantities = {}
                self.mainparent.quantity_arr = quantities
                self.mainparent.quantity_names = quant_names
                for i in range(nq):
                    self.mainparent.quantities[i] = quant_names[i]

                self.mainparent.nq = nq
                self.mainparent.nrecs = nrecs
                self.mainparent.nradii = nradii
                self.mainparent.np = np
                self.mainparent.nth = nth

                if (nrecs == 1):
                    tmp = numpy.empty((nth, np, nradii, nq, 2))
                    tmp[:,:,:,:,0] = data[:,:,:,:,0]
                    tmp[:,:,:,:,1] = data[:,:,:,:,0]
                    data = tmp
                    
                # size = (nth, np, nradii, nq, nrecs)
                self.mainparent.data = data

                self.mainparent.theta = theta
                self.mainparent.radii = radii
                self.mainparent.phi = phi

                pmin = numpy.amin(phi)
                pmax = numpy.amax(phi)
                thmin = numpy.amin(theta)
                thmax = numpy.amax(theta)
                self.mainparent.plot_xmin = pmin*180./numpy.pi
                self.mainparent.plot_xmax = pmax*180./numpy.pi
                # theta is standard spherical theta so when plotting,
                # lower theta is closer to top and larger theta is at bottom
                self.mainparent.plot_ymin = thmax*180./numpy.pi
                self.mainparent.plot_ymax = thmin*180./numpy.pi

                iq_init = 0
                nrad_init = nradii - 1
                nrec_init = 0
                self.mainparent.current_qnt = iq_init
                self.mainparent.current_rad = nrad_init
                self.mainparent.current_rec = nrec_init

                # update menus
                self.mainparent.UpdateMenu()

                # update status
                qnt = self.mainparent.quantities[iq_init]
                self.mainparent.statusbar.SetStatusText("Quantity: " + qnt, 1)
                self.mainparent.plot_title = "Shell Slice of " + qnt

                # set record slider
                if (nrecs == 1):
                    self.slider_r.SetMax(1)
                else:
                    self.slider_r.SetMax(nrecs-1)
                self.slider_r.SetValue(self.mainparent.current_rec)

                self.mainparent.statusbar.SetStatusText("Loaded: " + 
                                               self.mainparent.filename, 0)
                # update plot
                self.mainparent.plotpane.plot()

        dlg.Destroy()

    def Vorticity(self, event):
        if (self.mainparent.filename == ""):
            self.mainparent.ThrowError("Must Load a File First")
            return

        self.mainparent.ThrowWarning(\
                            "Vorticity Can Be Slow To Load....Please Wait", \
                            cancel=False)

        mainstatusbar = self.mainparent.statusbar
        mainstatusbar.SetStatusText("Calc Vorticity ...", 0)

        nth = self.mainparent.nth
        np  = self.mainparent.np
        nrad = self.mainparent.nradii
        irec = self.mainparent.current_rec
        iq = self.mainparent.current_qnt
        data = self.mainparent.data
        theta = self.mainparent.theta
        phi = self.mainparent.phi
        r = self.mainparent.radii

        iup = -1
        iut = -1
        for id, qnt in self.mainparent.quantities.iteritems():
            if (qnt == "Vtheta"):
                iut = id
            elif (qnt == "Vphi"):
                iup = id
        if (iut == -1):
            self.mainparent.ThrowError("Vtheta was not found")
            return
        if (iup == -1):
            self.mainparent.ThrowError("Vphi was not found")
            return

        wr = numpy.empty((nth, np, nrad))
        sin_t = numpy.sin(theta)
        up = data[:,:,:,iup,irec]
        ut = data[:,:,:,iut,irec]
        dpi = 1./(numpy.amax(phi) - numpy.amin(phi))
        dti = 1./(numpy.amax(theta) - numpy.amin(theta))

        for ir in range(nrad):

            for ip in range(np):
                tmp = sin_t*up[:,ip,ir]
                # derivative over theta, ibc=1 for dirichlet
                tmp2 = derivatives.compact_fd6(dti, tmp, 0., 0., 1, dtype=0)
                tmp2 /= r[ir]
                tmp2 /= sin_t

                wr[:,ip,ir] = tmp2[:]

            for it in range(nth):
                tmp = ut[it,:,ir]
                b1  = ut[it,np-4:np-1+1,ir]
                b2  = ut[it,0:3+1,ir]
                # derivative over phi, ibc=3 for periodic
                tmp2 = derivatives.compact_fd6(dpi, tmp, b1, b2, 3, dtype=0)
                tmp2 /= r[ir]*sin_t[it]

                wr[it,:,ir] = wr[it,:,ir] - tmp2[:]

        mainstatusbar.SetStatusText("Calc Vorticity ... Complete", 0)
        mainstatusbar.SetStatusText("Quantity: Vorticity", 1)

        self.mainparent.data_vort = wr
        self.mainparent.plot_vorticity = True
        self.mainparent.plot_enthalpy = False
        if (self.mainparent.tex):
            self.mainparent.plot_title = "Shell Slice of $\omega_r$"
        else:
            self.mainparent.plot_title = "Shell Slice of Vorticity"

        # plot it
        self.mainparent.plotpane.plot()

    def Compute_Enthalpy(self, event):
        if (self.mainparent.filename == ""):
            self.mainparent.ThrowError("Must Load a File First")
            return

        mainstatusbar = self.mainparent.statusbar
        mainstatusbar.SetStatusText("Calc Enthalpy ...", 0)

        data = self.mainparent.data
        r    = self.mainparent.radii
        irec = self.mainparent.current_rec
        nrad = self.mainparent.nradii
        nth  = self.mainparent.nth
        np   = self.mainparent.np
        cp   = constants.cp
        lsun = constants.lsol
        rsun = constants.rsol

        itemp = -1
        irho = -1
        iur = -1
        for id, qnt in self.mainparent.quantities.iteritems():
            if (qnt == "Temperature"):
                itemp = id
            elif (qnt == "Rho"):
                irho = id
            elif (qnt == "Vr"):
                iur = id

        if (itemp == -1):
            self.mainparent.ThrowError("Temperature was not found")
            return
        if (irho == -1):
            self.mainparent.ThrowError("Rho was not found")
            return
        if (iur == -1):
            self.mainparent.ThrowError("Vr was not found")
            return

        f_en = numpy.empty((nth, np, nrad))
        rho = data[:,:,:,irho,irec]
        temper = data[:,:,:,itemp,irec]
        ur = data[:,:,:,iur,irec]

        fpi = 4.*numpy.pi
        one_over_np_nth = 1./(np*nth)
        for ir in range(nrad):
            ruavg = one_over_np_nth*numpy.sum(rho[:,:,ir]*ur[:,:,ir])
            tavg  = one_over_np_nth*numpy.sum(temper[:,:,ir])

            rho_ur = rho[:,:,ir]*ur[:,:,ir] - ruavg
            temp_tmp = temper[:,:,ir] - tavg

            f_en[:,:,ir] = fpi*r[ir]*r[ir]*rho_ur[:,:]*cp*temp_tmp[:,:]/lsun

            print "Avg F_en at r(ir) = ",r[ir]/rsun, " is ", \
                  numpy.average(f_en[:,:,ir])

        mainstatusbar.SetStatusText("Calc Enthalpy ... Complete", 0)
        mainstatusbar.SetStatusText("Quantity: F_en", 1)

        self.mainparent.data_fen = f_en
        self.mainparent.plot_enthalpy = True
        self.mainparent.plot_vorticity = False
        self.mainparent.plot_title = "Shell Slice of F_en"

        # plot it
        self.mainparent.plotpane.plot()

    def Animate(self, event):
        if (self.mainparent.filename == ""):
            self.mainparent.ThrowError("Must Load a File First")
            return

        self.mainparent.statusbar.SetStatusText("Animating ...", 0)

        # set animation speed
        msg = r"Enter Animation Wait Time (sec)"
        dlg = wx.TextEntryDialog(self.parent, str(msg), "Enter a Value",
                                 defaultValue="0.05")
        if (dlg.ShowModal() == wx.ID_OK):
            result = str(dlg.GetValue())
            wait_time = float(result)
        else:
            wait_time = 0.05
        dlg.Destroy()

        question = "Do you want to save each animation frame?"
        dlg = wx.MessageDialog(self.parent, question, 
                               'Save Images for a Movie?', 
                               wx.YES_NO|wx.ICON_QUESTION)
        if (dlg.ShowModal() == wx.ID_YES):
            self.mainparent.plot_save_figure = True
            message = "Enter Basename for Each Image:"
            msg = wx.TextEntryDialog(self.parent, str(message), 
                                     "Enter a Basename",
                                     defaultValue="shell_slice_frame")
            if (msg.ShowModal() == wx.ID_OK):
                output_base = str(msg.GetValue())
            else:
                output_base = "shell_slice_frame"
            msg.Destroy()
        else:
            self.mainparent.plot_save_figure = False
            output_base = ""
        dlg.Destroy()

        nrad = self.mainparent.nradii

        for i in range(nrad):

            self.mainparent.current_rad = i

            if (defaults.eps):
                suffix = '.eps'
            else:
                suffix = '.png'
            self.mainparent.plot_save_output = output_base + str(i) + suffix

            self.mainparent.plotpane.plot()

            time.sleep(wait_time)

        # reset to displaying images after each frame was plotted
        self.mainparent.plot_save_figure = False

    def Quit(self, event):
        self.mainparent.statusbar.SetStatusText("Quitting ...", 0)
        self.mainparent.Close()

    def SliderRec(self, event):
        value = self.slider_r.GetValue()
        self.mainparent.current_rec = int(value)
        self.mainparent.statusbar.SetStatusText("New Rec Value: "+str(value),0)

        if (self.mainparent.plot_enthalpy):
            self.Compute_Enthalpy(event)
        elif (self.mainparent.plot_vorticity):
            self.Vorticity(event)
        else:
            self.mainparent.plotpane.plot()

    def SliderMin(self, event):
        value = self.slider_m.GetValue()
        self.mainparent.slider_min = float(value)
        self.mainparent.statusbar.SetStatusText("New Min Value: "+str(value),0)
        self.mainparent.plot_manual_slider = True
        self.mainparent.plotpane.plot()

    def SliderMax(self, event):
        value = self.slider_M.GetValue()
        self.mainparent.slider_max = float(value)
        self.mainparent.statusbar.SetStatusText("New Max Value: "+str(value),0)
        self.mainparent.plot_manual_slider = True
        self.mainparent.plotpane.plot()

    def Plot(self, event):
        self.mainparent.plotpane.plot()

    def Annotate_Plot(self, event):
        self.mainparent.plotpane.annotate_plot()

    def Change_Colormap(self, event):

        new_cmap = GetNewColorMap(self.parent, -1, "Choose a Colormap")
        new_cmap.ShowModal()
        if (self.mainparent.plot_rev_cmap):
            self.mainparent.plot_cmap += "_r"
        else:
            # reverse is False, so remove "_r" if it is there
            if (string.rfind(self.mainparent.plot_cmap, "_r") > 0):
                self.mainparent.plot_cmap = self.mainparent.plot_cmap[0:-2]

        c = self.mainparent.plot_cmap
        self.mainparent.statusbar.SetStatusText("New Colormap: "+str(c), 0)
        self.mainparent.plotpane.plot()

    def Clear_Plot(self, event):
        self.mainparent.plotpane.clear_plot()

#####################################################################
# Combo-Box class
#####################################################################
class GetNewColorMap(wx.Dialog):

    def __init__(self, parent, id, title):

        wx.Dialog.__init__(self, parent, id, title)

        self.parent = parent
        self.mainparent = self.GetParent().GetParent()

        self.choice = self.mainparent.plot_cmap

        MaxImgSize = min(self.mainparent.xs, self.mainparent.ys)

        img = "colormaps.png"
        if (not os.path.isfile(img)):
            self.mainparent.ThrowError("Image: "+img+" cannot be opened")
            return

        img = wx.Image(img, wx.BITMAP_TYPE_ANY)

        W = img.GetWidth()
        H = img.GetWidth()
        if W > H:
            newW = MaxImgSize
            newH = MaxImgSize * H / W
        else:
            newH = MaxImgSize
            newW = MaxImgSize * W / H
        img = img.Scale(newW, newH)

        image = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(newW, newH))
        image.SetBitmap(wx.BitmapFromImage(img))

        self.maps = ['Accent', 'Blues', 'BrBG', 'BuGn', 'BuPu', 'Dark2', 
                'GnBu', 'Greens', 'Greys', 'OrRd', 'Oranges', 'PRGn', 
                'Paired', 'Pastel1', 'Pastel2', 'PiYG', 'PuBu', 'PuBuGn', 
                'PuOr', 'PuRd', 'Purples', 'RdBu', 'RdGy', 'RdPu', 'RdYlBu', 
                'RdYlGn', 'Reds', 'Set1', 'Set2', 'Set3', 'Spectral', 'YlGn', 
                'YlGnBu', 'YlOrBr', 'YlOrRd', 'autumn', 'binary', 'bone', 
                'cool', 'copper', 'flag', 'gist_earth', 'gist_gray', 
                'gist_heat', 'gist_ncar', 'gist_rainbow', 'gist_stern', 
                'gist_yarg', 'gray', 'hot', 'hsv', 'jet', 'pink', 'prism', 
                'spectral', 'spring', 'summer', 'winter']

        self.cbox = wx.ComboBox(self, -1, choices=self.maps)
        self.cbox.SetValue(self.mainparent.plot_cmap)
        self.revbut = wx.CheckBox(self, -1, "Reverse Colorbar")
        self.revbut.SetValue(self.mainparent.plot_rev_cmap)
        okbut = wx.Button(self, wx.ID_OK, "OK")
        canbut = wx.Button(self, wx.ID_CANCEL, "Cancel")

        self.revbut.Bind(wx.EVT_CHECKBOX, self.On_Reverse)
        okbut.Bind(wx.EVT_BUTTON, self.On_OK)
        canbut.Bind(wx.EVT_BUTTON, self.On_Cancel)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        bsizer = wx.BoxSizer(wx.VERTICAL)

        bsizer.Add(self.cbox, 0, wx.ALL, border=5)
        bsizer.Add(self.revbut, 0, wx.ALL, border=5)
        bsizer.Add(okbut, 0, wx.ALL, border=5)
        bsizer.Add(canbut, 0, wx.ALL, border=5)

        sizer.Add(image, 0, wx.ALL, border=5)
        sizer.Add(bsizer, 0, wx.ALL, border=5)

        self.SetSizer(sizer)

        self.Fit()
        self.Centre()

    def On_Reverse(self, event):
        self.mainparent.plot_rev_cmap = self.revbut.GetValue()

    def On_OK(self, event):
        self.mainparent.plot_cmap = str(self.cbox.GetValue())
        self.Destroy()

    def On_Cancel(self, event):
        self.Destroy()

#####################################################################
# setup the plotting panel
#####################################################################
class plotpanel(wx.Panel):

    def __init__(self, parent, xs, ys):

        wx.Panel.__init__(self, parent)

        self.parent = parent

        self.mainparent = self.GetParent().GetParent()

        self.mydpi = defaults.dpi
        self.fact = 1.0 #0.95

        self.xs = xs
        self.ys = ys

        figx = self.fact*ys/self.mydpi
        figy = self.fact*ys/self.mydpi

        self.figure = plt.figure(figsize=(figx, figy), dpi=self.mydpi)

        self.mainparent.plot_xmin = None
        self.mainparent.plot_xmax = None
        self.mainparent.plot_ymin = None
        self.mainparent.plot_ymax = None

        self.canvas = FigureCanvas(self, -1, self.figure)
        self.toolbar = NavigationToolbar2Wx(self.canvas)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.EXPAND|wx.ALL|wx.GROW)
        self.sizer.Add(self.toolbar, 0, wx.EXPAND)
        self.SetSizer(self.sizer)
        self.Fit()

    def clear_plot(self):

        self.mainparent.statusbar.SetStatusText("Clearing Figure ...", 0)

        self.figure.clf()
        self.axes = self.figure.add_subplot(111)

        # display blank white square
        self.axes.get_xaxis().set_visible(False)
        self.axes.get_yaxis().set_visible(False)

        self.figure.canvas.draw()

    def plot(self):

        if (self.mainparent.filename != ""):

            self.mainparent.statusbar.SetStatusText("Plotting ...", 0)

            # shortcuts:
            xmin    = self.mainparent.plot_xmin
            ymin    = self.mainparent.plot_ymin
            xmax    = self.mainparent.plot_xmax
            ymax    = self.mainparent.plot_ymax
            title   = self.mainparent.plot_title
            xlabel  = self.mainparent.plot_xlabel
            ylabel  = self.mainparent.plot_ylabel
            cm      = matplotlib.cm.get_cmap(self.mainparent.plot_cmap)
            smin    = self.mainparent.slider_min
            smax    = self.mainparent.slider_max
            iq      = self.mainparent.current_qnt
            irad    = self.mainparent.current_rad
            irec    = self.mainparent.current_rec
            if (not self.mainparent.plot_enthalpy and \
                not self.mainparent.plot_vorticity):
                data  = self.mainparent.data[:,:,irad, iq, irec]
            elif (self.mainparent.plot_enthalpy):
                data = self.mainparent.data_fen[:,:,irad]
            elif (self.mainparent.plot_vorticity):
                data = self.mainparent.data_vort[:,:,irad]
            else:
                # should not ever be here
                self.mainparent.ThrowError("Bad Logic in Plot Routine")
                return
            theta   = self.mainparent.theta*180./numpy.pi
            radii   = self.mainparent.radii
            phi     = self.mainparent.phi
            cb_rng  = self.mainparent.plot_manual_colorbar_range
            sl_rng  = self.mainparent.plot_manual_slider
            m_vmin  = self.mainparent.plot_vmin
            m_vmax  = self.mainparent.plot_vmax
            cb_title= self.mainparent.plot_cmap_label

            self.figure.clf()
            self.axes = self.figure.add_subplot(111)
            self.axes.get_xaxis().set_visible(True)
            self.axes.get_yaxis().set_visible(True)

            if (self.mainparent.tex):
                # adjust the title and labels for Latex Compatibility
                for val in ["F_ks", "F_en", "F_ke", "F_ac", "F_sl"]:
                    ind  = string.find(title,  val)
                    if (ind > -1):
                        title = title[:ind+1] + "$_{" + \
                                title[ind+2:ind+len(val)+1] + "}$" + \
                                title[ind+len(val)+1:]
                        self.mainparent.plot_title = title
                        break

                for val in ["F_ks", "F_en", "F_ke", "F_ac", "F_sl"]:
                    ind  = string.find(xlabel,  val)
                    if (ind > -1):
                        xlabel = xlabel[:ind+1] + "$_{" + \
                                xlabel[ind+2:ind+len(val)+1] + "}$" + \
                                xlabel[ind+len(val)+1:]
                        self.mainparent.plot_xlabel = xlabel
                        break

                for val in ["F_ks", "F_en", "F_ke", "F_ac", "F_sl"]:
                    ind  = string.find(ylabel,  val)
                    if (ind > -1):
                        ylabel = ylabel[:ind+1] + "$_{" + \
                                ylabel[ind+2:ind+len(val)+1] + "}$" + \
                                ylabel[ind+len(val)+1:]
                        self.mainparent.plot_ylabel = ylabel
                        break

            maxdata = numpy.amax(data[:,:])
            mindata = numpy.amin(data[:,:])

            # the colorbar range was manually typed
            if (cb_rng):
                vmin = m_vmin
                vmax = m_vmax
                self.mainparent.plot_manual_colorbar_range = False

                smin = 100*(vmin - mindata)/(maxdata - mindata)
                smax = 100*(vmax - mindata)/(maxdata - mindata)

            # colorbar was set by slider
            elif (sl_rng):
                vmin = mindata + 0.01*smin*(maxdata - mindata)
                vmax = mindata + 0.01*smax*(maxdata - mindata)
                self.mainparent.plot_manual_slider = False

            # user has not set colorbar values, use defaults
            else:
                vmin = mindata + 0.01*smin*(maxdata - mindata)
                vmax = mindata + 0.01*smax*(maxdata - mindata)

            # check for invalid values
            if (vmin >= vmax):
                self.mainparent.ThrowError(\
                                    "Vmin =< Vmax, Setting Smin=25 & Smax=75")
                smin = 25
                smax = 75
                vmin = mindata + 0.01*smin*(maxdata - mindata)
                vmax = mindata + 0.01*smax*(maxdata - mindata)

            # store current value of vmin/vmax and slider positions
            self.mainparent.plot_vmin = vmin
            self.mainparent.plot_vmax = vmax
            self.mainparent.buttons.slider_m.SetValue(int(smin))
            self.mainparent.buttons.slider_M.SetValue(int(smax))
            self.mainparent.slider_min = float(smin)
            self.mainparent.slider_max = float(smax)

            extent = (xmin, xmax, ymin, ymax)
            cax = self.axes.imshow(data, interpolation='quadric', cmap=cm,
                                   extent=extent, origin='upper', 
                                   vmin=vmin, vmax=vmax, aspect='auto')

            cb = self.figure.colorbar(cax)
            cb.set_clim(vmin, vmax)
            cb.set_label(cb_title)

            # title, labels
            self.axes.set_title(title)
            self.axes.set_xlabel(xlabel)
            self.axes.set_ylabel(ylabel)

            # range
            self.axes.set_xlim(xmin, xmax)
            self.axes.set_ylim(ymin, ymax)

            # display proper figure
            if (not self.mainparent.plot_save_figure):
                self.figure.canvas.draw()
            else:
                self.figure.canvas.draw()
                if (defaults.eps):
                    fmt = 'eps'
                else:
                    fmt = 'png'
                plt.savefig(self.mainparent.plot_save_output, format=fmt)
                print "\tSaved frame: "+self.mainparent.plot_save_output

            self.mainparent.statusbar.SetStatusText("Plotting ... Complete", 0)

        else:
            self.mainparent.ThrowError("Must Load a File First")

    def annotate_plot(self):

        self.mainparent.statusbar.SetStatusText("Annotating ...", 0)

        #if (self.mainparent.tex):
        #    title = "Figure Attributes (TeX Enabled)"
        #else:
        #    title = "Figure Attributes (TeX Disabled)"
        title = "Figure Attributes (TeX Enabled)"
        mydlg = MyAttributeDialog(self.parent, -1, title)

        if (mydlg.ShowModal() == wx.ID_OK):

            # set values
            mydlg.SetValues()

            # update plot
            self.plot()

        mydlg.Destroy()

class MyAttributeDialog(wx.Dialog):

    def __init__(self, parent, id, title):

        wx.Dialog.__init__(self, parent, -1, title)

        self.mainparent = self.GetParent().GetParent()

        # define everything to put in the window
        title  = wx.StaticText(self, -1, "Title:")
        xlabel = wx.StaticText(self, -1, "X label:")
        ylabel = wx.StaticText(self, -1, "Y label:")
        xmin   = wx.StaticText(self, -1, "X min:")
        ymin   = wx.StaticText(self, -1, "Y min:")
        xmax   = wx.StaticText(self, -1, "X max:")
        ymax   = wx.StaticText(self, -1, "Y max:")
        cbmin  = wx.StaticText(self, -1, "Cbar min:")
        cbmax  = wx.StaticText(self, -1, "Cbar max:")
        cbtit  = wx.StaticText(self, -1, "Cbar Title:")

        t = self.mainparent.plot_title
        xl = self.mainparent.plot_xlabel
        yl = self.mainparent.plot_ylabel
        xm = str(self.mainparent.plot_xmin)
        ym = str(self.mainparent.plot_ymin)
        xM = str(self.mainparent.plot_xmax)
        yM = str(self.mainparent.plot_ymax)

        cbm = str(self.mainparent.plot_vmin)
        cbM = str(self.mainparent.plot_vmax)
        cbt = str(self.mainparent.plot_cmap_label)

        if (t == None):
            t = ""
        if (xl == None):
            xl = ""
        if (yl == None):
            yl = ""
        self.tc   = wx.TextCtrl(self, 1, t)
        self.xlc  = wx.TextCtrl(self, 2, xl)
        self.ylc  = wx.TextCtrl(self, 3, yl)
        self.xmc  = wx.TextCtrl(self, 4, xm)
        self.ymc  = wx.TextCtrl(self, 5, ym)
        self.xMc  = wx.TextCtrl(self, 6, xM)
        self.yMc  = wx.TextCtrl(self, 7, yM)
        self.cbmc = wx.TextCtrl(self, 8, cbm)
        self.cbMc = wx.TextCtrl(self, 9, cbM)
        self.cbtc = wx.TextCtrl(self, 10,cbt)

        ok = wx.Button(self, wx.ID_OK, 'OK')
        cancel = wx.Button(self, wx.ID_CANCEL, 'Cancel')

        # now do sizing
        hbox = wx.BoxSizer(wx.VERTICAL)
        titles = wx.BoxSizer(wx.HORIZONTAL)
        xlabels = wx.BoxSizer(wx.HORIZONTAL)
        ylabels = wx.BoxSizer(wx.HORIZONTAL)
        xminbox = wx.BoxSizer(wx.HORIZONTAL)
        yminbox = wx.BoxSizer(wx.HORIZONTAL)
        xmaxbox = wx.BoxSizer(wx.HORIZONTAL)
        ymaxbox = wx.BoxSizer(wx.HORIZONTAL)
        minbox = wx.BoxSizer(wx.VERTICAL)
        maxbox = wx.BoxSizer(wx.VERTICAL)
        range = wx.BoxSizer(wx.HORIZONTAL)
        crange = wx.BoxSizer(wx.HORIZONTAL)
        button = wx.BoxSizer(wx.HORIZONTAL)
        ctitles = wx.BoxSizer(wx.HORIZONTAL)
        cminbox = wx.BoxSizer(wx.HORIZONTAL)
        cmaxbox = wx.BoxSizer(wx.HORIZONTAL)

        titles.Add(title, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        titles.Add(self.tc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        xlabels.Add(xlabel, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        xlabels.Add(self.xlc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        ylabels.Add(ylabel, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        ylabels.Add(self.ylc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        xminbox.Add(xmin, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        xminbox.Add(self.xmc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        yminbox.Add(ymin, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        yminbox.Add(self.ymc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        xmaxbox.Add(xmax, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        xmaxbox.Add(self.xMc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        ymaxbox.Add(ymax, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        ymaxbox.Add(self.yMc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        cminbox.Add(cbmin, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        cminbox.Add(self.cbmc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        cmaxbox.Add(cbmax, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        cmaxbox.Add(self.cbMc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        ctitles.Add(cbtit, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        ctitles.Add(self.cbtc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)

        minbox.Add(xminbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        minbox.Add(yminbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        maxbox.Add(xmaxbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        maxbox.Add(ymaxbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)

        range.Add(minbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        range.Add(maxbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)

        crange.Add(cminbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        crange.Add(cmaxbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)

        button.Add(ok, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        button.Add(cancel, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)

        hbox.Add(titles, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        hbox.Add(xlabels, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        hbox.Add(ylabels, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        hbox.Add(range, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        hbox.Add(ctitles, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        hbox.Add(crange, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        hbox.Add(button, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)

        self.SetSizer(hbox)
        hbox.Fit(self)

        self.Center()

    def OnOK(self, event):

        self.Destroy()

    def OnCancel(self, event):

        self.Destroy()

    def SetValues(self):

        # adjust the proper value with the entered values
        # must convert to string first, GetValue() returns 'unicode'
        self.mainparent.plot_title = str(self.tc.GetValue())
        self.mainparent.plot_xlabel = str(self.xlc.GetValue())
        self.mainparent.plot_ylabel = str(self.ylc.GetValue())
        self.mainparent.plot_cmap_label = str(self.cbtc.GetValue())
        xmVal = str(self.xmc.GetValue())
        ymVal = str(self.ymc.GetValue())
        xMVal = str(self.xMc.GetValue())
        yMVal = str(self.yMc.GetValue())

        # did the vmin/vmax values change?
        if (self.cbmc.IsModified()):
            self.mainparent.plot_manual_colorbar_range = True
            self.mainparent.plot_vmin = float(str(self.cbmc.GetValue()))

        if (self.cbMc.IsModified()):
            self.mainparent.plot_manual_colorbar_range = True
            self.mainparent.plot_vmax = float(str(self.cbMc.GetValue()))

        if (self.mainparent.plot_title == "None"):
            self.mainparent.plot_title = ""

        if (self.mainparent.plot_xlabel == "None"):
            self.mainparent.plot_xlabel == ""

        if (self.mainparent.plot_ylabel == "None"):
            self.mainparent.plot_ylabel == ""

        if (xmVal == "None"):
            self.mainparent.plot_xmin = None
        else:
            self.mainparent.plot_xmin = float(xmVal)

        if (ymVal == "None"):
            self.mainparent.plot_ymin = None
        else:
            self.mainparent.plot_ymin = float(ymVal)

        if (xMVal == "None"):
            self.mainparent.plot_xmax = None
        else:
            self.mainparent.plot_xmax = float(xMVal)

        if (yMVal == "None"):
            self.mainparent.plot_ymax = None
        else:
            self.mainparent.plot_ymax = float(yMVal)

#####################################################################
# return a color and linestyle for many different plots
#####################################################################
def getLineColor(i):

    color = ['r', 'b',  'g',  'm', 'k', 'y', 'c']

    lines = ['-', '--', '-.', ':']

    c_ind = i % len(color)
    l_ind = i % len(lines)

    c = color[c_ind]
    l = lines[l_ind]

    return l, c

#####################################################################
#  Main program
#####################################################################
def main(tex=True):

    app = wx.App(redirect=False)
    frame = window(None, title="Plot Shell Slices", tex=tex)
    frame.Show()
    app.MainLoop()

    print "\n\n---Complete---\n"


#####################################################################
#  Define a usage program
#####################################################################
def usage():

    print "\nPython GUI to plot CSS shell slices data\n"
    print "Usage:\n"
    print "    ./plot_shell_slice.py [options]\n"
    print "    -h, --help             Display this help message\n"
    print "    --no-tex               Do not use TeX supported labels\n"


#####################################################################
#  If called as a command line program, this serves as the main prog
#####################################################################
if __name__ == "__main__":

    try:
        opts, args = getopt.getopt(sys.argv[1:], "h",
                          ["help", "no-tex"])

    except getopt.GetoptError:

        print "\n---ERROR: Unknown Command Line Option---\n"
        usage()
        sys.exit(2)

    # defaults
    tex = True

    for opt, arg in opts:

        if opt in ("-h", "--help"):
           usage()
           sys.exit(2)
        elif opt in ("--no-tex"):
           tex = False

    # call main program with appropriate arguments
    main(tex=tex)


