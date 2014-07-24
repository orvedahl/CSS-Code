#!/usr/bin/env python
#
# Plot CSS Azimuthal Averages as a GUI
#
#
# 2014-06-15 R. Orvedahl
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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx
#--------------------------------------------------------------------
from read_az_avg import *
import defaults
import constants
from az_averages import *
import streamfunction as sfunc
import curved_polar
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

    def Quantity(self, event):

        id = event.GetId()

        if (id not in self.quantities.keys()):
            self.ThrowError("Quantity Undefined")
        else:
            quant = self.quantities[id]
            self.statusbar.SetStatusText("Quantity: " + quant, 1)

            if (quant != "Omega"):
                self.quantitydata = self.data[:,:,id]

            else:
                self.quantitydata = self.omega[:,:]

            self.plot_1D_data = False
            self.loaded_psi = False
            self.plot_title = "Azimuthal Average of "+quant
            self.plot_cmap_label = ""
            self.plot_vmin = None
            self.plot_vmax = None
            self.slider_min = 0.0
            self.slider_max = 100.0
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

    def ThrowWarning(self, message):
        msg = wx.MessageDialog(None, message, 'WARNING', 
                               wx.OK|wx.ICON_EXCLAMATION)
        msg.ShowModal()

    def DefineAttributes(self):

        self.data = None           # hold data to be displayed (many-D)
        self.quantitydata = None   # hold specific quantity data (2D)

        self.quantity_arr = None   # array of quantity indices
        self.num_quants = 0        # number of quantities
        self.quantities = {}       # dictionary of ID:Name pairs
        self.quantity_names = None # names of quantities

        self.omega = None          # hold 2D omega data
        self.loaded_omega = False  # has omega been added to menu

        self.psi = None            # hold 2D streamfunction data
        self.loaded_psi = False    # do we plot psi

        self.xdata = None          # hold 1D data
        self.ydata = None          # 1D data (or 2D for multiple 1D plots)
        self.num_1d_data_sets = 0  # how many 1D data sets to plot
        self.labels = None         # labels for 1D data sets

        self.radius = None         # hold radius data
        self.theta = None          # hold theta data
        self.nr = 0                # number of elements in radius data
        self.nth = 0               # number of elements in theta data

        self.dirname = ""          # hold directory
        self.filename = ""         # hold filename
        self.fullfilename = ""     # hold full path to file

        self.slider_min = 0.0      # values of the sliders
        self.slider_max = 100.0    # values of the sliders

        # plot attributes
        self.plot_1D_data = False  # plot colorbar map or line plot
        self.plot_title = "Azimuthal Average"
        self.plot_xlabel = "Radius (cm)"
        self.plot_ylabel = "Theta (deg)"
        self.plot_xminc = None     # for the colorbar plot
        self.plot_xmaxc = None
        self.plot_yminc = None
        self.plot_ymaxc = None
        self.plot_xmino = None     # for the 1D line plot
        self.plot_xmaxo = None
        self.plot_ymino = None
        self.plot_ymaxo = None
        self.plot_cmap = defaults.cmap
        self.plot_rev_cmap = False
        self.plot_cmap_label = ""
        self.plot_vmin = None
        self.plot_vmax = None
        self.plot_manual_colorbar_range = False
        self.plot_manual_slider = False
        self.plot_aspect_ratio = 'auto'

    def SetMainMenuBar(self):

        # setup a menubar (hard coded IDs that are different make it easier
        #                  to handle each event)
        self.mainmenuBar = wx.MenuBar()

        fileMenu = wx.Menu()
        fileMenu.Append(-10, '&Load File(s)')
        fileMenu.Append(-20, '&Plot Data')
        fileMenu.Append(-30, '&Clear Plot')

        quit_item = wx.MenuItem(fileMenu, wx.ID_EXIT, '&Quit\tCtrl+w')
        fileMenu.AppendItem(quit_item)

        # empty on initialization
        quantity = wx.Menu()
        quantity.Append(-60, 'No File Loaded')

        # for a submenu
        #fileMenu.AppendMenu(wx.ID_ANY, '&Select Quantity', quantity)

        self.mainmenuBar.Append(fileMenu, '&File')
        self.mainmenuBar.Append(quantity, '&Select Quantity')
        self.SetMenuBar(self.mainmenuBar)

        # make each object do something
        self.Bind(wx.EVT_MENU, self.buttons.Load, id=-10)
        self.Bind(wx.EVT_MENU, self.buttons.Plot, id=-20)
        self.Bind(wx.EVT_MENU, self.buttons.Clear_Plot, id=-30)
        self.Bind(wx.EVT_MENU, self.OnQuit, quit_item)
        self.Bind(wx.EVT_MENU, self.Quantity, id=-60)

    def UpdateQuantitiesMenu(self):

        newquantity = wx.Menu()

        # add quantities and include action for each quantity
        for i in range(self.num_quants):
            newquantity.Append(i, self.quantities[i])
            self.Bind(wx.EVT_MENU, self.Quantity, id=i)

        # replace old "Select Quantity" menu with updated quantities
        self.mainmenuBar.Replace(1, newquantity, '&Select Quantity')

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
        self.Lfile = wx.Button(self, -1, "Load File(s)")
        self.Lfile.Bind(wx.EVT_BUTTON, self.Load)
        self.sizer.Add(self.Lfile, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.omegaBut = wx.Button(self, -1, "Plot Omega v R")
        self.omegaBut.Bind(wx.EVT_BUTTON, self.Omega)
        self.sizer.Add(self.omegaBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.windBut = wx.Button(self, -1, "Thermal Wind")
        self.windBut.Bind(wx.EVT_BUTTON, self.T_Wind)
        self.sizer.Add(self.windBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.streamBut = wx.Button(self, -1, "Streamlines")
        self.streamBut.Bind(wx.EVT_BUTTON, self.Streamlines)
        self.sizer.Add(self.streamBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.plotBut = wx.Button(self, -1, "Return to/Update Plot")
        self.plotBut.Bind(wx.EVT_BUTTON, self.Plot)
        self.sizer.Add(self.plotBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.clearBut = wx.Button(self, -1, "Clear Plot")
        self.clearBut.Bind(wx.EVT_BUTTON, self.Clear_Plot)
        self.sizer.Add(self.clearBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
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
        self.slider_min_default = 0
        self.slider = wx.Slider(self, -1, self.slider_min_default, 0, 100, 
                                wx.DefaultPosition,
                                wx.DefaultSize, wx.SL_LABELS|
                                wx.SL_HORIZONTAL|wx.SL_TOP, name="Min Value")
        self.slider.Bind(wx.EVT_SCROLL_CHANGED, self.SliderMin)
        self.sizer.Add(self.slider, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.slider_max_default = 100
        self.slider2 = wx.Slider(self, -1, self.slider_max_default, 0, 100, 
                                wx.DefaultPosition,
                                wx.DefaultSize, wx.SL_LABELS|
                                wx.SL_HORIZONTAL|wx.SL_BOTTOM, name="Max Value")
        self.slider2.Bind(wx.EVT_SCROLL_CHANGED, self.SliderMax)
        self.sizer.Add(self.slider2, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------

        #self.sizer.Add(self.sizer2, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        self.SetSizer(self.sizer)
        self.Fit()

    def Load(self, event):
        self.mainparent.statusbar.SetStatusText("Loading Files ...", 0)

        self.dirname = defaults.dir
        dlg = wx.FileDialog(self, "Select File(s)", self.dirname, 
                            "", "*", wx.FD_OPEN|wx.FD_MULTIPLE)
        if (dlg.ShowModal() == wx.ID_OK):
            self.mainparent.filename = dlg.GetFilename()
            files = []
            for p in dlg.GetPaths():
                files.append(str(p))
            self.mainparent.dirname = dlg.GetDirectory()
            self.mainparent.fullfilename = dlg.GetPath()
            self.mainparent.statusbar.SetStatusText("Loaded: " + 
                                           self.mainparent.filename, 0)
            if (len(files) == 1):
                self.mainparent.statusbar.SetStatusText("Current File: " + 
                                           files[0], 2)
            else:
                fnumbers = []
                for f in files:
                    ind = string.rfind(f, "_")
                    fnumbers.append(int(f[ind+1:]))
                # find lowest and highest
                lo = min(fnumbers)
                hi = max(fnumbers)
                ind = string.rfind(files[0], "_")
                allfiles = files[0][:ind+1]+str(lo)+"-"+str(hi)
                self.mainparent.statusbar.SetStatusText("Current Files: " + 
                                           allfiles, 2)

            # read 1st file
            data, nrecs, theta, radius, quantities, quant_names, ierr = \
               read_az_avg(files[0], rThetaQuantities=True)

            if (ierr):
                self.mainparent.ThrowError("File Does Not Appear To Be Binary")

            else:
                # average over files
                data = average_over_files(files, len(files), len(radius),
                                          len(theta), len(quant_names))

                self.mainparent.quantities = {}
                self.mainparent.quantity_arr = quantities
                self.mainparent.quantity_names = quant_names
                for i in range(len(quant_names)):
                    self.mainparent.quantities[i] = quant_names[i]

                self.mainparent.num_quants = len(quant_names)

                # if array has shape (x,y,z,1) convert it to (x,y,z)
                if (len(numpy.shape(data)) > 3):
                    print "---Reformating data array---"
                    print "\tOld shape:", numpy.shape(data)
                    data = data[:,:,:,0]
                    print "\tNew shape:",numpy.shape(data)
                self.mainparent.data = data
                iq_init = 0
                self.mainparent.quantitydata = data[:,:,iq_init]
                self.mainparent.theta = theta
                self.mainparent.radius = radius
                self.mainparent.nr = len(radius)
                self.mainparent.nth = len(theta)

                rmin = numpy.amin(radius)
                rmax = numpy.amax(radius)
                self.mainparent.plot_xminc = rmin
                self.mainparent.plot_xmaxc = rmax
                # theta is standard spherical theta so when plotting,
                # lower theta is closer to top and larger theta is at bottom
                # set limits in degrees (easier to interpret then radians)
                self.mainparent.plot_yminc = numpy.amax(theta)*180./numpy.pi
                self.mainparent.plot_ymaxc = numpy.amin(theta)*180./numpy.pi
                self.mainparent.plot_xmino = rmin
                self.mainparent.plot_xmaxo = rmax

                smin = self.mainparent.slider_min
                smax = self.mainparent.slider_max
                maxdata = numpy.amax(data)
                mindata = numpy.amin(data)
                delta = maxdata - mindata
                self.mainparent.plot_vmin = mindata + 0.01*smin*delta
                self.mainparent.plot_vmax = mindata + 0.01*smax*delta

                self.mainparent.loaded_omega = False
                self.mainparent.loaded_psi = False

                # update menu
                self.mainparent.UpdateQuantitiesMenu()

                # update status
                qnt = self.mainparent.quantities[iq_init]
                self.mainparent.statusbar.SetStatusText("Quantity: " + qnt, 1)
                self.mainparent.plot_title = "Azimuthal Average of " + qnt

                # update plot
                self.mainparent.plotpane.plot()

        dlg.Destroy()

    def Omega(self, event):
        if (self.mainparent.filename == ""):
            self.mainparent.ThrowError("Must Load a File First")
            return

        self.mainparent.statusbar.SetStatusText("Plotting Omega vs R ...", 0)
        old_quant = self.mainparent.statusbar.GetStatusText(1)
        self.mainparent.statusbar.SetStatusText("Previous "+old_quant, 1)

        msg = r"Enter $\Omega_0$ (rad/sec):"
        dlg = wx.TextEntryDialog(self.parent, str(msg), "Enter a Value",
                                 defaultValue="2.66e-6")
        if (dlg.ShowModal() == wx.ID_OK):
            result = str(dlg.GetValue())
            omega0 = float(result)
        else:
            omega0 = 2.66e-6
        dlg.Destroy()

        # find data-index associated with Vphi
        qphi = -2
        for id, qnt in self.mainparent.quantities.iteritems():
            if (qnt == "Vphi"):
                qphi = id
                break
        if (qphi == -2):
            self.mainparent.ThrowError("Vphi was not loaded")
            return

        # extract Vphi data
        vphi = self.mainparent.data[:,:,qphi]

        nr = self.mainparent.nr
        nth = self.mainparent.nth

        omega = numpy.empty((nth, nr))
        theta = self.mainparent.theta
        radius = self.mainparent.radius

        amp = 1.e9/(2.*numpy.pi) # converts rads to nHz
        for ir in range(nr):
            rinv = 1./radius[ir]
            for it in range(nth):
                sinthinv = 1./numpy.sin(theta[it])
                omega[it,ir] = amp*(omega0+vphi[it,ir]*rinv*sinthinv)

        self.mainparent.xdata = radius
        self.mainparent.num_1d_data_sets = 5
        ydata = numpy.empty((self.mainparent.num_1d_data_sets, nr))
        labels = []

        if (self.mainparent.tex):
            units = "$^\circ$ lat"
        else:
            units = " degrees lat"
        # fill the 1D data and labels
        ydata[0,:] = omega[nth/2,:]
        l = round(90. - theta[nth/2]*180./numpy.pi)
        labels.append(str(l)+units)

        ydata[1,:] = 0.5*(omega[nth/3,:]+omega[2*nth/3,:])
        l = round(90. - theta[nth/3]*180./numpy.pi)
        labels.append(str(l)+units)

        ydata[2,:] = 0.5*(omega[nth/4,:]+omega[3*nth/4,:])
        l = round(90. - theta[nth/4]*180./numpy.pi)
        labels.append(str(l)+units)

        ydata[3,:] = 0.5*(omega[nth/6,:]+omega[5*nth/6,:])
        l = round(90. - theta[nth/6]*180./numpy.pi)
        labels.append(str(l)+units)

        # Omega_0 in nHz
        ydata[4,:] = omega0*numpy.ones((nr))*1.e9/2./numpy.pi
        if (self.mainparent.tex):
            labels.append(str("$\Omega_0$"))
        else:
            labels.append(str("Omega0"))

        self.mainparent.omega = omega
        self.mainparent.ydata = ydata
        self.mainparent.labels = labels

        # update menu if omega has not been loaded already for this file
        if (not self.mainparent.loaded_omega):
            self.mainparent.quantities[self.mainparent.num_quants] = "Omega"
            self.mainparent.quantity_names += ["Omega"]
            self.mainparent.num_quants += 1
            self.mainparent.UpdateQuantitiesMenu()

        self.mainparent.plot_1D_data = True
        self.mainparent.loaded_omega = True
        self.mainparent.loaded_psi = False

        # plot it
        self.mainparent.plotpane.plot()

    def T_Wind(self, event):
        if (self.mainparent.filename == ""):
            self.mainparent.ThrowError("Must Load a File First")
            return

        self.mainparent.ThrowError("Thermal Wind Not Yet Supported")

        #self.mainparent.statusbar.SetStatusText("Plot Thermal Wind ...", 0)

    def Streamlines(self, event):
        if (self.mainparent.filename == ""):
            self.mainparent.ThrowError("Must Load a File First")
            return

        self.mainparent.statusbar.SetStatusText("Plot Streamlines ...", 0)

        nth = self.mainparent.nth
        nr = self.mainparent.nr
        radius = self.mainparent.radius
        theta = self.mainparent.theta
        data = self.mainparent.data

        iqrho = -1
        iqvr = -1
        iqvt = -1
        for id, qnt in self.mainparent.quantities.iteritems():
            if (qnt == "Rho"):
                iqrho = id
            elif (qnt == "Vr"):
                iqvr = id
            elif (qnt == "Vtheta"):
                iqvt = id
        if (iqrho == -1):
            self.mainparent.ThrowError("Rho was not found")
            return
        if (iqvr == -1):
            self.mainparent.ThrowError("Vr was not found")
            return
        if (iqvt == -1):
            self.mainparent.ThrowError("Vtheta was not found")
            return

        #rhovr = numpy.empty((nth, nr))
        #rhovr[:,:] = data[:,:,iqrho]*data[:,:,iqvr]
        vr = numpy.empty((nth, nr))
        vr[:,:] = data[:,:,iqvr]

        #rhovtheta = numpy.empty((nth, nr))
        #rhovtheta[:,:] = data[:,:,iqrho]*data[:,:,iqvt]
        vtheta = numpy.empty((nth, nr))
        vtheta[:,:] = data[:,:,iqvt]

        psi = numpy.empty((nth, nr))
        #psi = sfunc.streamfunction(rhovr, rhovtheta, radius, numpy.cos(theta))
        psi = sfunc.streamfunction(vr, vtheta, radius, numpy.cos(theta))

        self.mainparent.psi = psi
        self.mainparent.loaded_psi = True
        self.mainparent.plot_1D_data = False

        self.mainparent.plot_title = "Stokes Streamfunction"
        self.mainparent.statusbar.SetStatusText("Quantity: Streamlines", 1)

        self.mainparent.plotpane.plot()

    def Quit(self, event):
        self.mainparent.statusbar.SetStatusText("Quitting ...", 0)
        self.mainparent.Close()

    def SliderMin(self, event):
        value = self.slider.GetValue()
        self.mainparent.slider_min = float(value)
        self.mainparent.statusbar.SetStatusText("New Min Value: "+str(value), 0)
        self.mainparent.plot_manual_slider = True
        self.mainparent.plotpane.plot()

    def SliderMax(self, event):
        value = self.slider2.GetValue()
        self.mainparent.slider_max = float(value)
        self.mainparent.statusbar.SetStatusText("New Max Value: "+str(value), 0)
        self.mainparent.plot_manual_slider = True
        self.mainparent.plotpane.plot()

    def Plot(self, event):
        self.mainparent.plot_1D_data = False
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

        self.mydpi = 100
        self.fact = 1.0 #0.95

        self.xs = xs
        self.ys = ys

        figx = self.fact*ys/self.mydpi
        figy = self.fact*ys/self.mydpi

        self.figure = plt.figure(figsize=(figx, figy), dpi=self.mydpi)

        self.mainparent.plot_xminc = None # colorbar plot
        self.mainparent.plot_xmaxc = None
        self.mainparent.plot_yminc = None
        self.mainparent.plot_ymaxc = None
        self.mainparent.plot_xmino = None # omega plot
        self.mainparent.plot_xmaxo = None
        self.mainparent.plot_ymino = None
        self.mainparent.plot_ymaxo = None

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
            plot_1D  = self.mainparent.plot_1D_data
            xminc    = self.mainparent.plot_xminc
            yminc    = self.mainparent.plot_yminc
            xmaxc    = self.mainparent.plot_xmaxc
            ymaxc    = self.mainparent.plot_ymaxc
            xmino    = self.mainparent.plot_xmino
            ymino    = self.mainparent.plot_ymino
            xmaxo    = self.mainparent.plot_xmaxo
            ymaxo    = self.mainparent.plot_ymaxo
            title    = self.mainparent.plot_title
            xlabel   = self.mainparent.plot_xlabel
            ylabel   = self.mainparent.plot_ylabel
            cm       = matplotlib.cm.get_cmap(self.mainparent.plot_cmap)
            smin     = self.mainparent.slider_min
            smax     = self.mainparent.slider_max
            plot_psi = self.mainparent.loaded_psi
            if (not plot_psi):
                data = self.mainparent.quantitydata
            else:
                data = self.mainparent.psi
            theta    = self.mainparent.theta*180./numpy.pi
            radius   = self.mainparent.radius
            xdata    = self.mainparent.xdata
            ydata    = self.mainparent.ydata
            n1Dsets  = self.mainparent.num_1d_data_sets
            labels   = self.mainparent.labels
            cb_rng   = self.mainparent.plot_manual_colorbar_range
            sl_rng   = self.mainparent.plot_manual_slider
            m_vmin   = self.mainparent.plot_vmin
            m_vmax   = self.mainparent.plot_vmax
            cb_title = self.mainparent.plot_cmap_label
            aspect   = self.mainparent.plot_aspect_ratio

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

            maxdata = numpy.amax(data)
            mindata = numpy.amin(data)

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
            self.mainparent.buttons.slider.SetValue(int(smin))
            self.mainparent.buttons.slider2.SetValue(int(smax))
            self.mainparent.slider_min = float(smin)
            self.mainparent.slider_max = float(smax)

            if (not plot_1D):
                extent = (xminc, xmaxc, yminc, ymaxc)

                # FIXME: implement working version of curved plots
                if (True):
                    cax = self.axes.imshow(data, interpolation='quadric', 
                                           cmap=cm, extent=extent, 
                                           origin='upper', vmin=vmin, 
                                           vmax=vmax, aspect='auto')

                    cb = self.figure.colorbar(cax)
                    cb.set_clim(vmin, vmax)
                    cb.set_label(cb_title)

                # EXPERIMENTAL: curved plots
                else:
                    data2 = numpy.zeros((len(theta), len(radius)))
                    for i in range(len(theta)):
                        th=theta[i]
                        for j in range(len(radius)):
                            rad=radius[j]
                            data2[i,j] = rad*numpy.sin(th)
                    curved_polar.polar_image(data2, radius, theta, 
                                       radians=False, vmin=vmin, vmax=vmax,
                                       cmap=cm, aspect=aspect, add_cont=False,
                                       extent=None, cb_title=cb_title, 
                                       r_bcz=[0.95*constants.rsol])

                # title, labels
                self.axes.set_title(title)
                self.axes.set_xlabel(xlabel)
                self.axes.set_ylabel(ylabel)

                # range
                self.axes.set_xlim(xminc, xmaxc)
                self.axes.set_ylim(yminc, ymaxc)

            else:

                for i in range(n1Dsets):
                    l, c = getLineColor(i)
                    self.axes.plot(xdata, ydata[i,:], label=labels[i], 
                                   color=c,linestyle=l)

                self.axes.legend(loc="upper left")
                leg = self.figure.gca().get_legend()
                #plt.setp(leg.get_texts(), fontsize='small')
                leg.draw_frame(False)

                # title, labels
                self.axes.set_title("Rotation Profiles")
                self.axes.set_xlabel("Radius (cm)")
                if (self.mainparent.tex):
                    self.axes.set_ylabel(r"$\Omega$ (nHz)")
                else:
                    self.axes.set_ylabel(r"Omega (nHz)")

                # range
                self.axes.set_xlim(xmino, xmaxo)
                self.axes.set_ylim(ymino, ymaxo)


            # display proper figure
            if (not plot_1D):
                self.figure.canvas.draw()
            else:
                self.figure.canvas.draw()

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
        aspect = wx.StaticText(self, -1, "Aspect:")
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
        asp = str(self.mainparent.plot_aspect_ratio)
        if (self.mainparent.plot_1D_data):
            xm = str(self.mainparent.plot_xmino)
            ym = str(self.mainparent.plot_ymino)
            xM = str(self.mainparent.plot_xmaxo)
            yM = str(self.mainparent.plot_ymaxo)
        else:
            xm = str(self.mainparent.plot_xminc)
            ym = str(self.mainparent.plot_yminc)
            xM = str(self.mainparent.plot_xmaxc)
            yM = str(self.mainparent.plot_ymaxc)

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
        self.aspc = wx.TextCtrl(self, 4, asp)
        self.xmc  = wx.TextCtrl(self, 5, xm)
        self.ymc  = wx.TextCtrl(self, 6, ym)
        self.xMc  = wx.TextCtrl(self, 7, xM)
        self.yMc  = wx.TextCtrl(self, 8, yM)
        self.cbmc = wx.TextCtrl(self, 9, cbm)
        self.cbMc = wx.TextCtrl(self, 10,cbM)
        self.cbtc = wx.TextCtrl(self, 11,cbt)

        ok = wx.Button(self, wx.ID_OK, 'OK')
        cancel = wx.Button(self, wx.ID_CANCEL, 'Cancel')

        # now do sizing
        hbox = wx.BoxSizer(wx.VERTICAL)
        titles = wx.BoxSizer(wx.HORIZONTAL)
        xlabels = wx.BoxSizer(wx.HORIZONTAL)
        ylabels = wx.BoxSizer(wx.HORIZONTAL)
        aspects = wx.BoxSizer(wx.HORIZONTAL)
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
        aspects.Add(aspect, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        aspects.Add(self.aspc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
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
        hbox.Add(aspects, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
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
        aspect_val = str(self.aspc.GetValue())

        # set aspect ratio
        try:
            aspc = float(aspect_val)
        except ValueError:
            aspc = str(aspect_val)

        self.mainparent.plot_aspect_ratio = aspc

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

        if (self.mainparent.plot_1D_data):
            if (xmVal == "None"):
                self.mainparent.plot_xmino = None
            else:
                self.mainparent.plot_xmino = float(xmVal)

            if (ymVal == "None"):
                self.mainparent.plot_ymino = None
            else:
                self.mainparent.plot_ymino = float(ymVal)

            if (xMVal == "None"):
                self.mainparent.plot_xmaxo = None
            else:
                self.mainparent.plot_xmaxo = float(xMVal)

            if (yMVal == "None"):
                self.mainparent.plot_ymaxo = None
            else:
                self.mainparent.plot_ymaxo = float(yMVal)

        else:
            if (xmVal == "None"):
                self.mainparent.plot_xminc = None
            else:
                self.mainparent.plot_xminc = float(xmVal)

            if (ymVal == "None"):
                self.mainparent.plot_yminc = None
            else:
                self.mainparent.plot_yminc = float(ymVal)

            if (xMVal == "None"):
                self.mainparent.plot_xmaxc = None
            else:
                self.mainparent.plot_xmaxc = float(xMVal)

            if (yMVal == "None"):
                self.mainparent.plot_ymaxc = None
            else:
                self.mainparent.plot_ymaxc = float(yMVal)


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
    frame = window(None, title="Plot Azimuthal Averages", tex=tex)
    frame.Show()
    app.MainLoop()

    print "\n\n---Complete---\n"


#####################################################################
#  Define a usage program
#####################################################################
def usage():

    print "\nPython GUI to plot CSS azimuthal average data\n"
    print "Usage:\n"
    print "    ./plot_az_avgs.py [options]\n"
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


