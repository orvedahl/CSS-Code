#!/usr/bin/env python
#
# Plot CSS Shell Averages as a GUI
#
#
# 2014-06-30 R. Orvedahl
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
from read_shell_avg import *
import defaults
import constants
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

        # define a few attributes
        self.DefineAttributes()

        # setup Main window with a plot pane and a buttons pane
        self.mainpane = wx.Panel(self)
        self.plotpane = plotpanel(self.mainpane, self.xs, self.ys)
        self.buttons = buttonpanel(self.mainpane)

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

    #-----------------------------------------------
    # handle a "File-->Select Quantity" Event
    #-----------------------------------------------
    def Quantity(self, event):

        id = event.GetId()

        if (id not in self.quantities.keys()):
            self.ThrowError("Quantity Undefined")
        else:
            quant = self.quantities[id]
            self.statusbar.SetStatusText("Quantity: " + quant, 1)

            if (quant != "Omega"):
                self.plot_ylabel = quant

            else:
                if (self.tex):
                    self.plot_ylabel = "$\Omega / \Omega_0$"
                else:
                    self.plot_ylabel = "Omega / Omega0"

            self.current_qnt = id
            self.plot_title = "Shell Average of "+quant

            # reset y-axis
            self.plot_ymin = None
            self.plot_ymax = None

            self.plotpane.plot()

    #-----------------------------------------------
    # handle a "Change Plot Attributes" Event
    #-----------------------------------------------
    def Annotate(self, event):

        self.plotpane.annotate_plot()

    #-----------------------------------------------
    # handle a Quit Event
    #-----------------------------------------------
    def OnQuit(self, event):

        self.plotpane.Destroy()
        self.buttons.Destroy()
        self.mainpane.Destroy()
        self.Close()

    #-----------------------------------------------
    # display an error message
    #-----------------------------------------------
    def ThrowError(self, message):

        msg = wx.MessageDialog(None, message, 'ERROR', wx.OK|wx.ICON_ERROR)
        msg.ShowModal()

    #-----------------------------------------------
    # display a warning message
    #-----------------------------------------------
    def ThrowWarning(self, message):

        msg = wx.MessageDialog(None, message, 'WARNING', 
                               wx.OK|wx.ICON_EXCLAMATION)
        msg.ShowModal()

    #-----------------------------------------------
    # define various variables associated with the class
    #-----------------------------------------------
    def DefineAttributes(self):

        self.quantity_arr = None   # array of quantity indices
        self.nq = 0                # number of quantities
        self.nq_menu = 0           # number of quantities in menu
        self.quantities = {}       # dictionary of ID:Name pairs
        self.quantity_names = None # names of quantities
        self.current_qnt = 0       # current quantity
        self.data = None           # hold the 2D data
        self.data_avg = None       # hold average of data

        self.domega = None         # hold 2D omega data
        self.loaded_omega = False  # has omega been added to menu

        self.radius = None         # hold radius data
        self.nr = 0                # number of elements in radius data
        self.nrecs = 0             # number of records in data file

        self.dirname = ""          # hold directory
        self.filename = ""         # hold filename
        self.fullfilename = ""     # hold full path to file

        # plot attributes
        self.plot_title = "Shell Average"
        if (self.tex):
            self.plot_xlabel = "r/R$_\odot$"
        else:
            self.plot_xlabel = "r/Rsol"
        self.plot_ylabel = ""
        self.plot_xmin = None
        self.plot_xmax = None
        self.plot_ymin = None
        self.plot_ymax = None
        self.plot_rsol = constants.rsol
        self.show_legend = False
        self.leg_loc = "lower left"
        self.plot_avg = False
        self.plot_animation = False

        self.log10 = False
        self.perturb = False
        self.luminosity = False

    #-----------------------------------------------
    # setup the main "File" menu
    #-----------------------------------------------
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

    #-----------------------------------------------
    # update the "Select Quantity" Menu with updated quantities
    #-----------------------------------------------
    def UpdateQuantitiesMenu(self):

        newquantity = wx.Menu()

        # add quantities and include action for each quantity
        for i in range(self.nq_menu):
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
        self.animBut = wx.Button(self, -1, "Animate")
        self.animBut.Bind(wx.EVT_BUTTON, self.Animate)
        self.sizer.Add(self.animBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.pertBut = wx.ToggleButton(self, -1, "Perturb")
        self.pertBut.Bind(wx.EVT_TOGGLEBUTTON, self.Perturb)
        self.sizer.Add(self.pertBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        self.pertBut.SetValue(self.mainparent.perturb)
        #--------------------------------------------------------------
        self.lumBut = wx.ToggleButton(self, -1, "Luminosity")
        self.lumBut.Bind(wx.EVT_TOGGLEBUTTON, self.Luminosity)
        self.sizer.Add(self.lumBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        self.lumBut.SetValue(self.mainparent.luminosity)
        #--------------------------------------------------------------
        self.omegaBut = wx.Button(self, -1, "Omega")
        self.omegaBut.Bind(wx.EVT_BUTTON, self.Omega)
        self.sizer.Add(self.omegaBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.log10But = wx.ToggleButton(self, -1, "Log 10")
        self.log10But.Bind(wx.EVT_TOGGLEBUTTON, self.Log10)
        self.sizer.Add(self.log10But, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        self.log10But.SetValue(self.mainparent.log10)
        #--------------------------------------------------------------
        self.avgBut = wx.Button(self, -1, "Average")
        self.avgBut.Bind(wx.EVT_BUTTON, self.Average)
        self.sizer.Add(self.avgBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.plotBut = wx.Button(self, -1, "Return to/Update Plot")
        self.plotBut.Bind(wx.EVT_BUTTON, self.Plot)
        self.sizer.Add(self.plotBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.annotateBut = wx.Button(self, -1, "Change Plot Attributes")
        self.annotateBut.Bind(wx.EVT_BUTTON, self.Annotate_Plot)
        self.sizer.Add(self.annotateBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------
        self.quitBut = wx.Button(self, -1, "Quit")
        self.quitBut.Bind(wx.EVT_BUTTON, self.Quit)
        self.sizer.Add(self.quitBut, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        #--------------------------------------------------------------

        self.SetSizer(self.sizer)
        self.Fit()

    #-----------------------------------------------
    # handle "Load File" Event
    #-----------------------------------------------
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
            values, nrecs, radius, quantities, quant_names, ierr = \
                                read_shell_avg(files[0], rQuantities=True)

            if (ierr):
                self.mainparent.ThrowError(\
                         "File Does Not Appear To Be Binary\n"+files[0])

            else:
                self.mainparent.quantities = {}
                self.mainparent.quantity_arr = quantities
                self.mainparent.quantity_names = quant_names
                for i in range(len(quant_names)):
                    self.mainparent.quantities[i] = quant_names[i]

                nr = len(radius)
                nq = len(quant_names)

                self.mainparent.nq = nq
                self.mainparent.nq_menu = nq
                self.mainparent.nr = nr

                # size is (nr x nq x nfiles*nrecs)
                data = numpy.empty((nr, nq, len(files)*nrecs))
                data[:, 0:nq-1+1, 0:nrecs-1+1] = values

                if (len(files) > 1):
                    for i in range(1, len(files)):
                        values, nrs, ierr = read_shell_avg(files[i])
                        if (ierr):
                            self.mainparent.ThrowError(\
                              "File Does Not Appear To Be Binary\n"+files[i])
                        else:
                            data[:,0:nq-1+1,nrecs:nrecs+nrs-1+1] = values
                            nrecs += nrs

                self.mainparent.nrecs = nrecs

                self.mainparent.radius = radius
                rsol = self.mainparent.plot_rsol
                self.mainparent.data = data
                self.mainparent.plot_xmin = numpy.amin(radius/rsol)
                self.mainparent.plot_xmax = numpy.amax(radius/rsol)
                self.mainparent.plot_ymin = None
                self.mainparent.plot_ymax = None

                self.mainparent.loaded_omega = False

                # update menu
                self.mainparent.UpdateQuantitiesMenu()

                # update status
                iq_init = 0
                self.mainparent.current_qnt = iq_init
                qnt = self.mainparent.quantities[iq_init]

                self.mainparent.statusbar.SetStatusText("Quantity: " + qnt, 1)
                self.mainparent.plot_title = "Shell Average of " + qnt
                self.mainparent.plot_ylabel = qnt

                # update plot
                self.mainparent.plotpane.plot()

        dlg.Destroy()

    #-----------------------------------------------
    # handle "Perturb" Event
    #-----------------------------------------------
    def Perturb(self, event):

        perturb = self.pertBut.GetValue()

        self.mainparent.perturb = perturb
        self.mainparent.statusbar.SetStatusText("Perturb: "+str(perturb), 0)

        # update plot
        self.mainparent.plotpane.plot()

    #-----------------------------------------------
    # handle "Luminosity" Event
    #-----------------------------------------------
    def Luminosity(self, event):

        lum = self.lumBut.GetValue()

        self.mainparent.luminosity = lum
        self.mainparent.statusbar.SetStatusText("Luminosity: "+str(lum), 0)

        # update plot
        self.mainparent.plotpane.plot()

    #-----------------------------------------------
    # handle "Log10" Event
    #-----------------------------------------------
    def Log10(self, event):

        log10 = self.log10But.GetValue()

        self.mainparent.log10 = log10
        self.mainparent.statusbar.SetStatusText("Log 10: "+str(log10), 0)

        # update plot
        self.mainparent.plotpane.plot()

    #-----------------------------------------------
    # handle "Average" Event
    #-----------------------------------------------
    def Average(self, event):

        self.mainparent.statusbar.SetStatusText("Average ...", 0)

        if (self.mainparent.filename == ""):
            self.mainparent.ThrowError("Must Load a File First")
            return

        iq = self.mainparent.current_qnt
        nq = self.mainparent.nq
        nrecs = self.mainparent.nrecs
        data = self.mainparent.data
        domega = self.mainparent.domega

        if (iq == nq and domega == None):
            self.mainparent.ThrowError(\
                            "Must load Omega before trying to average it")
            return

        # Omega
        if (iq == nq):
            avg = numpy.sum(domega, axis=1)/nrecs
        else:
            avg = numpy.sum(data[:,iq,:], axis=1)/nrecs

        self.mainparent.data_avg = avg
        self.mainparent.plot_avg = True

        avg_title = " (avg over records)"
        plot_title = self.mainparent.plot_title
        if (string.find(plot_title, avg_title) < 0):
            self.mainparent.plot_title = plot_title + avg_title

        # update plot
        self.mainparent.plotpane.plot()

    #-----------------------------------------------
    # handle "Animate" Event
    #-----------------------------------------------
    def Animate(self, event):
        if (self.mainparent.filename == ""):
            self.mainparent.ThrowError("Must Load a File First")
            return

        self.mainparent.statusbar.SetStatusText("Animate ...", 0)

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

        iq      = self.mainparent.current_qnt
        nq      = self.mainparent.nq
        iqr     = self.mainparent.quantity_arr[iq]
        iql     = self.mainparent.quantity_arr[nq-1]
        nr      = self.mainparent.nr
        radius  = self.mainparent.radius
        rsol    = self.mainparent.plot_rsol
        if (iq != nq):
            data = self.mainparent.data[:,iq,:]
        else:
            data = self.mainparent.domega
        four_pi_r2 = 4.*numpy.pi*radius*radius
        nrecs   = self.mainparent.nrecs
        lum     = self.mainparent.luminosity
        perturb = self.mainparent.perturb

        if (iqr <= iql):
            if ((lum) and ((iqr >= 15) and (iqr <= 19))):
                lum = numpy.empty((nrecs, nr))
                lum[:,:] = 0.0
                minv = 0.0
                maxv = 0.0
                for ir in range(nrecs):
                    tmax = numpy.amax(lum[ir,:])
                    tmin = numpy.amin(lum[ir,:])
                    lum[ir,:] = four_pi_r2*data[:,ir]
                    if (tmax > maxv):
                        maxv = tmax
                    if (tmin < minv):
                        minv = tmin
                for i in range(nrecs):
                    self.mainparent.plot_animation = True
                    self.mainparent.plotpane.plot(iter=i, x=radius/rsol, \
                                           y=lum[i,:], alabel="rec: "+str(i))
                    time.sleep(wait_time)

            elif (perturb):
                minv = 0.0
                maxv = 0.0
                for ir in range(nrecs):
                    tmax = numpy.amax(data[:,ir] - data[:,0])
                    tmin = numpy.amin(data[:,ir] - data[:,0])
                    if (tmax > maxv):
                        maxv = tmax
                    if (tmin < minv):
                        minv = tmin
                for ir in range(nrecs):
                    self.mainparent.plot_animation = True
                    self.mainparent.plotpane.plot(iter=ir, x=radius/rsol, \
                                                  y=data[:,ir]-data[:,0], \
                                                  alabel="rec: "+str(ir)+"-0")
                    time.sleep(wait_time)

            else:
                for ir in range(nrecs):
                    self.mainparent.plot_animation = True
                    self.mainparent.plotpane.plot(iter=ir, x=radius/rsol, \
                                         y=data[:,ir], alabel="rec: "+str(ir))
                    time.sleep(wait_time)

        else:

            temp = self.mainparent.domega.copy()
            if (temp == None):
                self.mainparent.ThrowError("Must Load Omega Before Animating")
                return

            if (perturb):
                for i in range(1,nrecs):
                    temp[:,i] = temp[:,i] - temp[:,0]
                temp[:,0] = 0.0

            for ir in range(nrecs):
                self.mainparent.plot_animation = True
                if (perturb):
                    lab = "rec: "+str(ir)+"-0"
                else:
                    lab = "rec: "+str(ir)
                self.mainparent.plotpane.plot(iter=ir, x=radius/rsol, \
                                              y=temp[:,ir], alabel=lab)
                time.sleep(wait_time)


    #-----------------------------------------------
    # handle "Omega" Event
    #-----------------------------------------------
    def Omega(self, event):
        if (self.mainparent.filename == ""):
            self.mainparent.ThrowError("Must Load a File First")
            return

        self.mainparent.statusbar.SetStatusText("Plotting Omega ...", 0)

        msg = r"Enter $\Omega_0$ (rad/sec)"
        dlg = wx.TextEntryDialog(self.parent, str(msg), "Enter a Value",
                                 defaultValue="2.66e-6")
        if (dlg.ShowModal() == wx.ID_OK):
            result = str(dlg.GetValue())
            omega0 = float(result)
        else:
            omega0 = 2.66e-6
        dlg.Destroy()

        nr = self.mainparent.nr
        nrecs = self.mainparent.nrecs

        domega = numpy.empty((nr, nrecs))
        radius = self.mainparent.radius

        # find Vphi
        qphi = -2
        for id, qnt in self.mainparent.quantities.iteritems():
            if (qnt == "Vphi"):
                qphi = id
                break
        if (qphi == -2):
            self.mainparent.ThrowError("Vphi was not loaded")
            return

        for ir in range(nrecs):
            domega[:,ir] = omega0 + self.mainparent.data[:,qphi,ir]/radius

        domega[:,:] = domega[:,:]/omega0

        self.mainparent.domega = domega

        # update menu if omega has not been loaded already for this file
        if (not self.mainparent.loaded_omega):
            self.mainparent.quantities[self.mainparent.nq] = "Omega"
            self.mainparent.quantity_names += ["Omega"]
            self.mainparent.nq_menu += 1
            q = self.mainparent.quantity_arr
            # add Omega as the last quantity in the array
            q = list(q)
            q = q + [q[-1]+1]
            q = numpy.array(q)
            self.mainparent.quantity_arr = q
            self.mainparent.UpdateQuantitiesMenu()

        self.mainparent.plot_title = "Shell Average of Omega"
        if (self.mainparent.tex):
            self.mainparent.plot_ylabel = "$\Omega / \Omega_0$"
        else:
            self.mainparent.plot_ylabel = "Omega / Omega0"
        self.mainparent.loaded_omega = True
        self.mainparent.current_qnt = self.mainparent.nq

        self.mainparent.statusbar.SetStatusText("Quantity: Omega", 1)

        # plot it
        self.mainparent.plotpane.plot()

    #-----------------------------------------------
    # handle Quit Event
    #-----------------------------------------------
    def Quit(self, event):
        self.mainparent.statusbar.SetStatusText("Quitting ...", 0)
        self.mainparent.Close()

    #-----------------------------------------------
    # handle "Reload/Update Plot" Event
    #-----------------------------------------------
    def Plot(self, event):
        self.mainparent.plotpane.plot()

    #-----------------------------------------------
    # handle "Change Plot Attributes" Event
    #-----------------------------------------------
    def Annotate_Plot(self, event):
        self.mainparent.plotpane.annotate_plot()

    #-----------------------------------------------
    # handle "Clear Figure" Event
    #-----------------------------------------------
    def Clear_Plot(self, event):
        self.mainparent.plotpane.clear_plot()

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

    #-----------------------------------------------
    # display a blank white figure
    #-----------------------------------------------
    def clear_plot(self):

        self.mainparent.statusbar.SetStatusText("Clearing Figure ...", 0)

        self.figure.clf()
        self.axes = self.figure.add_subplot(111)

        # display blank white square
        self.axes.get_xaxis().set_visible(False)
        self.axes.get_yaxis().set_visible(False)

        self.figure.canvas.draw()

    #-----------------------------------------------
    # plot the various quantities
    #-----------------------------------------------
    def plot(self, iter=None, x=None, y=None, alabel=None):

        if (self.mainparent.filename != ""):

            self.mainparent.statusbar.SetStatusText("Plotting ...", 0)

            # shortcuts:
            xmin     = self.mainparent.plot_xmin
            ymin     = self.mainparent.plot_ymin
            xmax     = self.mainparent.plot_xmax
            ymax     = self.mainparent.plot_ymax
            rsol     = self.mainparent.plot_rsol
            title    = self.mainparent.plot_title
            xlabel   = self.mainparent.plot_xlabel
            ylabel   = self.mainparent.plot_ylabel
            radius   = self.mainparent.radius
            leg_loc  = self.mainparent.leg_loc
            show_leg = self.mainparent.show_legend
            nq       = self.mainparent.nq
            nq_menu  = self.mainparent.nq_menu
            nr       = self.mainparent.nr
            nrecs    = self.mainparent.nrecs
            iq       = self.mainparent.current_qnt
            if (iq < len(self.mainparent.data[0,:,0])):
                iqt = iq
                data = self.mainparent.data[:,iqt,:]
            else:
                # should be plotting domega ... so set data iq to 0
                iqt = 0
                data = self.mainparent.data[:,iqt,:]
            iqr      = self.mainparent.quantity_arr[iq]
            iql      = self.mainparent.quantity_arr[nq-1]  # last quantity
            lum_flag = self.mainparent.luminosity
            perturb  = self.mainparent.perturb
            log10    = self.mainparent.log10
            show_avg = self.mainparent.plot_avg
            animate  = self.mainparent.plot_animation

            if (show_avg):
                avg_data = self.mainparent.data_avg

            # adjust the title and labels for Latex Compatibility
            if (self.mainparent.tex):
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

            self.figure.clf()
            self.axes = self.figure.add_subplot(111)
            self.axes.get_xaxis().set_visible(True)
            self.axes.get_yaxis().set_visible(True)

            if (iqr <= iql and (not show_avg) and (not animate)):
                if ((lum_flag) and ((iqr >= 15) and (iqr <= 20))):
                    lum = numpy.empty((nrecs, nr))
                    lum[:,:] = 0.0
                    minv = 0.0
                    maxv = 0.0
                    for ir in range(nrecs):
                        lum[ir,:] = 4.0*numpy.pi*radius[:]*radius[:]*data[:,ir]
                        tmax = numpy.amax(lum[ir,:])
                        tmin = numpy.amin(lum[ir,:])
                        if (tmin < minv):
                            minv = tmin
                        if (tmax < maxv):
                            maxv = tmax

                    # plot
                    l, c = getLineColor(0)
                    lab = "rec = 0"
                    self.axes.plot(radius/rsol, lum[0,:], label=lab,
                                   color=c, linestyle=l)
                    # oplot
                    for i in range(1, nrecs):
                        l, c = getLineColor(i)
                        lab = "rec = "+str(i)
                        self.axes.plot(radius/rsol, lum[i,:], label=lab,
                                       color=c, linestyle=l)

                elif (perturb):
                    minv = 0.0
                    maxv = 0.0
                    for ir in range(nrecs):
                        tmax = numpy.amax(data[:,ir] - data[:,0])
                        tmin = numpy.amin(data[:,ir] - data[:,0])
                        if (tmin < minv):
                            minv = tmin
                        if (tmax < maxv):
                            maxv = tmax
                    # plot
                    l, c = getLineColor(0)
                    lab = "rec: "+str(nrecs-1)+"-0"
                    self.axes.plot(radius/rsol, data[:,nrecs-1]-data[:,0],
                                   label=lab, color=c, linestyle=l)
                    # oplot
                    for i in range(1, nrecs-1):
                        l, c = getLineColor(i)
                        lab = "rec: "+str(i)+"-0"
                        self.axes.plot(radius/rsol, data[:,i]-data[:,0],
                                       label=lab, color=c, linestyle=l)

                elif (log10):

                    # plot
                    l, c = getLineColor(0)
                    lab = "rec: 0"
                    self.axes.plot(radius/rsol, 
                                   numpy.log10(numpy.abs(data[:,0])), 
                                   label=lab, color=c, linestyle=l)
                    # oplot
                    for i in range(1, nrecs):
                        l, c = getLineColor(i)
                        lab = "rec: "+str(i)
                        self.axes.plot(radius/rsol, 
                                       numpy.log10(numpy.abs(data[:,i])),
                                       label=lab, color=c, linestyle=l)

                else:

                    # plot
                    l, c = getLineColor(0)
                    lab = "rec: 0"
                    self.axes.plot(radius/rsol, data[:,0], 
                                   label=lab, color=c, linestyle=l)
                    # oplot
                    for i in range(1, nrecs):
                        l, c = getLineColor(i)
                        lab = "rec: "+str(i)
                        self.axes.plot(radius/rsol, data[:,i],
                                       label=lab, color=c, linestyle=l)

            elif ((not show_avg) and (not animate)):

                if (self.mainparent.domega == None):
                    self.mainparent.ThrowError(\
                             "Must load Omega before trying to plot it")
                    return

                temp = self.mainparent.domega.copy()

                if (perturb):
                    for i in range(1, nrecs):
                        temp[:,i] = temp[:,i] - temp[:,0]
                    temp[:,0] = 0.0

                minv = numpy.amin(temp)
                maxv = numpy.amax(temp)

                #plot
                l, c = getLineColor(0)
                lab = "rec: 0"
                if (perturb):
                    lab = "rec: 0-0"
                self.axes.plot(radius/rsol, temp[:,0], label=lab, color=c,
                               linestyle=l)
                for i in range(1, nrecs):
                    l, c = getLineColor(i)
                    lab = "rec: "+str(i)
                    if (perturb):
                        lab += "-0"
                    self.axes.plot(radius/rsol, temp[:,i], label=lab,
                                   color=c, linestyle=l)

            # plot average
            elif (not animate):

                l, c = getLineColor(0)
                lab = "Avg over records"
                self.axes.plot(radius/rsol, avg_data, label=lab, color=c, 
                               linestyle=l)

            # animate
            else:
                self.axes.plot(x, y, label=alabel, color='r', linestyle='-')


            avg_title_add_on = " (avg over records)"
            ind = string.find(title, avg_title_add_on)
            if (ind > -1):
                self.mainparent.plot_title = title[:ind] # strip the add-on

            self.mainparent.plot_avg = False
            self.mainparent.plot_animation = False

            if (show_leg):
                self.axes.legend(loc=leg_loc)
                leg = self.figure.gca().get_legend()
                #plt.setp(leg.get_texts(), fontsize='small')
                leg.draw_frame(False)

            # title, labels
            self.axes.set_title(title)
            self.axes.set_xlabel(xlabel)
            self.axes.set_ylabel(ylabel)

            # range
            self.axes.set_xlim(xmin, xmax)
            self.axes.set_ylim(ymin, ymax)

            self.figure.canvas.draw()

        else:
            self.mainparent.ThrowError("Must Load a File First")

    #-----------------------------------------------
    # Change the plot attributes
    #-----------------------------------------------
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
        lloc   = wx.StaticText(self, -1, "Legend Loc:")

        t = self.mainparent.plot_title
        xl = self.mainparent.plot_xlabel
        yl = self.mainparent.plot_ylabel
        xm = str(self.mainparent.plot_xmin)
        ym = str(self.mainparent.plot_ymin)
        xM = str(self.mainparent.plot_xmax)
        yM = str(self.mainparent.plot_ymax)
        llt = str(self.mainparent.leg_loc)

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
        self.lltc = wx.TextCtrl(self, 8, llt)

        self.leg = wx.ToggleButton(self, -1, "Show Legend")
        self.leg.Bind(wx.EVT_TOGGLEBUTTON, self.OnShowLeg)
        self.leg.SetValue(self.mainparent.show_legend)
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
        button = wx.BoxSizer(wx.HORIZONTAL)
        ltitles = wx.BoxSizer(wx.HORIZONTAL)

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
        ltitles.Add(lloc, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        ltitles.Add(self.lltc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)

        minbox.Add(xminbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        minbox.Add(yminbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        maxbox.Add(xmaxbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        maxbox.Add(ymaxbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)

        range.Add(minbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        range.Add(maxbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=5)

        button.Add(self.leg, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        button.Add(ok, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        button.Add(cancel, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)

        hbox.Add(titles, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        hbox.Add(xlabels, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        hbox.Add(ylabels, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        hbox.Add(range, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        hbox.Add(ltitles, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)
        hbox.Add(button, 0, wx.EXPAND|wx.ALL|wx.GROW, border=5)

        self.SetSizer(hbox)
        hbox.Fit(self)

        self.Center()

    #-----------------------------------------------
    # button to show/not show the legend
    #-----------------------------------------------
    def OnShowLeg(self, event):
        show = self.leg.GetValue()
        self.mainparent.show_legend = show

    #-----------------------------------------------
    # OK button
    #-----------------------------------------------
    def OnOK(self, event):
        self.Destroy()

    #-----------------------------------------------
    # Cancel buttton
    #-----------------------------------------------
    def OnCancel(self, event):
        self.Destroy()

    #-----------------------------------------------
    # set the values that were changed by user
    #-----------------------------------------------
    def SetValues(self):

        # adjust the proper value with the entered values
        # must convert to string first, GetValue() returns 'unicode'
        self.mainparent.plot_title = str(self.tc.GetValue())
        self.mainparent.plot_xlabel = str(self.xlc.GetValue())
        self.mainparent.plot_ylabel = str(self.ylc.GetValue())
        self.mainparent.leg_loc = str(self.lltc.GetValue())
        xmVal = str(self.xmc.GetValue())
        ymVal = str(self.ymc.GetValue())
        xMVal = str(self.xMc.GetValue())
        yMVal = str(self.yMc.GetValue())

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
    frame = window(None, title="Plot Shell Averages", tex=tex)
    frame.Show()
    app.MainLoop()

    print "\n\n---Complete---\n"


#####################################################################
#  Define a usage program
#####################################################################
def usage():

    print "\nPython GUI to plot CSS shell average data\n"
    print "Usage:\n"
    print "    ./plot_shell_avgs.py [options]\n"
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


