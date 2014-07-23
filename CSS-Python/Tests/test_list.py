import wx

class MyCheckDialog(wx.Dialog):
     def __init__(self, Choices):
         wx.Dialog.__init__(self, None, -1, 'Figure Attributes')
         self.Choices = Choices

         # define everything to put in the window
         title = wx.StaticText(self, -1, "title")
         xlabel = wx.StaticText(self, -1, "xlabel")
         ylabel = wx.StaticText(self, -1, "ylabel")
         xmin = wx.StaticText(self, -1, "xmin")
         ymin = wx.StaticText(self, -1, "ymin")
         xmax = wx.StaticText(self, -1, "xmax")
         ymax = wx.StaticText(self, -1, "ymax")

         tc = wx.TextCtrl(self, -1, "")
         xlc = wx.TextCtrl(self, -1, "")
         ylc = wx.TextCtrl(self, -1, "")
         xmc = wx.TextCtrl(self, -1, "")
         ymc = wx.TextCtrl(self, -1, "")
         xMc = wx.TextCtrl(self, -1, "")
         yMc = wx.TextCtrl(self, -1, "")

         ok = wx.Button(self, wx.ID_OK, 'OK')
         cancel = wx.Button(self, -1, 'Cancel')

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

         titles.Add(title, 0, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         titles.Add(tc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         xlabels.Add(xlabel, 0, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         xlabels.Add(xlc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         ylabels.Add(ylabel, 0, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         ylabels.Add(ylc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         xminbox.Add(xmin, 0, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         xminbox.Add(xmc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         yminbox.Add(ymin, 0, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         yminbox.Add(ymc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         xmaxbox.Add(xmax, 0, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         xmaxbox.Add(xMc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         ymaxbox.Add(ymax, 0, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         ymaxbox.Add(yMc, 1, wx.EXPAND|wx.ALL|wx.GROW, border=2)

         minbox.Add(xminbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         minbox.Add(yminbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         maxbox.Add(xmaxbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         maxbox.Add(ymaxbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=2)

         range.Add(minbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         range.Add(maxbox, 1, wx.EXPAND|wx.ALL|wx.GROW, border=2)

         button.Add(ok, 0, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         button.Add(cancel, 0, wx.EXPAND|wx.ALL|wx.GROW, border=2)

         hbox.Add(titles, 0, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         hbox.Add(xlabels, 0, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         hbox.Add(ylabels, 0, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         hbox.Add(range, 0, wx.EXPAND|wx.ALL|wx.GROW, border=2)
         hbox.Add(button, 0, wx.EXPAND|wx.ALL|wx.GROW, border=2)

         self.SetSizer(hbox)
         hbox.Fit(self)

         self.Center() # make it come up on the center of the screen

     def GetChecked(self):
         Checked = []
         for (index, item) in enumerate(self.Choices):
             if self.clb.IsChecked(index):
                 Checked.append(item)
         return Checked


## And to use it:

a = wx.App()

myd = MyCheckDialog(['Title', 'xlabel', 'ylabel', 
                     'xmin', 'xmax', 'ymin', 'ymax'])
if myd.ShowModal() != wx.ID_OK:
     exit(1)
else:
     print "You checked:", myd.GetChecked()

