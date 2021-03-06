# -*- coding: utf-8 -*-
import numpy as np
import os,sys
import pylab as plt
from pprint import pprint

cmap = 'gray'
from matplotlib.widgets import Slider, Button

def write_minc(volume, vox_size, filename):
    """Write out a data volume in minimal MINC format, using h5py."""
    import h5py
    if os.path.exists(filename):
        os.remove(filename)
    out = h5py.File(filename,'w')
    d = out.create_group('dimensions')
    d['xspace'] = out.create_dataset('xspace',shape=(),dtype="<f8")
    d['xspace'].attrs['length'] = volume.shape[2]
    d['xspace'].attrs['step'] = vox_size[2]
    d['yspace'] = out.create_dataset('yspace',shape=(),dtype="<f8")
    d['yspace'].attrs['length'] = volume.shape[1]
    d['yspace'].attrs['step'] = vox_size[1]
    d['zspace'] = out.create_dataset('zspace',shape=(),dtype="<f8")
    d['zspace'].attrs['length'] = volume.shape[0]
    d['zspace'].attrs['step'] = vox_size[0]
    g0 = out.create_group('0')
    g0['image'] = out.create_dataset('image',shape=volume.shape,dtype="<f8")
    g0['image'].attrs['dimorder'] = 'zspace,yspace,xspace'
    g0['image'][()] = volume.astype(np.float32)
    g0['image-max'] = out.create_dataset('image-max',shape=(),dtype="<f8")
    g0['image-max'][()] = volume.max()
    g0['image-min'] = out.create_dataset('image-min',shape=(),dtype="<f8")
    g0['image-min'][()] = volume.min()
    im = out.create_group('im')
    im['0'] = g0
    info = out.create_group('info')
    main = out.create_group('minc-2.0')
    main['info'] = info
    main['image'] = im
    main['dimensions'] = d
    out.flush()
    out.close()

def read_minc(filename):
    """Read a data volume in MINC format, using h5py."""
    import h5py
    out = h5py.File(filename,'r')
    v = out['minc-2.0']['image']['0']['image'].value
    z = [out['minc-2.0']['dimensions'][s].attrs['step'] for s in ['zspace','yspace','xspace']]
    return v.astype('float32'),z
    

class display_many:
    """As display, above, but instead of selecting between data-sets, show them
    side-by side, and provide a graphical cursor. Note that selecting time-point
    of -1 shows time1-time0 subtracted image."""
    def __init__(self,*arg,**kwargs):
        """Optional parameters:
        - fig: the matplotlib figure to use ("current" if left blank)
        z-slice selector slider; otherwise use array index.
        """
        fig = kwargs.get('fig', None )
        if fig==None: fig = plt.gcf()
        self.fig=fig
        fig.clf()
        self.text = fig.text(0,.99,"",va='top')
        self.num = len(arg)
        self.row,self.col = [1,1]
        fig.subplots_adjust(left=0.05, bottom=0.25,top=0.95,right=0.95,wspace=0.1,hspace=0.15)
        args = []
        self.zs = []
        
        self.image_index = 0
        
        # initialize using an array of images
        for a in arg:
            args.append(a)
            self.zs.append(np.arange(len(a)))
        
        maxz = max([max(z) for z in self.zs])
        minz = min([min(z) for z in self.zs])
        self.data = args
        self.im = []
        self.ax = []
        
        for i in range(len(self.data)):
            self.ax.append(fig.add_subplot(self.row,self.col,i+1))
            self.im.append(self.ax[-1].imshow(self.data[i][0],origin='upper',cmap=cmap,aspect='equal'))
            
        self.cmap=cmap
        axcolor = 'lightgoldenrodyellow'
        self.sliceax = fig.add_axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
        self.minax  = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
        self.maxax  = fig.add_axes([0.25, 0.05, 0.65, 0.03], axisbg=axcolor)
        self.sli = 0
        self.sli_wid = Slider(self.sliceax, 'Slice', minz, maxz, valinit=minz,valfmt="%i")
        self.vmin = Slider(self.minax, 'min', 0, 100, valinit=0,valfmt="%i")
        self.vmax = Slider(self.maxax, 'max', 0, 100, valinit=100,valfmt="%i")
        self.sli_wid.on_changed(self.update)
        self.vmin.on_changed(self.update)
        self.vmax.on_changed(self.update)
        self.multi = XYMultiCursor(self.fig.canvas, self.ax, color='b', lw=1, ls=':')
        self.fig.canvas.mpl_connect('button_press_event',self.key)
        self.fig.canvas.draw()
        self.points = {}
        
        # defines the positions of the previous/next slice buttons
        axprev = plt.axes([0.81, 0.9, 0.1, 0.075])
        axnext = plt.axes([0.81, 0.8, 0.1, 0.075])
        self.bprev = Button(axprev, 'Previous')
        self.bnext = Button(axnext, 'Next')
        # previous/next slice select buttons
        self.bprev.on_clicked(self.next_image)
        self.bnext.on_clicked(self.prev_image)
        pprint(self.data[0][0])
        
    def update(self,val):
        if self.vmin.val<self.vmax.val:
            vmin,vmax = self.vmin.val,self.vmax.val
            cm = self.cmap
        elif self.vmin.val>self.vmax.val:
            vmax,vmin = self.vmin.val,self.vmax.val
            cm = self.cmap+'_r'
        for i in range(self.num):
            z = np.argmin((self.zs[i] - self.sli_wid.val)**2)
            im = self.data[i][z]
            print "self.data acessed - i: {} z: {}".format(i, z)

            self.im[i].set_data(im)
            self.im[i].set_cmap(cm)
            self.im[i].set_clim(im.min()+(im.max()-im.min())*vmin/100.,im.min()+(im.max()-im.min())*vmax/100.)
        self.fig.canvas.draw()
        
    def next_image(self, event):
        ind = self.image_index + 1
        self.image_index = ind % len(self.data[0])
        self.set_image(self.image_index)

    def prev_image(self, event):
        ind = self.image_index - 1
        self.image_index = ind % len(self.data[0])
        self.set_image(self.image_index)
        
    def set_image(self, image_index):
        if self.vmin.val<self.vmax.val:
            vmin,vmax = self.vmin.val,self.vmax.val
            cm = self.cmap
            
        if self.vmin.val>self.vmax.val:
            vmax,vmin = self.vmin.val,self.vmax.val
            cm = self.cmap+'_r'        
        
        im = self.data[0][image_index]
        self.im[0].set_data(im)
        self.im[0].set_cmap(cm)
        self.im[0].set_clim(im.min()+(im.max()-im.min())*vmin/100.,im.min()+(im.max()-im.min())*vmax/100.)   
        self.fig.canvas.draw()
    
    def key(self,event):
        if event.button != 2:
            return
        if event.inaxes in self.ax:
            z = np.argmin((self.zs[self.ax.index(event.inaxes)] - self.sli_wid.val)**2)
            self.raw = (z,event.ydata,event.xdata)
            print "\nRaw: %i, %.1d, %.1d"%self.raw
            sys.stdout.flush()

class XYMultiCursor:
    """
    Provide a XY line cursor shared between multiple axes (see matplotlib.MultiCursor)
    """
    def __init__(self, canvas = None, axes = None, **lineprops):

        self.canvas = canvas or plt.gcf().canvas
        self.axes = axes or plt.gcf().axes
        xmin, xmax = axes[-1].get_xlim()
        xmid = 0.5*(xmin+xmax)
        ymin, ymax = axes[-1].get_ylim()
        ymid = 0.5*(ymin+ymax)
        self.linesx = [ax.axvline(xmid, visible=False, **lineprops) for ax in axes]
        self.linesy = [ax.axhline(ymid, visible=False, **lineprops) for ax in axes]
        self.visible = True
        self.background = None
        self.needclear = False
        self.canvas.mpl_connect('motion_notify_event', self.onmove)

    def onmove(self, event):
        if not(event.inaxes in self.axes): return
        for line in self.linesx:
            line.set_xdata((event.xdata, event.xdata))
            line.set_visible(self.visible)
        for line in self.linesy:
            line.set_ydata((event.ydata, event.ydata))
            line.set_visible(self.visible)
        self.canvas.draw()

def test_display_many():
    import blood_tools
    scans_to_process=['/../3T_data/23apr2016_trio/3D SE blood 2_24']

    #scans_to_process=['27feb2016_trio/T1map_TRUFI_40DEG_74',
    #                  '10mar2016_trio/T1map_TRUFI_40DEG shorter_69',
    #                  '30mar2016_trio/T1map_TRUFI_40DEG shorter_78',
    #                  '23apr2016_trio/T1map_TRUFI_40DEG shorter_82',
    #                  '27apr2016_trio/T1map_TRUFI_40DEG shorter_82']
                      
               
    parent_folder= os.getcwd()
    for jj, scan in enumerate(scans_to_process): 
        folder = parent_folder + scan
        images,TI=blood_tools.read_dicoms(folder,['InversionTime'])
        
        display_many(images)
    
if __name__ == '__main__':
    test_display_many()