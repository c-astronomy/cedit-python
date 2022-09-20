import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from astropy.visualization import simple_norm

#import wx
plt.style.use(astropy_mpl_style)

import time

from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
#%matplotlib inline    #tmp, dont know how it works 

from astropy.utils.data import download_file
#image_file = download_file('http://data.astropy.org/tutorials/FITS-images/HorseHead.fits', cache=True)

#from astropy.visualization import SqrtStretch
from astropy.visualization import (SqrtStretch, MinMaxInterval, ImageNormalize)

from astropy.visualization.stretch import SinhStretch, LinearStretch


import numpy as np

hdu_list = fits.open('/home/kali/Git/cedit-python/images/crescent.fit')



#Read fits data
def read():
    
    #with fits.open('/home/kali/Git/cedit-python/images//HorseHead.fits') as hdu:
        #science_data1 = hdu[1].data
        #science_header1 = hdu[1].header

        #error_data2 = hdu['err',2].data
        #error_header2 = hdy['err',2].header

    hdu_list = fits.open('/home/kali/Git/cedit-python/images/crescent.fit')
    #hdu_list = fits.open(image_file)
    hdu_list.info()

    image_data = hdu_list[0].data

    print(type(image_data))
    print(image_data.shape)

    hdu_list.close()

    #to skip steps just do
    # image_data = fits.getdata(image_file)
    # print(type(image_data))
    # print(image_data,shape)

    plt.imshow(image_data, cmap='gray')
    plt.colorbar()

    #get stats
    print('Min:', np.min(image_data))
    print('Max:', np.max(image_data))
    print('Mean:', np.mean(image_data))
    print('Stdev:', np.std(image_data))




#Load image
def loadimage():
    image_file = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')
    image_data = fits.getdata(image_file, ext=0)
    #histogram = plt.hist(image_data.flatten(), bins='auto')
    plt.imshow(image_data, cmap='gray')
    
    plt.grid(False)
    plt.show()

#another load image
def loadimagetwo():
    hdu_list = fits.open('/home/kali/Git/cedit-python/images/crescent.fit')
    hdu_list.info()

    image_data = hdu_list[0].data

    print(type(image_data))
    print(image_data.shape)

    #hdu_list.close()

    #to skip steps just do
    # image_data = fits.getdata(image_file)
    # print(type(image_data))
    # print(image_data,shape)


    #image_hist = plt.hist(image_data.flatten(), bins='auto')
    #plt.imshow(image_data, cmap='gray')
    #ifits.mpl_normalize(image_data, stretch=LinearStretch,clip = False)
    #norm = simple_norm(image_data, 'sqrt')
    #stretch = SqrtStretch()
    #stretch([-1., 0., 0.5, 1., 1.5],clip = False)
    
    norm = ImageNormalize(image_data, interval = MinMaxInterval(), stretch = SqrtStretch())
    
    stretch = LinearStretch(slope = 0.5, intercept = 0.5) + SinhStretch() +  LinearStretch(slope = 2, intercept = -1)
    #norm = ImageNormalize(stretch=stretch, vmin = -5, vmax = 5)
    #stretch = SqrtStretch()
    #stretch = SqrtStretch()
    #stretch([-1., 0., 0.5, 1., 1.5],clip = False)
    

    #plt.imshow(image_data, cmap='gray', vmin=2E3, vmax=3E3)
    
    plt.imshow(image_data, cmap='gray', origin = 'lower',  norm=norm)
    
    plt.colorbar()
    plt.grid(False)
    plt.show()

    #get stats
    print('Min:', np.min(image_data))
    print('Max:', np.max(image_data))
    print('Mean:', np.mean(image_data))
    print('Stdev:', np.std(image_data))

    hdu_list.close()


#stacking
def stack():
    #image_file = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')
    #Load files
    base_url = 'http://data.astropy.org/tutorials/FITS-images/M13_blue_{0:04d}.fits'
    image_list = [download_file(base_url.format(n), cache=True)
    for n in range(1, 5+1)]
    image_concat = [fits.getdata(image) for image in image_list]

    #stack
    final_image = np.zeros(shape=image_concat[0].shape)
    for image in image_concat:
        final_image += image  #Might use other algo later

    #Short version
    # final_image = np.sum(image_concat, axis=0)

    #Comment all out, extra for stretch
    image_hist = plt.hist(final_image.flatten(), bins='auto')

    plt.imshow(final_image, cmap='gray', vmin=2E3, vmax=3E3)
    plt.colorbar()
    plt.grid(False)
    plt.show()

    #plt.imshow(image_data, cmap='gray')
    #plt.colorbar()
    #plt.grid(False)
    #plt.show()



#histogram
#def histogram():





#stretch
def stretch():
    image_data = fits.getdata(image_file, ext=0)
    image_hist = plt.hist(image_data.flatten(), bins='auto')
    #plt.imshow(image_data, cmap='gray', vmin=2E3, vmax=3E3)
    #histogram = plt.hist(image_data.flatten(), bins='auto')
    plt.imshow(image_data, cmap='gray')
    plt.colorbar()
    plt.grid(False)
    plt.show()



#write image
def write():
    outfile = 'stacked_M13_blue.fits'
    hdu = fits.PrimaryHDU(final_image)
    hdu.writeto(outfile, overwrite=True)



#plot
#def plotti():
    #image_data = fits.getdata(image_file, ext=0)
    #histogram = plt.hist(image_data.flatten(), bins='auto')
    #plt.imshow(image_data, cmap='gray')
    
    #plt.grid(False)
    #plt.show()
    #print("hej")

#read()
#stack()
#loadimage()
loadimagetwo()
