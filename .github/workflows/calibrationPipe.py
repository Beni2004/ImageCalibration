from astropy.io import fits
import numpy as np
import matplotlib as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
# from astropy.io import units
# from astropy.io import astropy.nddata import CCDData
# import ccdproc

class CalibrationPipe():
    def __init__(self, image, dark, bias, flat):
        self.image = image
        self.dark = dark
        self.bias = bias
        self.flat = flat
    
    def calibrate(self):
        imagedata = self.image[0].data
        darkdata = self.dark[0].data
        biasdata = self.bias[0].data
        flatdata = self.flat[0].data
        
        imagedata = self.dark_frame(imagedata, darkdata)
        
        imagedata = self.bias_frame(imagedata, biasdata)
        
        imagedata = self.flat_frame(imagedata, flatdata)
        
        self.image[0].data = imagedata
        self.image.writeto("deltafile.fts")
        # if imagedata[10][20] == a - darkdata[10][20]:
            # print("ok")
    
    def bias_frame(self, imagedata, biasdata):
        imagedata = np.subtract(imagedata, biasdata)
        return imagedata
    
    def flat_frame(self, imagedata, flatdata):
        imagedata = np.divide(imagedata, flatdata)
        return imagedata
    
    def dark_frame(self, imagedata, darkdata):
        imagedata = np.subtract(imagedata, darkdata)
        return imagedata
              
        
print("start")

a = fits.open('C:/Users/simon/OneDrive/Dokumente/Praktikum CSH/CalibrationPipe/gammafile.fts')
b = fits.open('C:/Users/simon/OneDrive/Dokumente/Praktikum CSH/CalibrationPipe/betafile.fts')

calibPip = CalibrationPipe(a, b)
calibPip.calibrate()

data = a[0].data
# print(data)

print("stop")