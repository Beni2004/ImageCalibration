from astropy.io import fits
import numpy as np

class CalibrationPipe():
    def __init__(self, image_path, dark, bias, flat):
        self.image = fits.open(image_path)
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
        self.image.writeto("zetafile.fts")
        
    def verifyHeaderDark(self):
        if self.image[0].header["DATE"] != self.dark[0].header["DATE"]:
            return False
    
    def verifyHeaderBias(self):
        if self.image[0].header["DATE"] != self.bias[0].header["DATE"]:
            return False
    
    def bias_frame(self, imagedata, biasdata):
        imagedata = np.subtract(imagedata, biasdata)
        return imagedata
    
    def flat_frame(self, imagedata, flatdata):
        a = imagedata[10][20]
        imagedata = np.divide(imagedata, flatdata)
        if imagedata[10][20] == a / flatdata[10][20]:
            print("mhm")
        return imagedata
    
    def dark_frame(self, imagedata, darkdata):
        imagedata = np.subtract(imagedata, darkdata)
        return imagedata
              
        
print("start")

a = 'C:/Users/simon/OneDrive/Dokumente/Praktikum CSH/CalibrationPipe/2021-04-02T19-10-26_M1_Clear_280s_Simon-H.fts'
b = fits.open('C:/Users/simon/OneDrive/Dokumente/Praktikum CSH/CalibrationPipe/HAT-P-10-001dark.fit')
c = fits.open('C:/Users/simon/OneDrive/Dokumente/Praktikum CSH/CalibrationPipe/Bias-001.fit')
d = fits.open('C:/Users/simon/OneDrive/Dokumente/Praktikum CSH/CalibrationPipe/HAT-P-10-001light.fit')

calibPip = CalibrationPipe(a, b, c, d)
calibPip.calibrate()

print("stop")