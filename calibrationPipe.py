from astropy.io import fits
import numpy as np
import statistics

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
        # self.image.writeto("zetafile.fts")
        
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
        imagedata = np.divide(imagedata, flatdata)
        return imagedata
    
    def dark_frame(self, imagedata, darkdata):
        imagedata = np.subtract(imagedata, darkdata)
        return imagedata
    
    def stack_images(self, images, calculation):
        if calculation == "median":
            master_image = images[0]
            l = []
            for i in images:
                data = i[0].data
                l.append(data)
                
            for y in np.array(master_image[0].data).tolist():
                for x in y:
                    m = []
                    for i in l:
                        m.append(i[y][x])
                    statistics.median(m)
                
            
            """for i in l:
                xList = []
                for j in range(len(i)):
                    yList = []
                    for k in range(len(i[j])):
                        yList.append(i[j][k])
                    xList.append(yList)"""
                
                    
                    
            
        elif calculation == "average":
            l = []
            for i in images:
                data = i[0].data
                l.append(data)
            
            for i in range(len(images) - 1):
                l[i] = np.add(l[i], l[i + 1])
        
            stacked_images = np.true_divide(l[i], len(images))
            print(stacked_images)      
        
print("start")

a = '2021-04-02T19-10-26_M1_Clear_280s_Simon-H.fts'
b = fits.open('HAT-P-10-001dark.fit')
c = fits.open('Bias-001.fit')
d = fits.open('HAT-P-10-001light.fit')

calibPip = CalibrationPipe(a, b, c, d)
calibPip.calibrate()

imgs = [fits.open(a), b]
calibPip.stack_images(imgs, "median")

print("stop")