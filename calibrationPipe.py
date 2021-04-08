from astropy.io import fits
import numpy as np
import statistics
import time
from multiprocessing import Process

t = time.time()

class CalibrationPipe():
    def __init__(self, image_path, dark, bias, flat):
        self.image = fits.open(image_path)
        self.dark = dark
        self.bias = bias
        self.flat = flat
        self.image.info()
        
    def calibrate(self):
        imagedata = self.image[0].data
        darkdata = self.dark[0].data
        biasdata = self.bias[0].data
        flatdata = self.flat[0].data
        
        amounty=len(imagedata)
        amountx=len(imagedata[0])
        
        imagedata = self.run_check(imagedata, amountx, amounty)
        
        imagedata = self.dark_frame(imagedata, darkdata, biasdata)
        
        imagedata = self.bias_frame(imagedata, biasdata)
        
        imagedata = self.flat_frame(imagedata, flatdata)
        
        self.image[0].data = imagedata
        self.image.writeto("output.fts")
    
    def close(self):
        self.image.close()
        
        
    def verifyDateDark(self):
        if self.image[0].header["DATE"] != self.dark[0].header["DATE"]:
            return False
    
    def verifyDateBias(self):
        if self.image[0].header["DATE"] != self.bias[0].header["DATE"]:
            return False
        
    def scale_dark(self):
        if self.dark[0].header["EXPTIME"] != self.image[0].header["EXPTIME"]:
            # ...
            pass
    
    def bias_frame(self, imagedata, biasdata):
        imagedata = np.subtract(imagedata, biasdata)
        return imagedata
    
    def flat_frame(self, imagedata, flatdata):
        imagedata = np.divide(imagedata, flatdata)
        return imagedata
    
    def dark_frame(self, imagedata, darkdata, biasdata):
        darkdata = np.subtract(darkdata, biasdata)
        imagedata = np.subtract(imagedata, darkdata)
        return imagedata
    
    def stack_images(self, images, calculation):
        master_image = images[0]
        master_image_data = master_image[0].data
        
        if calculation == "median":
            l = []
            for i in images:
                l.append(i[0].data.astype("uint"))
                
            for y in range(4096):
                for x in range(4096):
                    m = []
                    for i in l:
                        m.append(i[x][y])
                    master_image_data[x][y] = statistics.median(m)
                    # print("iteration")
            print(master_image_data)
            master_image[0].data = master_image_data
            master_image.writeto("deltafile.fts")
                       
        elif calculation == "average":
            l = []
            for i in images:
                data = i[0].data
                l.append(data)
            
            for i in range(len(images) - 1):
                l[i] = np.add(l[i], l[i + 1])
        
            stacked_images = np.true_divide(l[i], len(images))
            master_image[0].data = stacked_images
            master_image.writeto("gammafile.fts")
            print(stacked_images)
            
    def remove_hotpixel(self, x, y, imagedata):
        Sum=imagedata[y, x-1].astype("uint")+imagedata[y+1, x].astype("uint")+imagedata[y, x+1].astype("uint")+imagedata[y-1, x].astype("uint")
        #Sum=imagedata[y, x-1].astype("uint")+imagedata[y, x+1].astype("uint")
        average=Sum//4
        value=imagedata[y, x]
        if value>1.5*average:
            imagedata[y, x]=average
        return imagedata
        
        
    def run_check(self, imagedata, amountx, amounty):
        for y in range(amounty-1):
            for x in range(amountx-1):
                data = self.remove_hotpixel(x, y, imagedata)
        return data  
        
print("start")

image_path = '2021-04-02T19-10-26_M1_Clear_280s_Simon-H.fts'
b = fits.open('HAT-P-10-001dark.fit')
c = fits.open('Bias-001.fit')
d = fits.open('HAT-P-10-001light.fit')
e = fits.open("Reg_2021-03-31T23-55-09_M51_Red_280s_Benjamin-A.fit")

calibPip = CalibrationPipe(image_path, b, c, d)
calibPip.calibrate()


imgs = [fits.open(image_path), e]
calibPip.stack_images(imgs, "median")
print(time.time() - t)
calibPip.stack_images(imgs, "average")

calibPip.close()

print("stop")