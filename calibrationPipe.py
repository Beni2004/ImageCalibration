from astropy.io import fits
import numpy as np
import statistics

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
        master_image = images[0]
        master_image_data = master_image[0].data
        
        if calculation == "median":
            l = []
            for i in images:
                l.append(i[0].data.astype("uint"))
                
            for y in range(4095):
                for x in range(4095):
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
        
print("start")

a = '2021-04-02T19-10-26_M1_Clear_280s_Simon-H.fts'
b = fits.open('HAT-P-10-001dark.fit')
c = fits.open('Bias-001.fit')
d = fits.open('HAT-P-10-001light.fit')
e = fits.open("Reg_2021-04-03T19-27-40_M51_Red_200s_Simon-H.fit")

calibPip = CalibrationPipe(a, b, c, d)
calibPip.calibrate()

imgs = [fits.open(a), e]
calibPip.stack_images(imgs, "median")
calibPip.stack_images(imgs, "average")

print("stop")