from astropy.io import fits
import numpy as np
import statistics
import time
import multiprocessing as mp
import os
# from multiprocessing import Process
from queue import Queue
# import matplotlib.pyplot as plt

t = time.time()

class CalibrationPipe():
    def __init__(self, image_path, dark_path, bias_path, flat_path):
        with fits.open(image_path) as image:
            self.imagedata = image[0].data
            self.imagehdr = image[0].header
            
        with fits.open(dark_path) as dark:
            self.darkdata = dark[0].data
            self.darkhdr = dark[0].header
        
        with fits.open(bias_path) as bias:
            self.biasdata = bias[0].data
            self.biashdr = bias[0].header
        
        with fits.open(flat_path) as flat:
            self.flatdata = flat[0].data
            self.flathdr = flat[0].header
        
        self.imagepath = image_path
        # self.image.info()
        
    def run(self):
        print("ruuuuun!")
        imagedata = self.imagedata
        darkdata = self.darkdata
        biasdata = self.biasdata
        flatdata = self.flatdata
        
        imagedata = self.dark_frame(imagedata, darkdata, biasdata)
        
        imagedata = self.bias_frame(imagedata, biasdata)
        
        imagedata = self.flat_frame(imagedata, flatdata)
        
        amounty=len(imagedata)
        amountx=len(imagedata[0])
        
        # q1 = Queue()
        # q2 = Queue()
        
        print("Starting multiprocessing...")
        
        processes = []
        for i in [(imagedata, 0, 0, amountx, amounty//2), (imagedata, 0, amounty//2, amountx, amounty//2)]:
            p = mp.Process(target = self.run_check, args = i)
            processes.append(p)
            print("Iteration")
        
        [x.start() for x in processes]
            
        print("All processes should have started now")
        
        [x.join() for x in processes]
        
        # imagedata1 = q1.get()
        # imagedata2 = q2.get()
        
        """ctx = mp.get_context('spawn')
        
        q1 = ctx.Queue()
        q2 = ctx.Queue()
            
        p1 = ctx.Process(target=self.run_check, args=(imagedata, 0, 0, amountx, amounty//2, q1))
        p2 = ctx.Process(target=self.run_check, args=(imagedata, 0, amounty//2, amountx, amounty//2, q2))
        
        p1.start()
        p2.start()
        imagedata_1 = q1.get()
        imagedata_2 = q2.get()
        p1.join()
        p2.join()"""
        
        # imagedata = self.run_check(imagedata, 0, 0, amountx, amounty//2)
        
        # imagedata = self.run_check(imagedata, 0, amounty//2, amountx, amounty//2)
        
        with fits.open(self.imagepath) as image:
            image[0].data = imagedata
            image.writeto("output.fts", overwrite=True)
    
    """def close(self):
        self.image.close()"""
        
        
    def verifyDateDark(self):
        if self.imagehdr["DATE"] != self.darkhdr["DATE"]:
            return False
    
    def verifyDateBias(self):
        if self.imagehdr["DATE"] != self.biashdr["DATE"]:
            return False
        
    def scale_dark(self):
        if self.darkhdr["EXPTIME"] != self.imagehdr["EXPTIME"]:
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

            print(master_image_data)
            master_image[0].data = master_image_data
            master_image.writeto("deltafile.fts", overwrite=True)
                       
        elif calculation == "average":
            l = []
            for i in images:
                data = i[0].data
                l.append(data)
            
            for i in range(len(images) - 1):
                l[i] = np.add(l[i], l[i + 1])
        
            stacked_images = np.true_divide(l[i], len(images))
            master_image[0].data = stacked_images
            master_image.writeto("gammafile.fts", overwrite=True)
            print(stacked_images)
            
    def remove_hotpixel(self, x, y, imagedata):
        #Sum=imagedata[y, x-1].astype("uint")+imagedata[y+1, x].astype("uint")+imagedata[y, x+1].astype("uint")+imagedata[y-1, x].astype("uint")
        Sum=imagedata[y, x-1].astype("uint")+imagedata[y, x+1].astype("uint")
        av=Sum//2
        value=imagedata[y, x]
        if value > 1.25*av:
            imagedata[y, x]=av
        return imagedata
        
    def run_check(self, imagedata, startx, starty, amountx, amounty):
        print("A process was started")
        for y in range(amounty-1):
            for x in range(amountx):
                try:
                    data = self.remove_hotpixel(x + startx, y + starty, imagedata)
                except IndexError:
                    pass
        print(data)
        # q.put(data)  

if __name__ == '__main__':
    mp.set_start_method('spawn')

    print("start")

    image_path = '2021-04-02T19-10-26_M1_Clear_280s_Simon-H.fts'
    b = 'HAT-P-10-001dark.fit'
    c = 'Bias-001.fit'
    d = 'HAT-P-10-001light.fit'
    e = "Reg_2021-03-31T23-55-09_M51_Red_280s_Benjamin-A.fit"

    calibPip = CalibrationPipe(image_path, b, c, d)
    # calibPip.close()
    calibPip.run()

    """timeStats = []
    startTime = time.time()
    imgs = [c, c]
    calibPip.stack_images(imgs, "median")
    timeStats.append(time.time() - startTime)

    startTime = time.time()
    imgs = [c, c, c]
    calibPip.stack_images(imgs, "median")
    timeStats.append(time.time() - startTime)

    startTime = time.time()
    imgs = [c, c, c, c]
    calibPip.stack_images(imgs, "median")
    timeStats.append(time.time() - startTime)

    plt.plot([1,2,3], timeStats)
    plt.show()"""

    """now=time.time()
    imgs = [fits.open(image_path), e]
    calibPip.stack_images(imgs, "median")
    calibPip.stack_images(imgs, "average")"""


    # calibPip.close()

    # print('stashing: ' + str(time.time() - now))
    print('total: ' + str(time.time() - t))

    print("stop")