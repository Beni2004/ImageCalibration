from astropy.io import fits
import numpy as np
import statistics
import time
import multiprocessing as mp
import math
import matplotlib.pyplot as plt

t = time.time()

progress=0

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
        
        
    def run(self):
        print("ruuuuun!")
        imagedata = self.imagedata
        darkdata = self.darkdata
        biasdata = self.biasdata
        flatdata = self.flatdata
        
        amounty=len(imagedata)
        amountx=len(imagedata[0])
        
        # darkdata = self.scale_dark(darkdata, amountx, amounty)
        
        q1 = mp.Queue()
        q2 = mp.Queue()
        q3 = mp.Queue()
        q4 = mp.Queue()
        
        processes = []
        for i in [(darkdata, 0, 0, amountx, amounty//4, q1), (darkdata, 0, amounty//4, amountx, amounty//4, q2), (darkdata, 0, amounty//2, amountx, amounty//4, q3), (darkdata, 0, 3*amounty//4, amountx, amounty//4, q4)]:
            p = mp.Process(target = self.scale_dark, args = i)
            processes.append(p)
            p.start()
        
        darkdata1 = q1.get()
        darkdata2 = q2.get()
        darkdata3 = q3.get()
        darkdata4 = q4.get()
        [x.join() for x in processes]
        
        darkdata = np.concatenate((darkdata1[0:amounty//4], darkdata2[amounty//4:amounty//2], darkdata3[amounty//2:3*amounty//4], darkdata4[3*amounty//4:amounty]), axis=0)
        
        # imagedata = self.dark_frame(imagedata, darkdata, biasdata)
        
        # imagedata = self.bias_frame(imagedata, biasdata)
        
        # imagedata = self.flat_frame(imagedata, flatdata)
        
        q1 = mp.Queue()
        q2 = mp.Queue()
        q3 = mp.Queue()
        q4 = mp.Queue()
        
        print("Starting multiprocessing...")
        t_start = time.time()
        processes = []
        for i in [(imagedata, 0, 0, amountx, amounty//4, q1), (imagedata, 0, amounty//4, amountx, amounty//4, q2), (imagedata, 0, amounty//2, amountx, amounty//4, q3), (imagedata, 0, 3*amounty//4, amountx, amounty//4, q4)]:
            p = mp.Process(target = self.run_check, args = i)
            processes.append(p)
            print("Iteration")
        
        [x.start() for x in processes]
            
        print("All processes should have started now")
        imagedata1 = q1.get()
        imagedata2 = q2.get()
        imagedata3 = q3.get()
        imagedata4 = q4.get()
        [x.join() for x in processes]
        
        imagedata = np.concatenate((imagedata1[0:amounty//4], imagedata2[amounty//4:amounty//2], imagedata3[amounty//2:3*amounty//4], imagedata4[3*amounty//4:amounty]), axis=0)
        
        """plt.imshow(imagedata)
        plt.colorbar(orientation='vertical')
        plt.show()"""
        
        print("imgdat1:", imagedata1)
        print("imgdat2:", imagedata2)
        print("removal time:", time.time()-t_start)
        
        with fits.open(self.imagepath) as image:
            image[0].data = imagedata
            image.writeto("output.fts", overwrite=True)
        
        
    def verifyDateDark(self):
        if self.imagehdr["DATE"] != self.darkhdr["DATE"]:
            return False
    
    def verifyDateBias(self):
        if self.imagehdr["DATE"] != self.biashdr["DATE"]:
            return False
        
    def scale_dark(self, darkdata, startx, starty, amountx, amounty, q):
        if self.darkhdr["EXPTIME"] != self.imagehdr["EXPTIME"]:
            scaling_factor = 1
            delta_t = math.sqrt((self.darkhdr["EXPTIME"] - self.imagehdr["EXPTIME"]) ** 2)
            
            for y in range(amounty):
                for x in range(amountx):
                    darkdata[startx + x][starty + y] += (scaling_factor * delta_t)
            
            print(darkdata)
            q.put(darkdata)
        
        return darkdata
    
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
            
    def remove_hotpixel(self, x, y, imagedata):
        #Sum=imagedata[y, x-1].astype("uint")+imagedata[y+1, x].astype("uint")+imagedata[y, x+1].astype("uint")+imagedata[y-1, x].astype("uint")
        Sum=imagedata[y, x-1].astype("uint")+imagedata[y, x+1].astype("uint")
        av=Sum//2
        value=imagedata[y, x]
        if value > 1.25*av:
            imagedata[y, x]=av
        return imagedata
        
    def run_check(self, imagedata, startx, starty, amountx, amounty, q):
        global progress
        total_area=amountx*amounty
        print("A process was started")
        for y in range(amounty-1):
            for x in range(amountx):
                try:
                    data = self.remove_hotpixel(x + startx, y + starty, imagedata)
                except IndexError:
                    pass
                area=(x-startx)*(y-starty)
                percent=100*area/total_area
                if percent - progress >=5:
                    progress=round(percent)
                    print(str(progress) + " %")
        print("Process finished")
        q.put(data)
        
class ImageCombiner():
    def __init__(self):
        pass
        
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
            master_image.writeto("masterfile.fts", overwrite=True)
                       
        elif calculation == "average":
            l = []
            for i in images:
                data = i[0].data
                l.append(data)
            
            for i in range(len(images) - 1):
                l[i] = np.add(l[i], l[i + 1])
        
            stacked_images = np.true_divide(l[i], len(images))
            master_image[0].data = stacked_images
            master_image.writeto("masterfile.fts", overwrite=True)
            print(stacked_images)

if __name__ == '__main__':
    mp.set_start_method('spawn')

    print("start")
    print(str(progress) + " %")

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