#Contributors: Benjamin A. and Simon H.
#For Questions, feel free to contact imagecalibration-contributors@outlook.com
#If something isn't working as it's supposed to, try:
    #check if you have all the modules installed
    #check if you declared your filepaths correctly
    #Have a look at line 284
#While the program is running, it only prints the progress of the hotpixel removal, not the progress of the whole process.
from astropy.io import fits
import numpy as np
import statistics
import time
import multiprocessing as mp
#import math
#import matplotlib.pyplot as plt

#if debug is 'True' it prints some helpful stuff, otherwise not.
debug = True

t = time.time()

#number of processes executed at the same time
n=32

progress=0

#wether or not the program should remove hotpixels or not, depends on wether you need it or not (needs a lot of time, ～1.5-2 min.).
remove_hotpixels = True

#wether or not the program should do the calibration, i.e. bias, dark, flat, depends on wether you need it and have the necessary files.
calibrate_with_bias = True
calibrate_with_dark = False
calibrate_with_flat = True

class CalibrationPipe():
    #opens all images with their path
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
        
    #Executes all processes, mainly the multiprocessing for removing hotpixels    
    def run(self):
        if debug:
            print("ruuuuun!")
        imagedata = self.imagedata
        darkdata = self.darkdata
        biasdata = self.biasdata
        flatdata = self.flatdata
        
        amounty=len(imagedata)
        amountx=len(imagedata[0])
        
        if calibrate_with_dark:
            if debug:
                print("Starting multiprocessing...")
            q = mp.Queue()
            
            processes = []
            runs = self.make_processes(darkdata, amountx, amounty, n, q)
            for i in runs:
                p = mp.Process(target = self.scale_dark, args = i)
                processes.append(p)
                if debug:
                    print("Iteration")
            
            for x in processes:
                x.start()
                
            for i in range(n):
                tupel = q.get()
                current = tupel[0]
                pos = tupel[1]
                darkdata = np.concatenate((darkdata[0:pos*amounty//n], current[pos*amounty//n:(pos+1)*amounty//n], darkdata[(pos+1)*amounty//n:amounty]), axis=0)
            
            [x.join() for x in processes]
            
            imagedata = self.dark_frame(imagedata, darkdata, biasdata)
        
        if calibrate_with_bias:
            imagedata = self.bias_frame(imagedata, biasdata)
            
        if calibrate_with_flat:
            imagedata = self.flat_frame(imagedata, flatdata)
        
        if remove_hotpixels:
            if debug:
                print("Starting multiprocessing...")
            q = mp.Queue()
            t_start = time.time()
            processes = []
            runs = self.make_processes(imagedata, amountx, amounty, n, q)
            for i in runs:
                p = mp.Process(target = self.clean_image, args = i)
                processes.append(p)
                if debug:
                    print("Iteration")
            
            for x in processes:
                x.start()
                
            if debug:
                print("All processes should have started now")

            for i in range(n):
                tupel = q.get()
                current = tupel[0]
                pos = tupel[1]
                imagedata = np.concatenate((imagedata[0:pos*amounty//n], current[pos*amounty//n:(pos+1)*amounty//n], imagedata[(pos+1)*amounty//n:amounty]), axis=0)
            
            if debug:
                print("removal time:", time.time()-t_start)
        
        with fits.open(self.imagepath) as image:
            image[0].data = imagedata
            image.writeto("output.fts", overwrite=True)
            
    def make_processes(self, imagedata, amountx, amounty, number, q):
        runs = []
        for i in range(number):
            process = (imagedata, 0, i*amounty//number, amountx, amounty//number, q, i)
            runs.append(process)
        return runs
        
        
    def verifyDateDark(self):
        if self.imagehdr["DATE"] != self.darkhdr["DATE"]:
            return False
        return True
    
    def verifyDateBias(self):
        if self.imagehdr["DATE"] != self.biashdr["DATE"]:
            return False
        return True
        
    def scale_dark(self, darkdata, startx, starty, amountx, amounty, q, p):
        #Function that scales the dark.
        #ARGUMENTS: The dark as nparray,
        #the x and y coordinate of the bottom left corner of the area of the dark that should be processed,
        #the x and y scale of the area that should be processed and the multiprocessing queue
        
        if debug:
            print("A process was started")
        
        expt = self.imagehdr["EXPTIME"]
        scaling_factor = (-0.000629169*(expt ** 2) + 0.147650091*expt + 1023.13674955)/self.darkhdr["EXPTIME"]
        #This calculates the average pixelvalue that a dark with given exposuretime should have divided by the actual exposuretime of the dark
        #The result is the scaling factor. Each pixel then has to be multiplied with that value.
            
        for y in range(amounty):
            for x in range(amountx):
                darkdata[startx + x][starty + y] *= scaling_factor
            
        #if debug:
            #print(darkdata)
        q.put((darkdata, p))
        if debug:
            print("Process finished")
    
    def bias_frame(self, imagedata, biasdata):
        #Subtracts the BIAS-frame from the image
        #Takes two nparrays as arguments
        imagedata = np.subtract(imagedata, biasdata)
        return imagedata
    
    def flat_frame(self, imagedata, flatdata):
        #Divides the image with the flat-frame
        #This function should be called after subtracting the dark and the BIAS
        #Takes two nparrays as arguments
        imagedata = np.divide(imagedata, flatdata)
        return imagedata
    
    def dark_frame(self, imagedata, darkdata, biasdata):
        #Subtracts the BIAS-frame from the darkframe and then subtracts the result from the image
        #In other words: Calibrates the dark with the BIAS-frame and then subtracts the dark.
        #Takes two nparrays as arguments
        darkdata = np.subtract(darkdata, biasdata)
        imagedata = np.subtract(imagedata, darkdata)
        return imagedata
            
    #Removes a hotpixel, if there is one. Looks at the pixels left and right of the pixel.
    #The brightness of a hotpixel is set to the average of its two neighbors if the following conditions are true:
    #The hotpixel needs to be at least 5 % brighter than the average of the neighbors.
    #The two neighboring pixels have to be darker than 5000 (out of the max brightness 256^2 - 1).
    def remove_hotpixel(self, x, y, imagedata):
        #Sum=imagedata[y, x-1].astype("uint")+imagedata[y+1, x].astype("uint")+imagedata[y, x+1].astype("uint")+imagedata[y-1, x].astype("uint")
        left=imagedata[y, x-1].astype("uint")
        right=imagedata[y, x+1].astype("uint")
        Sum=left+right
        av=Sum//2
        value=imagedata[y, x]
        if value > 1.05*av and av < 5000:
            imagedata[y, x]=av
        return imagedata
    
    #Runs the hotpixel check for every pixel in the image, except if the maximum y-coordinate is reached.
    #Four processes are running simultaniously due to multiprocessing. Hence the y-axis is divided four times.
    #It calculates the progress of the current process. If the process has gone further than 5 % since the last print, it prints the progress.
    def clean_image(self, imagedata, startx, starty, amountx, amounty, q, p):
        global progress
        total_area=amountx*amounty
        if debug:
            print("A process was started")
        for y in range(amounty):
            for x in range(amountx):
                try:
                    data = self.remove_hotpixel(x + startx, y + starty, imagedata)
                except IndexError:
                    pass
                area=(x-startx)*(y-starty)
                percent=100*area/total_area
                if percent - progress >=5:
                    progress=round(percent)
                    if debug:
                        print(str(progress) + " %")
        q.put((data, p))
        if debug:
            print("Process finished")
        
class ImageCombiner():
    #This class contains a function to combine images
    def __init__(self):
        pass
        
    def stack_images(self, images, calculation):
        #Function to combine images
        #The first argument must be a list with the paths of the images as strings,
        #the second one is a string that must be either "median" or "average", depending on the desired calculation method.
        master_image_path = images[0]
        with fits.open(master_image_path) as master:
            master_image_data = master[0].data
        
        if calculation == "median":
            l = []
            for i in images:
                with fits.open(i) as image:
                    l.append(image[0].data.astype("uint"))
                
            for y in range(4096):
                for x in range(4096):
                    m = []
                    for i in l:
                        m.append(i[x][y])
                    master_image_data[x][y] = statistics.median(m)

            if debug:
                print(master_image_data)
                
            with fits.open(master_image_path) as master:
                master[0].data = master_image_data
                master.writeto("masterfile.fts", overwrite=True)
                       
        elif calculation == "average":
            l = []
            for i in images:
                with fits.open(i) as image:
                    data = image[0].data
                    l.append(data)
            
            for i in range(len(images) - 1):
                l[i] = np.add(l[i], l[i + 1])
        
            stacked_images = np.true_divide(l[i], len(images))
            if debug:
                print(stacked_images)
            
#runs the multiprocessing
            with fits.open(master_image_path) as master:
                master[0].data = master_image_data
                master.writeto("masterfile.fts", overwrite=True)
                
if __name__ == '__main__':
    mp.set_start_method('fork') #Maybe try both methods 'fork' and 'spawn', it can lead to the program being faster or just working at all.
    #Depending on your operating system, only one of the methods may work.
    if debug:
        print("start")
        print(str(progress) + " %")
    
    image_path = '2021-04-02T19-10-26_M1_Clear_280s_Simon-H.fts' #relative path of the main image to be edited
    #If only the name of the file is given, the program and the file have to be in the same folder.
    #The path can also be given by the files absolute path: 'folder/folder/folder/.../file', also for Windows use normal slashes
    
    if calibrate_with_dark:
        dark_path = 'HAT-P-10-001dark.fit' #The dark image's path as a string
    else:
        dark_path = image_path
    if calibrate_with_bias:
        bias_path = 'Bias-001.fit' #The bias image's path as a string
    else:
        bias_path = image_path
    if calibrate_with_flat:
        flat_path = 'HAT-P-10-001light.fit' #The light image's path as a string
    else:
        flat_path = image_path

    calibPip = CalibrationPipe(image_path, dark_path, bias_path, flat_path)
    
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

    if debug:
        #print('stashing: ' + str(time.time() - now))
        print('total: ' + str(time.time() - t))
        print("stop")