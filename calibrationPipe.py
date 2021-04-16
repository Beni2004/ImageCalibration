#Contributors: Benjamin A. and Simon H.
#For Questions, feel free to contact imagecalibration-contributors@outlook.com
#put in your file paths starting from line 41
#If something isn't working as it's supposed to, try:
    #Have a look at #multiprocessing on lines 20-26
    #check if you have all the necessary modules correctly installed (astropy and numpy)
    #check if you declared your filepaths correctly (image_path and all the ones you set calibrate_with_[frame] to True)
#While the program is running, it only prints the progress of the hotpixel removal, not the progress of the whole process.

from astropy.io import fits #non-default module, installation required
import numpy as np #non-default module, installation required
import statistics #default module, no installation required
import time #default module, no installation required
import multiprocessing as mp #default module, no installation required

#debug
#if debug is 'True' it prints some helpful stuff to the shell, otherwise not.
debug = True #Set this to True or False

#multiprocessing
n=16 #number of processes executed at the same time
#If something doesn't work, try lower numbers like 2 or 4. 2^n numbers are favorable.

method = 'fork'
#Maybe try both methods: 'fork' and 'spawn', it can lead to the program being faster or just working at all.
#Depending on your operating system, only one of the methods may work.

#boolean inputs
#wether or not the program should remove hotpixels or not, depends on wether you need it or not (would require most of the time, ～1.5-2 min.).
#If you have a dark frame to subtract, this may not be necessary
remove_hotpixels = True #Set this to True or False

#wether or not the program should do the calibration, i.e. bias, dark, flat, depends on wether you need it and have the necessary files.
calibrate_with_dark = False #Set this to True or False
calibrate_with_bias = False #Set this to True or False
calibrate_with_flat = False #leave False normally, don't change this

#Defining files
#If only the name of the file is given, the program and the file have to be in the same folder.
#The path can also be given by the file's absolute path: 'folder/folder/folder/.../file' (also for Windows OS use normal slashes, not backslashes)

image_path = 'Master_Blue_M1.fit' #Put your (relative) path of the main image to be edited

if calibrate_with_dark:
    dark_path = 'dark.fit' #Put your dark frame's (relative) path as a string
else:
    dark_path = None

if calibrate_with_bias:
    bias_path = 'bias.fit' #Put your bias frame's (relative) path as a string
else: 
    bias_path = None

if calibrate_with_flat:
    flat_path = 'flat.fit' #Put your flat frame's (relative) path as a string
else:
    flat_path = None
    
t = time.time()
progress=0 #starting variables to see the progress and the total time needed in the end.

class CalibrationPipe():
    #opens all images with their path if they are given. If not, they are not opened.
    def __init__(self, image_path, dark_path=None, bias_path=None, flat_path=None):
        with fits.open(image_path) as image:
            self.imagedata = image[0].data
            self.imagehdr = image[0].header
        
        if dark_path != None:
            with fits.open(dark_path) as dark:
                self.darkdata = dark[0].data
                self.darkhdr = dark[0].header
                self.doDark = True
        else:
            self.doDark = False
            self.darkdata = None
            self.darkhdr = None
        
        if bias_path != None:
            with fits.open(bias_path) as bias:
                self.biasdata = bias[0].data
                self.biashdr = bias[0].header
                self.doBias = True
        else:
            self.doBias = False
            self.biasdata = None
            self.biashdr = None
        
        if flat_path != None:
            with fits.open(flat_path) as flat:
                self.flatdata = flat[0].data
                self.flathdr = flat[0].header
                self.doFlat = True
        else:
            self.doFlat = False
            self.flatdata = None
            self.flathdr = None
        
        self.imagepath = image_path
        self.average_dark_pixelvalue = None
        
    #Executes all processes, mainly the multiprocessing for scaling the dark frame and removing hotpixels    
    def run(self):
        if debug:
            print("ruuuuun!")
        
        amounty=len(self.imagedata)
        amountx=len(self.imagedata[0])
        
        # Calculates the average pixelvalue of the darkframe
        if self.doDark:
            y_averages = []
            for x in range(amountx):
                y_values = []
                for y in range(amounty):
                    y_values.append(self.darkdata[y][x])
                y_averages.append(statistics.mean(y_values))
            
            self.average_dark_pixelvalue = np.average(y_averages)
            
            #defining the processes
            if debug:
                print("Starting multiprocessing...")
            q = mp.Queue()
            processes = []
            runs = self.make_processes(self.darkdata, amountx, amounty, n, q)
            for i in runs:
                p = mp.Process(target = self.scale_dark, args = i)
                processes.append(p)
                if debug:
                    print("Iteration")
            
            #starting the processes
            for x in processes:
                x.start()
                
            if debug:
                print("All processes should have started now")
                
            #Get the resulting values from the calculation and define the new dark frame
            for i in range(n):
                tupel = q.get()
                current = tupel[0]
                pos = tupel[1]
                self.darkdata = np.concatenate((self.darkdata[0:pos*amounty//n], current[pos*amounty//n:(pos+1)*amounty//n], self.darkdata[(pos+1)*amounty//n:amounty]), axis=0)
            
            for x in processes:
                x.join()
                x.terminate()
            
            #subtract the bias from the dark frame and subtract the result from the main image
            self.imagedata = self.dark_frame(self.imagedata, self.darkdata, self.biasdata)
        
        #subtract the bias frame from the main image
        if self.doBias:
            self.imagedata = self.bias_frame(self.imagedata, self.biasdata)
            
        #subtract the flat frame from the main image
        if self.doFlat:
            self.imagedata = self.flat_frame(self.imagedata, self.flatdata)
        
        #removes the hotpixels in the whole image
        if remove_hotpixels:
            #defining the processes
            if debug:
                print("Starting multiprocessing...")
            q = mp.Queue()
            t_start = time.time()
            processes = []
            runs = self.make_processes(self.imagedata, amountx, amounty, n, q)
            for i in runs:
                p = mp.Process(target = self.clean_image, args = i)
                processes.append(p)
                if debug:
                    print("Iteration")
            
            #starting the processes
            for x in processes:
                x.start()
                
            if debug:
                print("All processes should have started now")

            #Get the resulting values from the calculation and define the new main image
            for i in range(n):
                tupel = q.get()
                current = tupel[0]
                pos = tupel[1]
                self.imagedata = np.concatenate((self.imagedata[0:pos*amounty//n], current[pos*amounty//n:(pos+1)*amounty//n], self.imagedata[(pos+1)*amounty//n:amounty]), axis=0)
            
            for x in processes:
                x.join()
                x.terminate()
            
            if debug:
                print("removal time:", time.time()-t_start)
        
        #Writes the resulting data to a new file called output.fit in the same folder as the program.
        #If there already exists a file with this name in this folder, it will be overwritten.
        with fits.open(self.imagepath) as image:
            image[0].data = self.imagedata
            image.writeto("output.fit", overwrite=True)
            
    #creates the processes for the multiprocessing
    #takes an nparray as the image, the length of the process in the x and y direction, the number of processes (an integer) and the multiprocessing-queue.
    #returns a list with all processes
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
        #the number of the process and the average pixelvalue of the darkframe
        if debug:
            print("A process was started")
        
        expt = self.imagehdr["EXPTIME"]
        #darkexpt = self.darkhdr["EXPTIME"]
        scaling_factor = (-0.000629169*(expt ** 2) + 0.147650091*expt + 1023.13674955)/self.average_dark_pixelvalue
        #This calculates the average pixelvalue that a dark with given exposuretime should have divided by the actual average pixelvalue of the dark
        #The result is the scaling factor. Each pixel then has to be multiplied with that value.
        
        for y in range(amounty):
            for x in range(amountx):
                darkdata[starty + y][startx + x] *= scaling_factor
                # darkdata[starty + y][startx + x] = 30000
        
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
        if self.doBias:
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
        if value > 1.05*av and av < 7000:
            imagedata[y, x]=av
        return imagedata
    
    #Runs the hotpixel check for every pixel in the image, except if the maximum y-coordinate is reached.
    #Many (n) processes are running simultaniously due to multiprocessing. Hence the y-axis is divided n-times.
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
    #This class contains a function to combine images, currently not in use.
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
            
            with fits.open(master_image_path) as master:
                master[0].data = master_image_data
                master.writeto("masterfile.fts", overwrite=True)

#runs the multiprocessing                
if __name__ == '__main__':
    mp.set_start_method(method)
    if debug:
        print("start")
        print(str(progress) + " %")

    calibPip = CalibrationPipe(image_path, dark_path, bias_path, flat_path)
    
    calibPip.run()

    if debug:
        #print('stashing: ' + str(time.time() - now))
        print('total: ' + str(time.time() - t))
        print("stop")