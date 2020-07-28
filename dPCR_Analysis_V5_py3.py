# Christian D. Ahrberg - BNTL Sogang University - 2017
#
# Program opening two different image files. The first image file contains a reference dye for identifying wells
# The second image file contains the flourescent images. In the flourescent images the well intensities
# are measured at the positions identified from the first image. Afterwards concentrations are callculated
# using Poisson statistics.
#
# Changelog:
# ===================
# 2017/10/27 Normalisation to reference channel added
# 2017/11/20 Removed values with normalised flu above 1
# 2017/11/20 Outputs file with callculated concentrations and threshold
# 2018/10/17 Variables modification
# 2019/10/10 Output to SQL file
#
###########################################################################################

# Importing required packages

import numpy as np
import cv2
import plotly.graph_objs as go
import plotly as py
import math
from operator import truediv
import sqlite3

# ==============================================================================================================
# Defining Variables can be changed by user

Nimg        = 1    # Number of images
wellheight  = 20    # Heigth of microwells in um
wellradius  = 10    # Radius of microwells in um

# Variables for circle detection
dp1         = 1     # Ratio of accumulator, one means same resolution as input image, bigger numbers mean image is reduced 
minDist1    = 8    # Minimum distance between two circle centers in pixels
param11     = 8     # Threshold passed to Canny edge detector
param21     = 3     # Accumulatoir threshold for circles, the lower the more false circles are recognised
minRadius1  = 1     # Minimum radius of found circles in pixcels
maxRadius1  = 2     # Maximum radius of found circles in pixcels
spec        = 4     # Which histogram to be plotted 4 = green / 5 = red / 6 = normalized    

# Defining appendixes required for opening images
name1 = "R"         # Appendix for reference dye images
name2 = "G"         # Appendix for fluorescence images
name3 = ".jpg"      # Appendix for image filetype
tag   = "-marked"   # Appendix for images in which found circles are marked

# ===============================================================================================================
# Defining constants that should not be changed by the user

Npos = 0.0 # Counter for positive wells, float so probability can be calculated
Nneg = 0.0 # Counter for negative wells, float so probability can be calculated

# ===============================================================================================================

# Function for Callculating flourescence intensity inside of cicle

def CircleIntensity(centerx,centery,radius,image, color):
    # Callculates the intensity indside of a circle
    #
    # Parameters:
    # centerx = x-coordinte of circle
    # centery = y-coordinate of circle
    # radius = Radius of circle
    # image = image for callculating intensity (brightness saved as 8 bit, i.e. from 0 to 255
    #
    # Returns:
    # Intensity = average intensity value of circle in the tested image

    if color == "B":
        coval = 0
    elif color == "G":
        coval = 1
    elif color == "R":
        coval = 2
    else:
        print("Color not RGB, assuming color is green")
        coval = 1
    
    # Definging required parameters
    npixels = 0.0     # Count for pixels to find average
    brightness = 0.0  # Count for brightness of pixel

    # Creating square around circle
    #if (centerx >=2 and centery >=2 and centerx <= image.shape[0]-1
    for x in range(centerx-radius-2,centerx+radius+2):        # Varying through x
            for y in range(centery-radius-2,centery+radius+2):       # Varying through y
                 if (x <= image.shape[1]-1 and y <= image.shape[0]-1 and x>=0 and y>=0):  # Making sure coordinate in image
                     pixeldistance = math.sqrt((centerx - x)**2 + (centery - y)**2)     # Pythagoras to find radius from iterated pixcel to center of circle
                     if pixeldistance <  radius:                                        # If Pixel is in circle add to intensity callculation
                         pixel = image[y,x]
                         brightness = brightness + float(pixel[coval])/255                  # Updating total brightness
                         npixels = npixels + 1                                          # Updating total pixcel count
    if npixels == 0:
        npixels = 1 # Preventing error, division by zero
        
    Intensity = brightness / npixels                                                    # Callculating average intesnity of circle

    if Intensity == 0:      # Preventing division by zero
        Intensity = 0.00000001 

    return Intensity

# ==============================================================================================================

# Function for plotting histograms
def Histogram(Data):

    # Making histogram
    data = [
        go.Histogram(
            x=Data, xbins=dict(start=0, size=0.001, end=1.01) 
        )
    ]

    layout = go.Layout(
        title='Histogram of flourescence intensity of wells',
        xaxis=dict(title='Flourescence intensity',
                   range=[0,1.01]),
        yaxis=dict(title='Well count'), 
        )

    fig = go.Figure(data=data, layout=layout)

    py.offline.plot(fig)

# ===============================================================================================================

# Function for esimating Concentration

def ConcCallculation(pHat,Npart,Vol):
    # Function callculating the concentration of outcome of dPCR
    #
    # Parameters:
    # pHat  =   estimated propability of positive partition
    # Npart =   total number of partitions
    # Vol   =   Volume of partitions in uL
    #
    # Returns:
    # C_est = callculated concentration in #particles/uL
    # C_low = lower confidence intervall of calculated concentration (95% z-distribution) in # particles/uL
    # C_upp = upper confidence intervall of calculated concentration (95% z-distribution) in # particles/uL
    #
    #######################################################

    # Defingin constants
    zc = 1.96   # 95% confidence intervall z-distribution

    # Callculation of confidence interval on pHat
    pHat_Dev = zc * math.sqrt((pHat * (1-pHat))/Npart)  # Deviation on expected result
    p_hat_low = pHat - pHat_Dev  # Lower bound of p_hat
    p_hat_upp = pHat + pHat_Dev  # Upper bound of p_hat

    # Callculating mean number of molecules per patition including 95%
    # confidence intervall
    lambda1 = -math.log(1-pHat)     # average number of molecules per division as per Poission distribution
    lambda_low = -math.log(1-p_hat_low)  # lower bound of average number of molecules per division
    lambda_upp = -math.log(1-p_hat_upp)  # upper bound of average number of molecules per division

    # Callculating concentrations in mol/uL from lambda values including
    # confidence intervalls
    C_est = lambda1 / Vol       # Esitmated concentration
    C_low = lambda_low / Vol    # Estimated lower bound of concentration
    C_upp = lambda_upp / Vol    # Estimated higher bound of concentration

    return C_est, C_low, C_upp

# ==============================================================================================================

# Main code file calling the previous defined functions

# Creating database for Saving results
conn = sqlite3.connect('OutputData.db')
c = conn.cursor()
c.execute('DROP TABLE IF EXISTS Wells')
c.execute('DROP TABLE IF EXISTS Concentration')
c.execute('CREATE TABLE IF NOT EXISTS Wells (ImmageNumber, Xcoordinate, Ycoordinate, Radius, GreenIntensity, RedIntensity, NormIntensity)')
c.execute('CREATE TABLE IF NOT EXISTS Concentration (Threshold, PosWells, NegWells, Propability, LowBound, UppBound, Conc)')
conn.commit()

#################################################################################################################
# Loop going through all the images taken, identifying circles, drawing cicles, and measuring intensities

for i in range(0,Nimg):
    
    # Print current image for debugging purposes
    print("Currently on image number %d" % (i+1))

    # Creating filenames for images to be opened
    nameRef = name1 + str(i+1) + name3
    nameReftag = str(i+1) + name1 + tag + name3
    nameFlu = name2 + str(i+1) + name3
    nameFlutag = str(i+1) + name2 + tag + name3
    print(nameRef)


    # Converting images for detection of circles and measurement of flourescence
    # Reference dye image
    imgRef  = cv2.imread(nameRef,0)     # Opening image
    imgRef  = cv2.medianBlur(imgRef,5)  # Smothening image data
    imgRef2 = cv2.imread(nameRef)       # Image for marking circles in
    imgRef3 = cv2.imread(nameRef)       # Measuring reference dye intesity in
    #imgRef3 = cv2.medianBlur(imgRef,5)  # Smothening image data
    # Fluorescence dye image
    imgFlu  = cv2.imread(nameFlu)       # Opening image
    #imgFlu  = cv2.medianBlur(imgFlu,5)  # Smothening image data
    imgFlu2 = cv2.imread(nameFlu)       # Image for marking circles in

    # Fitting circles using Hough transform from package cv2
    # function calls as follows: cv2.HoughCircles(image, method, dp, minDist, circles, param1, param2, minRadius, maxRadius)
    circles = cv2.HoughCircles(imgRef,cv2.HOUGH_GRADIENT,dp1,minDist1,param1=param11,param2=param21,minRadius=minRadius1,maxRadius=maxRadius1)

    # Drawing fitted circles to reference dye image
    circles = np.uint16(np.around(circles)) # Preparing data for plotting
    for j in circles[0,:]:
        # draw the outer circle
        #cv2.circle(imgRef2,(j[0],j[1]),j[2],(0,255,0),2)
        # draw the center of the circle
        cv2.circle(imgRef2,(j[0],j[1]),1,(255,255,255),1)
    # Saving image with marked circles
    cv2.imwrite(nameReftag,imgRef2)

    # Drawing fitted circles to reference dye image
    circles = np.uint16(np.around(circles)) # Preparing data for plotting
    for j in circles[0,:]:
        # draw the outer circle
        #cv2.circle(imgFlu2,(j[0],j[1]),j[2],(0,255,0),2)
        # draw the center of the circle
        cv2.circle(imgFlu2,(j[0],j[1]),1,(255,255,255),1)
    # Saving image with marked circles
    cv2.imwrite(nameFlutag,imgFlu2)

    # Measuring intensity in marked cicles
    intensity = []  # Empty list to save circle intensities into
    for j in circles[0,:]:
        intensity.append(CircleIntensity(j[0],j[1],j[2],imgFlu,name2))

    # Measuring intensity in marked cicles for reference dye
    intensity_ref = []  # Empty list to save circle intensities into
    for j in circles[0,:]:
        intensity_ref.append(CircleIntensity(j[0],j[1],j[2],imgRef3,name1))

    intensity_norm = list(map(truediv, intensity, intensity_ref))

    
 
    # Saving data to output database
    k=0     # Running index to sync intensity list with circle list
    for j in circles[0,:]:
        c.execute('''INSERT INTO Wells VALUES (?, ?, ?, ?, ?, ?, ?)''', (str(i+1) ,str(j[0]) ,str(j[1]) ,str(j[2]) ,str(intensity[k]) ,str(intensity_ref[k]) ,str(intensity_norm[k]))) 
        k = k+1

conn.commit()

#################################################################################################################
# Loop going through all the images taken, identifying circles, drawing cicles, and measuring intensities

# Importing data from database
data = c.execute('''Select * FROM Wells''')
data = [float(item[spec]) for item in data.fetchall()]

# Mapping data
plotdata = list(map(float, data))

# Plotting histogram
Histogram(plotdata)

# Defining manual threshold
Threshold = eval(input('Please enter fluorescence threshold for positive call: '))

# Counting positive and negative partitions
for x in plotdata:
    if x > 1.0:                   # Throwing out wrong counts
        Nneg = Nneg + 1
    elif (x < 1.0) & (x > Threshold):         # If well is positive average brightness is low
        Npos = Npos + 1
    elif x <= Threshold:        # Wells contianing now particle/cell have a higher intensity
        Nneg = Nneg + 1

# Calculating volume of well
VolWell = wellheight * math.pi * wellradius ** 2 # Volume of a well in cubic micrometer
VolWell = VolWell * math.pow(10,-9)

# Callculation of concentrations according to Poission distribution
pHat = Npos / (Npos + Nneg) # Estimated propability of positive well
Npart = Npos + Nneg # Total number of detected wells
C_est, C_low, C_upp = ConcCallculation(pHat,Npart,VolWell) # Callculation of concentrations according to predifined function
Stdev = C_est - C_low
CV = 100 * Stdev / C_est

print(('Number of Poitive calls:                 {0:.2f}' .format(Npos)))
print(('Number of Negative calls:                {0:.2f}' .format(Nneg)))
print(('Propability of positive partition:       {0:.2f} \n' .format(pHat)))
# Outputting results
print(('Estimated concentration:                 {0:.3f} copies/uL' .format(C_est)))
print(('Standard deviation:                      {0:.3f} copies/uL' .format(Stdev)))
print(('Coefficient of variation (CV):           {0:.3f} %' .format(CV)))
print(('Lower bound of 95% Confidence intervall: {0:.3f} copies/uL' .format(C_low)))
print(('Upper bound of 95% Confidence intervall: {0:.3f} copies/uL' .format(C_upp)))

c.execute(''' INSERT INTO Concentration VALUES (?, ?, ?, ?, ?, ?, ?)''', (Threshold , Npos, Nneg, pHat, int(C_low), int(C_upp), int(C_est)))
conn.commit()
conn.close()



    
