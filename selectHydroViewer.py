########## HYDROMETEOR VIEWER ##########
# IMPORT PACKAGES
import sys
import numpy as np
import xarray as xr
from PIL import Image
#import scipy.misc as smp
#
# IMAGE BUFFER & PLOTTING FUNCTIONS
def probe_defaults(probeName):
    invalidSlice = np.array([-1., -1., -1., -1., -1., -1., -1., -1.])
    if probeName=='2DS' or probeName=='HVPS':
        boundary = np.array([43690, 43690, 43690, 43690, 43690, 43690, 43690, 43690])
        boundaryTime = 0;
        bufferShape = [1700, 128]
    elif probeName=='CIP' or probeName=='PIP':
        boundary = np.array([170, 170, 170, 170, 170, 170, 170, 170])
        boundaryTime = 0;
        bufferShape = [1700, 64]
        
    return boundary, boundaryTime, invalidSlice, bufferShape

def load_partData(inFile, probeName):
    ds1 = xr.open_dataset(inFile)
    # check if rectangular & elliptical fits are part of particle data
    iRecEll = 0
    for varname, da in ds1.data_vars.items():
        if varname=='image_RectangleL':
            iRecEll = 1

    time = (ds1['Time'].values).astype(int)
    frame = (ds1['parent_rec_num'].values).astype(int)
    partNum = (ds1['particle_num'].values).astype(int)
    length = (ds1['image_length'].values).astype(int)
    longestY = (ds1['image_longest_y'].values).astype(int)
    width = (ds1['image_width'].values).astype(int)
    dmax = ds1['image_diam_minR'].values
    dmax[np.isnan(dmax)] = -1 # set dummy value for particle size if NaN
    if iRecEll==1:
        drec = ds1['image_RectangleL'].values
        dell = ds1['image_EllipseL'].values
    else:
        drec = np.empty(len(time))*np.nan
        dell = np.empty(len(time))*np.nan
    area = ds1['image_area'].values
    perim = ds1['image_perimeter'].values
    reject = (ds1['image_auto_reject'].values).astype(int)
    habit = (ds1['holroyd_habit'].values).astype(int)
    if probeName=='2DS' or probeName=='HVPS':
        tempTime = ds1['Time_in_seconds'] # time in TAS clock cycles
        intArr = np.zeros(len(tempTime.values))
        intArr[1:] = np.diff(tempTime.values) # time difference between particles
    else:
        intArr = ds1['inter_arrival'].values
    
    return time, frame, partNum, length, longestY, width, dmax, drec, dell, area, perim, reject, habit, intArr

def get_partInds(time, dmax, reject, habit, intArr, startTime, endTime, minD, maxD, intArrThresh, rejStatus, habitStatus):
    if rejStatus=='None':
        if habitStatus=='None':
            partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                (intArr>=intArrThresh))[0]
        else:
            if isinstance(habitStatus, int):
                habitStatus = [habitStatus]

            if len(habitStatus)==1:
                partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                    (intArr>=intArrThresh) & (habit==habitStatus[0]))[0]
            elif len(habitStatus)==2:
                partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                    (intArr>=intArrThresh) & ((habit==habitStatus[0]) | (habit==habitStatus[1])))[0]
            elif len(habitStatus)==3:
                partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                    (intArr>=intArrThresh) & ((habit==habitStatus[0]) | (habit==habitStatus[1]) |
                                                             (habit==habitStatus[2])))[0]
            elif len(habitStatus)==4:
                partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                    (intArr>=intArrThresh) & ((habit==habitStatus[0]) | (habit==habitStatus[1]) |
                                                             (habit==habitStatus[2]) | (habit==habitStatus[3])))[0]
    else:
        if isinstance(rejStatus, int):
                rejStatus = [rejStatus]

        if habitStatus=='None':
            if len(rejStatus)==1:
                partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                    (intArr>=intArrThresh) & (reject==rejStatus[0]))[0]
            elif len(rejStatus)==2:
                partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                    (intArr>=intArrThresh) & ((reject==rejStatus[0]) | (reject==rejStatus[1])))[0]
            elif len(rejStatus)==3:
                partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                    (intArr>=intArrThresh) & ((reject==rejStatus[0]) | (reject==rejStatus[1]) |
                                                             (reject==rejStatus[2])))[0]
            elif len(rejStatus)==4:
                partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                    (intArr>=intArrThresh) & ((reject==rejStatus[0]) | (reject==rejStatus[1]) |
                                                             (reject==rejStatus[2]) | (reject==rejStatus[3])))[0]
        else:
            if isinstance(habitStatus, int):
                habitStatus = [habitStatus]

            if len(habitStatus)==1:
                if len(rejStatus)==1:
                    partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                        (intArr>=intArrThresh) & (habit==habitStatus[0]) & (reject==rejStatus[0]))[0]
                elif len(rejStatus)==2:
                    partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                        (intArr>=intArrThresh) & (habit==habitStatus[0]) &
                                        ((reject==rejStatus[0]) | (reject==rejStatus[1])))[0]
                elif len(rejStatus)==3:
                    partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                        (intArr>=intArrThresh) & (habit==habitStatus[0]) &
                                        ((reject==rejStatus[0]) | (reject==rejStatus[1]) | (reject==rejStatus[2])))[0]
                elif len(rejStatus)==4:
                    partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                        (intArr>=intArrThresh) & (habit==habitStatus[0]) &
                                        ((reject==rejStatus[0]) | (reject==rejStatus[1]) | (reject==rejStatus[2]) |
                                         (reject==rejStatus[3])))[0]
            elif len(habitStatus)==2:
                if len(rejStatus)==1:
                    partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                        (intArr>=intArrThresh) & ((habit==habitStatus[0]) | (habit==habitStatus[1])) &
                                        (reject==rejStatus[0]))[0]
                elif len(rejStatus)==2:
                    partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                        (intArr>=intArrThresh) & ((habit==habitStatus[0]) | (habit==habitStatus[1])) &
                                        ((reject==rejStatus[0]) | (reject==rejStatus[1])))[0]
                elif len(rejStatus)==3:
                    partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                        (intArr>=intArrThresh) & ((habit==habitStatus[0]) | (habit==habitStatus[1])) &
                                        ((reject==rejStatus[0]) | (reject==rejStatus[1]) | (reject==rejStatus[2])))[0]
                elif len(rejStatus)==4:
                    partInds = np.where((time>=startTime) & (time<endTime) & (dmax>=minD) & (dmax<maxD) &
                                        (intArr>=intArrThresh) & ((habit==habitStatus[0]) | (habit==habitStatus[1])) &
                                        ((reject==rejStatus[0]) | (reject==rejStatus[1]) | (reject==rejStatus[2]) |
                                         (reject==rejStatus[3])))[0]
                    
    return partInds

def get_imageData(inFile, frameStart):
    ds = xr.open_dataset(inFile)
    data = ds['data'][frameStart-1:frameStart][0].values
    
    return data

def get_slice_endpoints(probeName, buf, boundary, boundaryTime, invalidSlice): # get slice start/end indices, and also for the boundaries
    numPart = 0 # number of particles in buffer
    startInd = [] # index of the start of a particle
    endInd = [] # index of the end of a particle
    boundaryInd = [] # index of the particle boundary (series of 8 consecutive '43690' values)
    
    j = 0
    while (buf[j,0] != -1) and (j+1 < buf.shape[0]):
        if (np.array_equal(buf[j,:],boundary)) and ((buf[j+1,0]==boundaryTime) or (probeName=='CIP') or
                                                   (probeName=='PIP')): # particle boundary
            boundaryInd.append(j)
            if j>0:
                endInd.append(j-1) # index of particle end before the particle boundary (current particle)
                numPart = numPart + 1 # for particle preceeding boundary
            startInd.append(j+2) # index of particle start after the particle boundary (next particle)
        j = j + 1
    
    boundaryInd = np.array(boundaryInd, dtype='int')
    if (boundaryInd.size>0) or (boundaryInd[0]>0): # no boundaries or first boundary after first slice in buffer
        startInd = np.insert(startInd, 0, 0)
    startInd = np.array(startInd, dtype='int')
    
    if boundaryInd[-1]+3 < buf.shape[0]: # particle occurs after last boundary if it's 3rd from last slice in record
        endInd.append(buf.shape[0]-1) # last particle ends at end of buffer
        endInd = np.array(endInd, dtype='int') 
    else:
        endInd = np.array(endInd, dtype='int') # boundary found at very end of record -- no more particles in buffer
        
    return numPart, boundaryInd, startInd, endInd
        
def image_buffer(buf, probeName, boundaryInd): # generate matrix of 1's and 0's from buffer
    if probeName=='2DS' or probeName=='HVPS':
        #boundaryData = np.tile([2,2,1,1], 32) # alternate 1's and 2's for boundary slice (white & cyan pixels)
        boundaryData = np.tile([1,1,1,1], 32) # array of white pixels for boundary slice
        buf[buf==-1] = 0 # change invalid values to 0 (unshadowed segment)
        #buf[buf==-1] = 255 # change invalid values to 255 (unshadowed segment)
        buf = 65535 - buf # 0: shadowed; 1: unshadowed
        
        # convert decimal to binary (8 image words for each slice)
        imageData = np.ones([1700,128]) # set up image buffer (1 means unshadowed pixels - RGB(255,255,255) is white

        for x in np.arange(buf.shape[0]):
            tempBuf = np.array([np.binary_repr(int(buf[x,0]),16), np.binary_repr(int(buf[x,1]),16),
                                np.binary_repr(int(buf[x,2]),16), np.binary_repr(int(buf[x,3]),16),
                                np.binary_repr(int(buf[x,4]),16), np.binary_repr(int(buf[x,5]),16),
                                np.binary_repr(int(buf[x,6]),16), np.binary_repr(int(buf[x,7]),16)])
            sliceBuf = []
            for y in np.arange((buf.shape[1])*16):
                sliceBuf.append(tempBuf[(np.floor(y/16)).astype(int)][np.mod(y,16)])
            sliceBuf = np.asarray(sliceBuf, dtype='int')
            imageData[x,:] = sliceBuf
    elif probeName=='CIP' or probeName=='PIP':
        #boundaryData = np.tile([2,2,1,1], 16) # alternate 1's and 2's for boundary slice (white & cyan pixels)
        boundaryData = np.tile([1,1,1,1], 16) # array of white pixels for boundary slice
        #buf[buf==-1] = 0 # change invalid values to 0 (unshadowed segment)
        buf[buf==-1] = 255 # change invalid values to 1 (unshadowed segment)
        
        # convert decimal to binary (8 image words for each slice)
        imageData = np.ones([1700,64]) # set up image buffer (1 means unshadowed pixels - RGB(255,255,255) is white

        for x in np.arange(buf.shape[0]):
            tempBuf = np.array([np.binary_repr(int(buf[x,0]),8), np.binary_repr(int(buf[x,1]),8),
                                np.binary_repr(int(buf[x,2]),8), np.binary_repr(int(buf[x,3]),8),
                                np.binary_repr(int(buf[x,4]),8), np.binary_repr(int(buf[x,5]),8),
                                np.binary_repr(int(buf[x,6]),8), np.binary_repr(int(buf[x,7]),8)])
            sliceBuf = []
            for y in np.arange((buf.shape[1])*8):
                sliceBuf.append(tempBuf[(np.floor(y/8)).astype(int)][np.mod(y,8)])
            sliceBuf = np.asarray(sliceBuf, dtype='int')
            imageData[x,:] = sliceBuf

    imageData[boundaryInd,:] = boundaryData # write in boundary slice
    
    if boundaryInd[-1]+1 < buf.shape[0]:
        imageData[boundaryInd+1,:] = 1
    else:
        imageData[boundaryInd[0:-2]+1,:] = 1
    
    return(imageData)

# ===== MAIN SCRIPT ===== #
'''
Main script that calls subroutines to plot particles meeting user-specified criteria.

Parameters
----------
imageFile: str
    Path to decompressed data generated from read_binary_* script in UIOPS.
particleFile: str
    Path to particle-by-particle data generated from imgProc_sm script in UIOPS.
plotDirectory: str
    Path to save image strips to file.
campaign: str
    Name of project (e.g., 'olympex'). Allows for project-specific conditional statements to be added to plotting routines.
probeName: str
    Probe type ('2DS', 'HVPS', 'CIP', 'PIP'). Used in determining image decryption methods specific to the manufacturer.
startTime: str
    Flight time to begin the plotting job in HHMMSS.
endTime: str
    Flight time to end the plotting job in HHMMSS.
rejStatus: float tuple
    Array of rejection status values (from image_auto_reject variable) of which to plot particles.
        48 - accepted; 97 - aspect ratio > 6; 116 - aspect ratio > 5 + image touching edge;
        112 - < 25% shadowed diodes in rectangle encompassing particle; 104, 72, 117 - hollow particle;
        115 - split image; 122 - zero area image; 102 - zero area image
minD: float
    Minimum particle size (mm) of which to plot particles.
maxD: float
    Minimum particle size (mm) of which to plot particles.
intArrThresh: float
    Minimum inter-arrival time (s) of which to plot particles.
habitStatus: float tuple
    Array of habit type to use in particle plotting.
        77 - zero image; 67 - center-out image; 116 - tiny; 111- oriented; 108 - linear; 97 - aggregate;
        103 - graupel; 115 - sphere; 104 - hexagonal; 105 - irregular; 100 - dendrite
'''
imageFile = sys.argv[1]
particleFile = sys.argv[2]
plotDirectory = sys.argv[3]
campaign = sys.argv[4]
probeName = sys.argv[5]
startTime = sys.argv[6]
endTime = sys.argv[7]
rejStatus = sys.argv[8]
minD = sys.argv[9]
maxD = sys.argv[10]
intArrThresh = sys.argv[11]
habitStatus = sys.argv[12]

#map = np.empty((3,3),dtype=np.uint8)
#map[0,:] = [0, 0, 0]; map[1,:] = [255, 255, 255]; map[2,:] = [0, 255, 255] # [black, white, cyan]
[boundary, boundaryTime, invalidSlice, bufferShape] = probe_defaults(probeName)

[time, frame, partNum, length, longestY, width, dmax, drec, dell, area, perim, reject, habit, intArr] = load_partData(
    particleFile, probeName) # load Particle Data

# ---------- Initialize Particle Criteria ----------
if startTime=='None':
    startTime = int(time[0])
else:
    startTime = int(startTime)

if endTime=='None':
    endTime = int(time[-1])
else:
    endTime = int(endTime)

if rejStatus!='None':
    rejStatus = np.array(rejStatus.split(','), dtype=np.uint8)

if minD=='None':
    minD = 0.
else:
    minD = np.float(minD)

if maxD=='None':
    maxD = 999.
else:
    maxD = np.float(maxD)

if intArrThresh=='None':
    intArrThresh = -999.
else:
    intArrThresh = np.float(intArrThresh)
    
if habitStatus!='None':
    habitStatus = np.array(habitStatus.split(','), dtype=np.uint8)

print('----- USER DEFINED VALUES -----')
print('startTime = {}; endTime = {}\nminD = {} mm; maxD = {} mm\nintArrThresh = {} s; rejStatus = {}; habitStatus = {}'
      .format(startTime, endTime, minD, maxD, intArrThresh, rejStatus, habitStatus))
print('-------------------------------')
# --------------------------------------------------

# particle indices meeting user criteria
partInds = get_partInds(time, dmax, reject, habit, intArr, startTime, endTime, minD, maxD, intArrThresh, rejStatus,
                        habitStatus)
frameInds = frame[partInds] # frame number of particles meeting user criteria
partnumInds = partNum[partInds]
timeInds = time[partInds]
print('Found {} particles that match your criteria.'.format(len(partInds)))
print('Constructing image buffers. Estimating ~ {} files will be generated.'.format(
    int(np.ceil((np.sum(length[partInds])+len(partInds))/1700))))

uniqueFrames = np.unique(frameInds)

imgPointer = 0
imgNum = 1
imgNew = np.ones(bufferShape) # set up image buffer (1's mean unshadowed pixels)

#
# LOOP THROUGH FRAMES FOR PLOTTING
for iter in range(len(uniqueFrames)):
    partSubinds = partInds[frameInds==uniqueFrames[iter]] # good particles for current frame in plotting loop
    partnumSubinds = partnumInds[frameInds==uniqueFrames[iter]]
    timeSubinds = timeInds[frameInds==uniqueFrames[iter]]
    data = get_imageData(imageFile, uniqueFrames[iter])
    [partCount, boundaryInd, partStart, partEnd] = get_slice_endpoints(probeName, data, boundary, boundaryTime, invalidSlice)
    img = image_buffer(data, probeName, boundaryInd)

    for particles in range(len(partSubinds)):
        #img_sub = np.array(img[partStart[partnumSubinds[particles]-1]:partEnd[partnumSubinds[particles]-1]+2,:],dtype=np.uint8)
        img_sub = np.array(255*img[partStart[partnumSubinds[particles]-1]:partEnd[partnumSubinds[particles]-1]+2,:],
                           dtype=np.uint8)

        if imgPointer+img_sub.shape[0]<bufferShape[0]: # current particle fits within buffer
            if imgPointer==0: # indicates we're starting a new image buffer to plot
                imgTimeStart = timeSubinds[particles]
            imgNew[imgPointer:imgPointer+img_sub.shape[0],:] = img_sub
            imgPointer = imgPointer + img_sub.shape[0];
        else: # current particle DOES NOT fit within buffer, save current one and begin new one
            imgTimeEnd = timeSubinds[particles]
            #image = smp.toimage(imgNew.T, cmin=0, cmax=3, pal=map, mode='P')
            fileStr = '{}{}fr{:02d}.{}_{}.png'.format(plotDirectory, probeName, imgNum, imgTimeStart, imgTimeEnd)
            Image.fromarray(imgNew.T.astype(np.uint8)).save(fileStr) # using PIL instead of scipy.misc.toimage
            #image.save(fileStr)
            imgPointer = 0; imgNum = imgNum + 1; imgNew = np.ones(bufferShape)

            imgTimeStart = timeSubinds[particles]
            imgNew[imgPointer:imgPointer+img_sub.shape[0],:] = img_sub
            imgPointer = imgPointer + img_sub.shape[0]

        if np.remainder(imgNum+1,1000)==0:
            print('Finished plotting buffer # {}.'.format(imgNum))