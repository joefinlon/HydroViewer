########## HYDROMETEOR VIEWER ##########
# IMPORT PACKAGES
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
#mpl.use('PDF')
from matplotlib import colors 
import matplotlib.pyplot as plt
import xarray as xr
import glob
#
# IMAGE BUFFER & PLOTTING FUNCTIONS
def get_slice_endpoints(buf): # get particle slice start/end indices, and the indices for the particle boundary
    boundary = np.array([43690, 43690, 43690, 43690, 43690, 43690, 43690, 43690])
    boundaryTime = 0;
    invalidSlice = np.array([-1., -1., -1., -1., -1., -1., -1., -1.])
    numPart = 0 # number of particles in buffer
    startInd = [] # index of the start of a particle
    endInd = [] # index of the end of a particle
    boundaryInd = [] # index of the particle boundary (series of 8 consecutive '43690' values)
    
    j = 0
    while (buf[j,0] != -1) and (j+1 < buf.shape[0]):
        if (np.array_equal(buf[j,:],boundary)) and (buf[j+1,0]==boundaryTime): # particle boundary
            boundaryInd.append(j)
            if j>0:
                endInd.append(j-1) # index of particle end before the particle boundary (previous particle)
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
    
def buffer_integrity(partCount, boundaryInd, partStart, partEnd):
    print('There are {} particles and {} boundaries.'.format(partCount, len(boundaryInd)))
    if np.any(partEnd-partStart<0): # particle end indices NOT >= start indices [ERROR]
        print('Error with particle start/end indices!')
    else:
        print('Particle start/end indices look OK.')
        
def image_buffer(buf, boundaryInd): # generate matrix of 1's and 0's from buffer
#     boundaryData = np.tile([2,1,1,1], 32) # alternate 1's and 2's for boundary slice (white & cyan pixels)
    boundaryData = np.tile([2,2,1,1], 32) # alternate 1's and 2's for boundary slice (white & cyan pixels)
    buf[buf==-1] = 0 # change invalid values to 0 (unshadowed segment)
    buf = 65535 - buf # 0: shadowed; 1: unshadowed
    
    # convert decimal to binary (8 image words for each slice)
    imageData = np.ones([1700,128]) # set up image buffer (1's mean unshadowed pixels)

    for x in np.arange(buf.shape[0]):
        tempBuf = np.array([np.binary_repr(buf[x,0],16), np.binary_repr(buf[x,1],16), np.binary_repr(buf[x,2],16),
                            np.binary_repr(buf[x,3],16), np.binary_repr(buf[x,4],16), np.binary_repr(buf[x,5],16),
                            np.binary_repr(buf[x,6],16), np.binary_repr(buf[x,7],16)])
        sliceBuf = []
        for y in np.arange((buf.shape[1])*16):
            sliceBuf.append(tempBuf[(np.floor(y/16)).astype(int)][np.mod(y,16)])
        sliceBuf = np.asarray(sliceBuf, dtype='int')
        imageData[x,:] = sliceBuf
    
    imageData[boundaryInd,:] = boundaryData
    if boundaryInd[-1]+1 < buf.shape[0]:
        imageData[boundaryInd+1,:] = 1.
    else:
        imageData[boundaryInd[0:-2]+1,:] = 1.
    
    return(imageData)

def annotate_particle_incriments(boundaryInd): # get indices along buffer for every 5 particles (optional plotting)
    partInds = boundaryInd[np.arange(3,(boundaryInd.size)-1, 5)]+1
    
    return partInds

#
#
#
# USER INPUTS
# User inputs go here. A few notes:
# campaign: 'plows', 'mc3e', 'gcpex', 'pecan', 'olympex', etc.
# date: 'yyyymmdd' (for file output path)
# probeName: '2DC', '2DP', 'CIP', 'PIP', '2DS', 'HVPS'
# imageMode: 0 - all particles shaded black; 1 - rejected particles shaded blue; 2 - particles shaded by reject code; 3 - particles shaded by habit
# annotate: 1 - a red asterisk is overlaid every 5 particles (useful for matching w/ PBP data); 0 - no annotations
# inFile: file path to decompressed image file (ouput data from readbinary*.m script)
def initialize_inputs(campaign, date, probeName, imageMode, annotate, chunkNum, inFile):
    cmap1 = colors.ListedColormap(['black', 'white', 'cyan'])
    bounds=[0,1,2,3]
    norm = colors.BoundaryNorm(bounds, cmap1.N)
    #
    ds = xr.open_dataset(inFile)
    numFrames = ds.dims['time']
    if 60000*chunkNum > numFrames:
        frameStart = np.arange(60000*(chunkNum-1)+1,numFrames-4,6)
    else:
        frameStart = np.arange(60000*(chunkNum-1)+1,60000*chunkNum-4,6)
    numFrames = 6*len(frameStart)
    print('Now starting iteration #{}.'.format(chunkNum))
    #
    # LOOP THROUGH FRAMES FOR PLOTTING
    for iter in range(len(frameStart)):
        yr = ds['year'][frameStart[iter]:frameStart[iter]+6]
        mon = ds['month'][frameStart[iter]:frameStart[iter]+6]
        day = ds['day'][frameStart[iter]:frameStart[iter]+6]
        hr = ds['hour'][frameStart[iter]:frameStart[iter]+6]
        minute = ds['minute'][frameStart[iter]:frameStart[iter]+6]
        sec = ds['second'][frameStart[iter]:frameStart[iter]+6]
        msec = ds['millisec'][frameStart[iter]:frameStart[iter]+6]
        data = ds['data'][frameStart[iter]:frameStart[iter]+6]
        data = data.values
    
        boundary = np.array([43690, 43690, 43690, 43690, 43690, 43690, 43690, 43690])
        boundaryTime = 0;
        invalidSlice = np.array([-1., -1., -1., -1., -1., -1., -1., -1.])
    
        [partCount1, boundaryInd, partStart, partEnd] = get_slice_endpoints(data[0])
        img1 = image_buffer(data[0], boundaryInd)
        if annotate==1:
            partInds1 = annotate_particle_incriments(boundaryInd)

        [partCount2, boundaryInd, partStart, partEnd] = get_slice_endpoints(data[1])
        img2 = image_buffer(data[1], boundaryInd)
        if annotate==1:
            partInds2 = annotate_particle_incriments(boundaryInd)

        [partCount3, boundaryInd, partStart, partEnd] = get_slice_endpoints(data[2])
        img3 = image_buffer(data[2], boundaryInd)
        if annotate==1:
            partInds3 = annotate_particle_incriments(boundaryInd)

        [partCount4, boundaryInd, partStart, partEnd] = get_slice_endpoints(data[3])
        img4 = image_buffer(data[3], boundaryInd)
        if annotate==1:
            partInds4 = annotate_particle_incriments(boundaryInd)

        [partCount5, boundaryInd, partStart, partEnd] = get_slice_endpoints(data[4])
        img5 = image_buffer(data[4], boundaryInd)
        if annotate==1:
            partInds5 = annotate_particle_incriments(boundaryInd)

        [partCount6, boundaryInd, partStart, partEnd] = get_slice_endpoints(data[5])
        img6 = image_buffer(data[5], boundaryInd)
        if annotate==1:
            partInds6 = annotate_particle_incriments(boundaryInd)
        
        fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, sharex=True, sharey=True, figsize=(16,8))

        ax1.imshow((img1[0:1701,:]).T, cmap=cmap1, norm=norm, aspect='auto')
        #ax1.pcolor((img1[0:1701,:]).T, cmap=cmap1, norm=norm)
        #ax1.set_aspect('auto')
        fileStr = '{:02d}{:02d}{:02d}_fr{}'.format(hr.values[0], minute.values[0], sec.values[0], frameStart[iter])
        titleStr = 'Frame #{}: {:02d}{:02d}{:02d}.{:03d} | {} particles'.format(frameStart[iter], hr.values[0],
                                                                                minute.values[0],sec.values[0],
                                                                                msec.values[0], partCount1)
        if annotate==1:
            for x in range(len(partInds1)):
                ax1.text(partInds1[x], 3, '*', color='red', horizontalalignment='left', verticalalignment='top')
        ax1.set_title(titleStr, size=10)
        ax1.set_xticks([]), ax1.set_yticks([])

        ax2.imshow((img2[0:1701,:]).T, cmap=cmap1, norm=norm, aspect='auto')
        #ax2.pcolor((img2[0:1701,:]).T, cmap=cmap1, norm=norm)
        #ax2.set_aspect('auto')
        titleStr = 'Frame #{}: {:02d}{:02d}{:02d}.{:03d} | {} particles'.format(frameStart[iter]+1, hr.values[1],
                                                                                minute.values[1],sec.values[1],
                                                                                msec.values[1], partCount2)
        if annotate==1:
            for x in range(len(partInds2)):
                ax2.text(partInds2[x], 3, '*', color='red', horizontalalignment='left', verticalalignment='top')
        ax2.set_title(titleStr, size=10)
        ax2.set_xticks([]), ax2.set_yticks([])

        ax3.imshow((img3[0:1701,:]).T, cmap=cmap1, norm=norm, aspect='auto')
        #ax3.pcolor((img3[0:1701,:]).T, cmap=cmap1, norm=norm)
        #ax3.set_aspect('auto')
        titleStr = 'Frame #{}: {:02d}{:02d}{:02d}.{:03d} | {} particles'.format(frameStart[iter]+2, hr.values[2],
                                                                                minute.values[2],sec.values[2],
                                                                                msec.values[2], partCount3)
        if annotate==1:
            for x in range(len(partInds3)):
                ax3.text(partInds3[x], 3, '*', color='red', horizontalalignment='left', verticalalignment='top')
        ax3.set_title(titleStr, size=10)
        ax3.set_xticks([]), ax3.set_yticks([])

        ax4.imshow((img4[0:1701,:]).T, cmap=cmap1, norm=norm, aspect='auto')
        #ax4.pcolor((img4[0:1701,:]).T, cmap=cmap1, norm=norm)
        #ax4.set_aspect('auto')
        titleStr = 'Frame #{}: {:02d}{:02d}{:02d}.{:03d} | {} particles'.format(frameStart[iter]+3, hr.values[3],
                                                                                minute.values[3],sec.values[3],
                                                                                msec.values[3], partCount4)
        if annotate==1:
            for x in range(len(partInds4)):
                ax4.text(partInds4[x], 3, '*', color='red', horizontalalignment='left', verticalalignment='top')
        ax4.set_title(titleStr, size=10)
        ax4.set_xticks([]), ax4.set_yticks([])

        ax5.imshow((img5[0:1701,:]).T, cmap=cmap1, norm=norm, aspect='auto')
        #ax5.pcolor((img5[0:1701,:]).T, cmap=cmap1, norm=norm)
        #ax5.set_aspect('auto')
        titleStr = 'Frame #{}: {:02d}{:02d}{:02d}.{:03d} | {} particles'.format(frameStart[iter]+4, hr.values[4],
                                                                                minute.values[4],sec.values[4],
                                                                                msec.values[4], partCount5)
        if annotate==1:
            for x in range(len(partInds5)):
                ax5.text(partInds5[x], 3, '*', color='red', horizontalalignment='left', verticalalignment='top')
        ax5.set_title(titleStr, size=10)
        ax5.set_xticks([]), ax5.set_yticks([])

        ax6.imshow((img6[0:1701,:]).T, cmap=cmap1, norm=norm, aspect='auto')
        #ax6.pcolor((img6[0:1701,:]).T, cmap=cmap1, norm=norm)
        #ax6.set_aspect('auto')
        titleStr = 'Frame #{}: {:02d}{:02d}{:02d}.{:03d} | {} particles'.format(frameStart[iter]+5, hr.values[5],
                                                                                minute.values[5],sec.values[5],
                                                                                msec.values[5], partCount6)
        if annotate==1:
            for x in range(len(partInds6)):
                ax6.text(partInds6[x], 3, '*', color='red', horizontalalignment='left', verticalalignment='top')
        ax6.set_title(titleStr, size=10)
        ax6.set_xticks([]), ax6.set_yticks([])
    
        if np.remainder(iter+1,500)==0:
            print('Finished plotting frame # {} / {}'.format(iter+1, numFrames/6))
        
        outFile = '/data/gpm/a/shared/finlon2/{}/images/{}/{}_imgmod{}_annot{}/{}.{}.{}.png'.format(campaign, date, probeName,
                                                                                                    imageMode,annotate, date,
                                                                                                    fileStr, probeName)
        plt.savefig(outFile, dpi=300)
        #plt.savefig(outFile, format='pdf')
        plt.clf(), plt.cla(), plt.close('all') # close figure so that it doesn't get displayed