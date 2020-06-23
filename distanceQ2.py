from myPYTHON import *
from mwispDBSCAN import MWISPDBSCAN


doFITS=myFITS()
doMWdbscan= MWISPDBSCAN()

#do this on server


class disQ2(object):

    rawDataPath ="/home/qzyan/WORK/diskMWISP/MWISPData/G100150/"
    tmpPath="/home/qzyan/WORK/diskMWISP/MWISPData/G105150Tmp/" #tmpFiles,need a large data area

    if os.path.isdir(rawDataPath):
        pass

    else: #

        rawDataPath = "./dataPath/"
        tmpPath = "./tmpPath/"



    mergeFITS= rawDataPath+"mergedCube.fits"
    rawCO12FITSLocal = rawDataPath  + "100_150_U_local.fits"
    rawCO12RMSfits = rawDataPath + "100_150_U_rms.fits"
    rawCO12FITSPer = rawDataPath  + "100_150_U_per.fits"

    cropRawCO12FITSLocal = rawDataPath  + "100_150_U_localCrop.fits"
    cropRawCO12FITSPer =  rawDataPath  + "100_150_U_PerCrop.fits"
    cropRawCO12RMSfits = rawDataPath + "100_150_U_rmsCrop.fits"


    lRange=[104.75,150.25]

    bRange=[ -5.25, 5.25 ]

    #the lRange, and bRange are used to



    def __init__(self):
        #check server environments

        pass


    def prepareData(self ):
        """
        produce cropFITS according to the lRange and bRange
        :return:
        """

        doFITS.cropFITS( self.rawCO12FITSLocal, Lrange=self.lRange,Brange=self.bRange, overWrite=True,outFITS=self.cropRawCO12FITSLocal)
        doFITS.cropFITS( self.rawCO12FITSPer, Lrange=self.lRange,Brange=self.bRange, overWrite=True,outFITS=self.cropRawCO12FITSPer )

        doFITS.cropFITS2D( self.rawCO12RMSfits, Lrange=self.lRange,Brange=self.bRange, overWrite=True,outFITS=self.cropRawCO12RMSfits )


    def processDBSCAN(self,coFITS):



        doMWdbscan.rawCOFITS =  coFITS
        doMWdbscan.rmsFITS = self.cropRawCO12RMSfits

        doMWdbscan.setDBSCANParameters( cutoff_sigma=2,minPts=4,connectivity=1)
        doMWdbscan.processPath = self.tmpPath

        doMWdbscan.computeDBSCAN()
        doMWdbscan.getCatFromLabelArray(doClean=True)

    def find_nearestIndex(self,a, a0):
        "Element in nd array `a` closest to the scalar value `a0`"



        idx = np.abs(a - a0).argmin()
        return  idx

    def mergeByVaxis(self,fits1,fits2,outPut="mergedCube.fits"):
        """
        #takes two fits files, and merge them together, to see if the SASMA can process this large data
        we merge the local and the perseus arm file,

        :param fits1:
        :param fits2:
        :param outPut:
        :return:
        """

        #find the fits, that has the lowerest velocity, and append is to the

        data1,head1=doFITS.readFITS(fits1)

        data2,head2=doFITS.readFITS(fits2)


        spec1, vaxis1 = doFITS.getSpectraByIndex(data1, head1,0,0)
        spec2, vaxis2 = doFITS.getSpectraByIndex(data2, head2,0,0)



        if vaxis1[0]< vaxis2[0]:

            lowData,  lowHead=data1,head1
            highData, highHead=data2, head2

            lowSpec,lowVaxis=spec1,vaxis1
            highSpec, highVaxis=spec2,vaxis2

        else:


            lowData,  lowHead=data2, head2
            highData, highHead=data1, head1

            lowSpec,lowVaxis=spec2,vaxis2
            highSpec, highVaxis=spec1,vaxis1


        #process low and high



        maxVlow= lowVaxis[-1]

        minVHigh = highVaxis[0]


        #find and middle position

        middleV= ( maxVlow+  minVHigh )/2.

        mergeIndexLow= self.find_nearestIndex(lowVaxis,middleV)

        mergeV=lowVaxis[mergeIndexLow]

        part1data  =  lowData[0:mergeIndexLow]

        mergeIndexHigh= self.find_nearestIndex(highVaxis,mergeV)

        part2data = highData[mergeIndexHigh:]

        if highVaxis[mergeIndexHigh] !=mergeV:
            print "Two fits has different cooridnate at velocity, cannot do direct merge, exist... "
            return

        print "Merging at {:.3f} km/s".format(mergeV)

        mergeData= np.vstack([part1data,part2data])

        fits.writeto( self.rawDataPath+outPut, mergeData , header=lowHead, overwrite=True)
