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

        self.rawDataPath = "./dataPath/"
        self.tmpPath = "./tmpPath/"




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




