from myPYTHON import *
from mwispDBSCAN import MWISPDBSCAN
import matplotlib.patheffects as path_effects
from mwispGaia import dbscanDis
import matplotlib as mpl
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
doFITS=myFITS()
doMWdbscan= MWISPDBSCAN()




#do this on server
def box( centerL, centerB, lSize, bSize, dummy=0):
    """
    return lRange and B Range
    """

    lSize = lSize / 3600.
    bSize = bSize / 3600.

    return [centerL - lSize / 2., centerL + lSize / 2.], [centerB - bSize / 2., centerB + bSize / 2.]


class disQ2(object):

    rawDataPath ="/home/qzyan/WORK/diskMWISP/MWISPData/G100150/"
    tmpPath="/home/qzyan/WORK/diskMWISP/MWISPData/G105150Tmp/" #tmpFiles,need a large data area

    if os.path.isdir(rawDataPath):
        pass

    else: #

        rawDataPath = "./dataPath/"
        tmpPath = "./tmpPath/"

    intPath="./intPath/"



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

        try:
            os.mkdir(self.intPath)

        except:
            pass


    def prepareData(self ):
        """
        produce cropFITS according to the lRange and bRange
        :return:
        """

        doFITS.cropFITS( self.rawCO12FITSLocal, Lrange=self.lRange,Brange=self.bRange, overWrite=True,outFITS=self.cropRawCO12FITSLocal)
        doFITS.cropFITS( self.rawCO12FITSPer, Lrange=self.lRange,Brange=self.bRange, overWrite=True,outFITS=self.cropRawCO12FITSPer )

        doFITS.cropFITS2D( self.rawCO12RMSfits, Lrange=self.lRange,Brange=self.bRange, overWrite=True,outFITS=self.cropRawCO12RMSfits )


    def processDBSCAN(self,coFITS, cutoff_sigma=2,minPts=4,connectivity=1):



        doMWdbscan.rawCOFITS =  coFITS
        doMWdbscan.rmsFITS = self.cropRawCO12RMSfits

        doMWdbscan.setDBSCANParameters( cutoff_sigma=cutoff_sigma,minPts=minPts,connectivity=connectivity)
        doMWdbscan.processPath = self.tmpPath

        doMWdbscan.computeDBSCAN()
        doMWdbscan.getCatFromLabelArray(doClean=True)
        doMWdbscan.produceCleanFITS()







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



    def getCleanFITS(self ):
        """

        :param rawDBSCANLabel:
        :param cleanTB:
        :return:
        """


        doMWdbscan.labelFITSName=self.tmpPath +  "mergedCubedbscanS2P4Con1.fits"
        doMWdbscan.cleanCatName=self.tmpPath + "mergedCubedbscanS2P4Con1_Clean.fit"

        doMWdbscan.produceCleanFITS()



    def getIntFITS(self):
        """
        :return:
        """

        rawCOFITS = self.mergeFITS
        labelFITS = self.tmpPath + "mergedCubedbscanS2P4Con1_Clean.fits"
        tableFile = self.tmpPath + "mergedCubedbscanS2P4Con1_Clean.fit"

        doMWdbscan.produceCloudIntFITS(rawCOFITS,labelFITS,tableFile, outputPath = self.intPath )



    def drawLargestCloud(self):

        """
        Draw a figure to discuss the largest moleuclar cloud
        :return:
        """
        intFITS="/home/qzyan/WORK/myDownloads/disQ2/intPath/Cloud227384_int.fits"
        maskFITS="/home/qzyan/WORK/myDownloads/disQ2/intPath/Cloud227384_mask.fits"

        fig = plt.figure(1, figsize=(35,8))
        rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica'], "size": 25 })
        # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

        rc('text', usetex=True)
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        dataCO,headCO= doFITS.readFITS( intFITS )
        dataMask,headMask= doFITS.readFITS( maskFITS )



        WCSCrop=WCS(headCO)
        
        grid_helper = pywcsgrid2.GridHelper(wcs=WCSCrop)

        # AxesGrid to display tow images side-by-side
 

        grid = ImageGrid(fig, (1, 1, 1), nrows_ncols=(1, 1),
                         cbar_mode="single", cbar_pad="0.3%",
                         cbar_location="right", cbar_size="1%",
                         axes_class=(pywcsgrid2.Axes, dict(header=WCSCrop)))
        
        
        axCO =         grid[0]  #pywcsgrid2.subplot(111, header=WCSCrop)
        cb_axes = grid.cbar_axes[0]  # colorbar axes


        cmapCO = plt.cm.bone
        cmapCO.set_bad('black')
        rms=1.8
        #im=axCO.imshow(  np.log(dataCO) , origin='lower', cmap=cmapCO,  norm=LogNorm(vmin=rms, vmax=rms*35),  interpolation='none')
        im=axCO.imshow(  np.log(dataCO) , origin='lower', cmap=cmapCO,   vmin = np.log(2), vmax = np.log(100),   interpolation='none')



        axCO.contour( dataMask ,[1], origin='lower',colors="red",linewidths= 0.3,alpha=0.5 )

        cb = cb_axes.colorbar(im)
        # cb_axes.axis["right"].toggle(ticklabels=True)
        cb_axes.set_ylabel(r"K km s$^{-1}$")
        cb_axes.set_xlabel("")
        tickesArray = np.asarray([ 2.5,5, 10, 20, 40, 80  ])
        # tickesArray=tickesArray**2
        cb.ax.set_yticks(np.log(tickesArray))
        #cb.ax.set_yticklabels(map(str, tickesArray))
        cb.ax.set_yticklabels( ["2.5","5", "10", "20", "40", "80"    ])

        #display Regions
        lRange1,bRange1=  box(137.8638900,3.9525955,6975.001,6000.001,0)
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$930^{+21}_{-19}$ pc"  )

        lRange1,bRange1=  box(127.0109528,-2.0675728,8028.376,6397.580,0)
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$926^{+17}_{-16}$ pc"  )

        lRange1,bRange1=  box(107.8594495,1.1121367,3104.504,3777.649,0)
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$796^{+37}_{-36}$ pc"  )

        lRange1,bRange1=  box(119.5058629,1.3143347,9328.786,8213.167,0)
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$695^{+26}_{-22}$ pc"  )




        lRange1,bRange1=  box(115.8929188,1.2735295,9793.059,9067.038,0)
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$846^{+10}_{-9}$ pc"  )



        lRange1,bRange1=  box(141.4359593,2.9489094,7384.201,5980.601,0)
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$1084^{+25}_{-27}$ pc"  )



        lRange1,bRange1=   box(133.5175787,2.8789620,4002.031,3893.189,0)   # 327.0+/-4 pc
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$346^{+27}_{-31}$ pc"  )

        lRange1,bRange1=    box(134.6664648,3.4999186,4002.031,3893.189,0)     # 327.0+/-4 pc
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$987^{+66}_{-70}$ pc"  )

        #########

        lRange1,bRange1=    box(129.5236198,3.2894556,4402.778,4541.667,0)
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$341^{+10}_{-9}$ pc"  )

        lRange1,bRange1=     box(128.1183343,2.5776500,4062.500,3000.000,0) 
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$1000^{+94}_{-79}$ pc"  )




        lRange1,bRange1=    box(122.0284110,2.4384642,4701.968,4430.700,0)
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$695^{+65}_{-62}$ pc"  )



        lRange1,bRange1=     box(138.3326235,2.1103346,6859.436,6254.772,0) 
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$994^{+34}_{-35}$ pc"  )

        lRange1,bRange1=     box(112.1445049,0.9304930,6952.280,5525.973,0) 
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$813^{+13}_{-12}$ pc"  )


        ############# 12
        lRange1,bRange1=     box(129.3362061,0.2729787,6944.433,6138.733,0)
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$927^{+20}_{-20}$ pc"  )

        lRange1,bRange1=      box(109.2206959,0.7334707,4303.439,4018.776,0)
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$760^{+18}_{-17}$ pc"  )

        lRange1,bRange1=     box(143.9981197,1.8963103,4672.695,3898.212,0)
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$927^{+54}_{-52}$ pc"  )


        lRange1,bRange1=      box(124.9077290,-2.9872898,3740.027,5338.542,0)
        self.showDistance( axCO,lRange1,bRange1, distanceStr= r"$953^{+42}_{-42}$ pc"  )


        ######## cloud around

        lRange,bRange =      box(123.2772635,-1.8503299,5409.934,4559.838,0)
        self.showDistance( axCO,lRange, bRange, distanceStr= r"$1127^{+24}_{-23}$ pc"  )

        lRange,bRange =       box(116.1756397,-2.2181734,7392.940,6729.796,0)
        self.showDistance( axCO,lRange, bRange, distanceStr= r"$360^{+9}_{-9}$ pc"  )

        lRange,bRange =       box(105.7252458,4.3841751,5802.083,5550.001,0)
        self.showDistance( axCO,lRange, bRange, distanceStr= r"$850^{+12}_{-12}$ pc"  )


        lRange,bRange =      box(136.4438973,-1.7631206,6871.908,6061.892,0)
        self.showDistance( axCO,lRange, bRange, distanceStr= r"$905^{+23}_{-22}$ pc"  )


        lRange,bRange =       box(139.9548056,-3.4351206,10911.600,9565.316,0)
        self.showDistance( axCO,lRange, bRange, distanceStr= r"$636^{+24}_{-22}$ pc"  )

        lRange,bRange =      box(146.7817001,-4.1019701,8675.000,7700.000,0)
        self.showDistance( axCO,lRange, bRange, distanceStr= r"$472^{+20}_{-16}$ pc"  )

        lRange,bRange =        box(148.1588470,2.6170538,14687.500,14427.083,0)
        self.showDistance( axCO,lRange, bRange, distanceStr= r"$216^{+14}_{-14}$ pc"  )

        lRange,bRange =         box(143.9353241,4.0589674,9563.079,8015.046,0)
        self.showDistance( axCO,lRange, bRange, distanceStr= r"$860^{+36}_{-40}$ pc"  )

        lRange,bRange =         box(134.6871099,0.2076384,5726.755,4963.188,0)
        self.showDistance( axCO,lRange, bRange, distanceStr= r"$834^{+42}_{-40}$ pc"  )

        lRange, bRange =  box(149.4269647,-1.5103541,4398.892,6453.813,0)
        self.showDistance(axCO, lRange, bRange, distanceStr=r"$1124^{+45}_{-40}$ pc")

        lRange, bRange =  box(124.8339230,-0.8728112,4036.458,3732.639,0) 
        self.showDistance(axCO, lRange, bRange, distanceStr=r"$912^{+46}_{-40}$ pc")

        lRange, bRange =   box(112.5211329,-2.8129080,14998.071,5774.981,0) 
        self.showDistance(axCO, lRange, bRange, distanceStr=r"$350^{+34}_{-22}$ pc")
        #RX13
        lRange, bRange =    box(105.6950933,-1.3000795,4360.651,4197.853,0) 
        self.showDistance(axCO, lRange, bRange, distanceStr=r"$717^{+50}_{-52}$ pc")
        #RX14
        lRange, bRange =    box(107.1700387,-0.1954013,3410.259,3296.801,0)
        self.showDistance(axCO, lRange, bRange, distanceStr=r"$690^{+41}_{-51}$ pc")

        #RX15
        lRange, bRange =     box(121.2507312,-3.8218493,6540.557,2692.579,0) 
        self.showDistance(axCO, lRange, bRange, distanceStr=r"$1083^{+31}_{-29}$ pc")

        #RX16
        lRange, bRange =    box(120.7435592,-1.4979573,5208.333,5121.527,0) 
        self.showDistance(axCO, lRange, bRange, distanceStr=r"$1120^{+67}_{-76}$ pc")

        #RX17
        lRange, bRange =      box(110.7551618,3.8044724,4667.928,4477.947,0)
        self.showDistance(axCO, lRange, bRange, distanceStr=r"$819^{+35}_{-35}$ pc")

        #RX18
        lRange, bRange =       box(134.0244930,-4.0537287,3906.250,3342.014,0) 
        self.showDistance(axCO, lRange, bRange, distanceStr=r"$842^{+72}_{-59}$ pc")




        #if 0:
        axCO.set_ticklabel1_type("absdeg",   locs=[105,110, 115, 120,125,130,135,140,145,150 ])
        axCO.set_ticklabel2_type("absdeg",   locs=[-5,-4,-3,-2,-1,0,1,2,3,4,5 ])

        #axCO.set_ticklabel_type("absdeg", "absdeg")
        axCO.axis[:].major_ticks.set_color("w")
        axCO.set_xlabel(r"Galactic Longitude ($^{\circ}$)")
        axCO.set_ylabel(r"Galactic Latitude ($^{\circ}$)")
        cloudName = "G125.1+02.6"

        ###########


        at = AnchoredText("G125.1+02.6", loc=4, frameon=False, prop={"color": "white"}  )
        axCO.add_artist(at)


        axCO["gal"].set_xticks( [120,130] )

        text = axCO["gal"].text( 142.5083,  1.6, "Cam OB1", color="white", alpha=0.8, fontsize=16,
                                     horizontalalignment='center', verticalalignment='center')
        text.set_path_effects([path_effects.Stroke(linewidth=2, foreground="black"), path_effects.Normal()])
        #add color bar



        print "Saving pdf..."

        #saveName=self.intFigurePath+"{}_{}CloudInt".format(calCode,ID)
        plt.savefig("largestCloud.pdf" , bbox_inches="tight"  )
        print "Saving png..."
        plt.savefig("largestCloud.png" , bbox_inches="tight" , dpi=300)



    def showDistance(self, axCO,lRange,bRange,distanceStr=""):
        """
        display the region box with green color
        :param axCO:
        :param lRange:
        :param bRange:
        :param distanceStr:
        :return:
        """
        lw = 2

        axCO["gal"].plot([lRange[0], lRange[0]], bRange, 'g--' ,   lw=1.0)
        axCO["gal"].plot([lRange[1], lRange[1]], bRange, 'g--', lw= 1.0 )

        axCO["gal"].plot(lRange, [bRange[0], bRange[0]], 'g--', lw= 1.0 )
        axCO["gal"].plot(lRange, [bRange[1], bRange[1]], 'g--', lw= 1.0 )
        text=axCO["gal"].text(np.mean(lRange) ,np.mean(bRange), distanceStr ,va="center" ,ha="center",color='green',fontsize=15)
        text.set_path_effects([path_effects.Stroke(linewidth= lw , foreground="black"), path_effects.Normal()])

    def testDistance(self, useAV=True):
        """
        A function that deal with distance
        :return:
        """
        print "Tesing distance"

        doDisDBSCAN = dbscanDis("disQ2", useAV)

        print doDisDBSCAN.maskPath
        print doDisDBSCAN.CO12FITS

        if 0: #distance tested
            pass
            #################### GOOD
            #lRange, bRange = doFITS.box(144.7837309,4.2655304,18729.167 ,5958.333 ,0)
            #useAV=False
            #doDisDBSCAN.calDisByID(360300, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=800, useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)


            #################### GOOD
            #lRange, bRange =  box(139.8176495,-3.2126518,13083.333,10000.000,0)
            #useAV=False
            #doDisDBSCAN.calDisByID(259585, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1000, useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)


            #################### The largest molecular cloud, need further inspect
            #lRange, bRange =  box(139.0814700,3.4082078,13475.000,8675.000,0) #region1  distance detected 958 pc
            #lRange, bRange =  box(110.4221234,1.3890335,9437.500,7562.500,0)   #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID(227384, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1500, useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)

            #################### Good
            #lRange, bRange =    box(136.3514000,-1.8665333,8418.240,6432.480,0)  #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID(260426, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1500, useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)

            #################### AG Failed
            #lRange, bRange = box(110.1143002,-0.8840061,6059.896,5341.339,0)  #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID(58773, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=3000, useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)


            #################### Good
            #lRange, bRange =  box(148.8437794,3.2313697,9994.338,9470.260,0)  #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID(539067, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=500, useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)


            #################### GooD
            #lRange, bRange =  box(112.3148208,-2.5667783,17020.833,12875.000,0)  #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID(412961 , useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1000, useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)

            #################### GOOD
            #lRange, bRange =  box(106.5445307,4.1271120,12378.472,7760.417,0)   #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID(264174 , useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1500, useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)


            #################### GOOD
            #lRange, bRange = box(115.5554637,-3.0557037,13644.000,11628.000,0)  #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID(402638  , useAV=useAV, lRange=lRange, bRange=bRange, cutDis=800, useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)



            #################### GOOD
            #lRange, bRange = box(122.4631200,-0.3014400,12312.000,13478.400,0)  #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID(234734  , useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1500, useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)

            #################### AG Faild
            #lRange, bRange = box(111.4064660,-3.1147597,8229.167,5017.361,0)  #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID(99653  , useAV=useAV, lRange=lRange, bRange=bRange, cutDis=3000 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)

            #################### AG Failed, AV detected
            #lRange, bRange =    box(118.1244078,2.9951216,7595.486,4065.394,0)    #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID( 238938  , useAV=useAV, lRange=lRange, bRange=bRange, cutDis=2000 ,   useHighAv=True, useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)




            #################### AV  dected with high av values  AG Faild
            #lRange, bRange =   box(143.2343728,0.8772509,9979.167,8520.833,0)    #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID( 454540  , useAV=useAV, lRange=lRange, bRange=bRange, cutDis=2000 , useBaseline=True, SL= 5 ,  useHighAv=True,   lExpand=0.5, bExpand=0.5, NL=3)




            #################### GOOD
            #lRange, bRange =   box(144.2706059,4.1616455,9461.806,6822.917,0)   #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID( 281664  , useAV=useAV, lRange=lRange, bRange=bRange, cutDis=2000 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)



             #################### AG Faild, truncated, Dected with HighAVvalues # AG Faild, truncated
            #lRange, bRange =   box(145.4477058,-0.2700691,8593.750,8038.194,0)  #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID( 264711  , useAV=useAV, lRange=lRange, bRange=bRange, cutDis=2000 , useHighAv=True,  useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)



            #################### AG GOOD
            #lRange, bRange =    box(105.7726078,0.7864653,7173.515,7685.909,0)   #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID( 403058  , useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1000 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)


            #################### AG Faild.AV detected
            #lRange, bRange =     box(114.7872894,3.6735803,2603.974,2361.031,0)    #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID( 222339  , useAV=useAV, lRange=lRange, bRange=bRange, cutDis= 2500 ,  useHighAv=True,  useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)


            #################### AG Truncated. AV faild
            #lRange, bRange =     box(133.1233857,0.4086375,3370.949,3170.868,0)   #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID( 89568   , useAV=useAV, lRange=lRange, bRange=bRange, cutDis= 3000 , useHighAv=True,  useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)




            #################### # AG  GooD
            lRange, bRange =    box(148.1421364,0.0953276,7263.937,6761.590,0)      #region2  749 pc, do they belong the same cloud?
            useAV= False
            doDisDBSCAN.calDisByID( 295121   , useAV=useAV, lRange=lRange, bRange=bRange, cutDis= 1500 ,  useHighAv=False, useBaseline=True, SL=  7 ,    lExpand=0.5, bExpand=0.5, NL=2 )







            #################### AV GOOD
            #lRange, bRange =    box(147.8746820,-4.4293344,7414.641,5401.235,0)    #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID( 215739   , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True,  cutDis= 1500 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)


            #################### AV AG failed
            #lRange, bRange =     box(148.5820281,-0.1852218,9270.823,3758.435,0)      #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID( 181901   , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True,  cutDis= 2500 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)

            #################### AV dected
            #lRange, bRange =  box(109.7146764,2.5990621,4687.500,3833.912,0)    #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID( 458317   , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True,  cutDis= 1500 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)


            #################### Failed
            #lRange, bRange =   box(144.4192462,-0.8854382,2612.204,2220.374,0)   #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID(  181719    , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True,  cutDis= 3500 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)

            #################### AG GOOD
            #lRange, bRange =   box(148.9838107,2.5939286,7306.134,6727.431,0)    #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID(  445141    , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False,  cutDis= 1000 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)


            #################### AV and AG faild
            #lRange, bRange =   box(123.3870426,-0.8092852,4687.500,4007.523,0)    #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID(  113910    , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False,  cutDis= 3500 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)

            ####################Both AV and AG failed
            #lRange, bRange = box(142.2383456,-0.4439023,3448.110,3243.152,0)   #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID(  260428    , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True,  cutDis= 3000 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)


            #################### AV and AG faild
            #lRange, bRange =   box(143.5557472,-2.4222247,3399.884,2922.454,0)   #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID(  258072    , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True,  cutDis= 3000 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)

            #################### AG GOOD
            #lRange, bRange =  box(141.3882957,-3.5589113,5791.667,6208.333,0)   #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID(  252392    , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False,  cutDis= 2000 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)




            #################### AG GOOD
            #lRange, bRange =     box(146.2001185,-4.3333877,7052.951,6510.417,0)    #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID(  325105    , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True,  cutDis= 1000 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)



            #################### AG GOOD
            #lRange, bRange =      box(144.2079849,-3.8978160,7042.904,5505.723,0)     #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID(  339429    , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True,  cutDis= 1000 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)

            #################### Failed W5?
            #lRange, bRange =       box(138.4173378,1.5443265,2686.161,2011.714,0)  #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID(  158599    , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False,  cutDis= 2500 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=3)

            #################### Jump unclear
            #lRange, bRange =         box(110.3977986,3.6680938,8391.204,3436.053,0)   #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID(  521860    , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True,  cutDis= 1500 , useBaseline=True, SL=  5   ,    lExpand=0.5, bExpand=0.5, NL=3)


            #################### > 2 kpc few background
            #lRange, bRange =          box(136.7522895,1.2777690,4462.516,3382.470,0)    #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID(  139159    , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True,  cutDis= 3500 , useBaseline=True, SL=  5   ,    lExpand=0.5, bExpand=0.5, NL=3)

            #################### AG GOOD
            #lRange, bRange =         box(131.6478138,-1.5676965,4042.219,3305.443,0)    #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID(  267103    , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True,  cutDis= 1500 , useBaseline=True, SL=  5   ,    lExpand=0.5, bExpand=0.5, NL=3)


            #################### Too close?
            #lRange, bRange =         box(147.7904842,-1.1587058,6452.546,4325.810,0)     #region2  749 pc, do they belong the same cloud?
            #useAV=False
            #doDisDBSCAN.calDisByID(   460233     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True,  cutDis= 1000 , useBaseline=True, SL=  4   ,    lExpand=0.5, bExpand=0.5, NL=3)



            #################### AV GOOD
            #lRange, bRange =        box(147.3862134,-3.1331406,6215.278,6319.444,0)     #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(   495686     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False,  cutDis= 800 , useBaseline=True, SL=  4   ,    lExpand=0.5, bExpand=0.5, NL= 2  )

            #################### AG GOOD
            #lRange, bRange =        box(143.3018896,2.9311760,5895.544,4605.517,0)     #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(   285283     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False,  cutDis= 2000 , useBaseline=True, SL=  5   ,    lExpand=0.5, bExpand=0.5, NL= 2  )

            #################### AG GOOD
            #lRange, bRange =        box(130.3226687,0.6311805,5058.836,3284.144,0)    #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(   293789     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False,  cutDis= 2000 , useBaseline=True, SL=  5   ,    lExpand=0.5, bExpand=0.5, NL= 2  )


            ####################  AV GOOD
            #lRange, bRange =       box(112.8452579,-0.7326747,3522.941,2984.313,0)   #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(   109773     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False,  cutDis= 3000 , useBaseline=True, SL=  5   ,    lExpand=0.5, bExpand=0.5, NL= 2  )


            #################### FFF
            #lRange, bRange =   box(117.2308905,2.3645485,2662.439,3868.072,0)           #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(   354740      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False,  cutDis= 2000 , useBaseline=True, SL=  6   ,    lExpand=0.5, bExpand=0.5, NL= 2  )


            #################### FFF, too weak
            #lRange, bRange =   box(134.9593570,-0.4954535,4211.342,3918.306,0)          #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(   127518      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False,  cutDis= 2500 , useBaseline=True, SL= 5   ,    lExpand=0.5, bExpand=0.5, NL= 2  )

            #################### FFFF
            #lRange, bRange =    box(130.4141411,4.6115989,5234.455,4571.357,0)         #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(   267827      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False ,  cutDis= 2500 , useBaseline=True, SL= 5   ,    lExpand=0.5, bExpand=0.5, NL= 2  )



            ####################  AV GOOD
            #lRange, bRange =     box(133.2848100,-2.6336915,4448.785,3592.785,0)          #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(   277600      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False ,  cutDis= 2000 , useBaseline=True, SL= 5   ,    lExpand=0.5, bExpand=0.5, NL= 2  )




            #################### FFF, far
            #lRange, bRange =      box(120.9329469,1.2343663,3982.728,2944.107,0)           #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(   56974      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False ,  cutDis= 3000 , useBaseline=True, SL= 5   ,    lExpand=0.5, bExpand=0.5, NL= 2  )



            #################### FFF, far
            #lRange, bRange =       box(116.9512307,-3.3349109,3908.259,3536.523,0)          #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(   164403      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False ,  cutDis= 3000 , useBaseline=True, SL= 5   ,    lExpand=0.5, bExpand=0.5, NL= 2  )


            #################### FFF, far
            #lRange, bRange =        box(107.7562989,1.9320080,4210.069,3298.611,0)         #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(   278708      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,  cutDis= 2000 , useBaseline=True, SL= 6 ,    lExpand=0.5, bExpand=0.5, NL= 2  )

            #################### FFF, far
            #lRange, bRange =        box(144.5973527,-2.5004614,4579.730,3725.740,0)         #region2  749 pc, do they belong the same cloud?
            #useAV=True
            #doDisDBSCAN.calDisByID(   247603      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,  cutDis= 2000 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL= 2  )

            #################### AG GOOD
            #lRange, bRange =         box(129.9984357,-2.2159297,7660.791,3642.015,0)        #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(   263681      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,  cutDis= 2000 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL= 2  )


            #################### FFF, weak or far, no detection
            #lRange, bRange =         box(105.4267963,0.2900904,4719.193,3425.533,0)        #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(   72973      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,  cutDis= 3000 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL= 2  )


            #################### FFF, locates a cluster, give up
            #lRange, bRange =         box(133.4649456,-1.6283099,4581.404,5172.164,0)        #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(   454189      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 1000 , useBaseline=True, SL= 3 ,    lExpand=0.5, bExpand=0.5, NL=1   )


            #################### AG GOOD
            #lRange, bRange =        box(115.9509875,-0.5314069,5309.606,3964.120,0)      #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(   371984      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 800 , useBaseline=True, SL= 3 ,    lExpand=0.5, bExpand=0.5, NL=2   )

            #################### Too few foreground stars FFF
            #lRange, bRange =        box(135.8951785,4.7437531,6677.028,3586.199,0)      #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(   446495      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 800 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL=2   )



            #################### two componennts.... FFF
            #lRange, bRange =      box(110.9392031,3.8172988,4960.676,2372.194,0)      #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(   363305      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 2000 , useBaseline=True, SL= 6 ,    lExpand=0.5, bExpand=0.5, NL=2   )


            #################### Jump unclear, FFF
            #lRange, bRange =       box(114.7798879,0.9069096,5129.051,5129.051,0)     #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(   457569      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 2500 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL= 3   )

            #################### FFF, unclear jump
            #lRange, bRange =       box(132.9653125,1.5062032,5132.311,2930.357,0)      #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(   422005      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 2500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### AG GOOD
            #lRange, bRange =        box(146.2876565,-2.3809948,6293.396,3098.473,0)      #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(   275042      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 2000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### AG GOOD
            #lRange, bRange =         box(106.8876165,0.9620629,6510.417,2742.814,0)       #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(   277035      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 1000 , useBaseline=True, SL=  4  ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            #################### AG GOOD
            #lRange, bRange =         box(129.1436150,-0.1023395,5497.685,3954.475,0)     #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(   163931      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 2000 , useBaseline=True, SL=  5  ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            #################### AG GOOD
            #lRange, bRange =         box(130.9951651,-0.9276344,5987.976,3154.739,0)      #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(   312594      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 1500 , useBaseline=True, SL=  5  ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### Far, FFFFF
            #lRange, bRange =         box(108.4459180,-1.1094894,5684.893,2988.964,0)     #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(   85790      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 2500 , useBaseline=True, SL=  5  ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### AG GOOD
            #lRange, bRange =         box(149.5081494,-1.1220855,3566.663,4738.806,0)      #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(   296464      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 2000 , useBaseline=True, SL=  5  ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### AG GOOD
            #lRange, bRange =        box(120.9062689,-1.4743259,4895.557,4005.984,0)       #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(    237970      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 2000 , useBaseline=True, SL=  5  ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### Far, FFF
            #lRange, bRange =        box(115.6790574,-1.5794234,6317.515,2965.856,0)       #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(     131040      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 3000 , useBaseline=True, SL=  5  ,    lExpand=0.5, bExpand=0.5, NL= 2   )




            #################### AG GOOD
            #lRange, bRange =        box(121.9820486,-0.4433275,5280.671,5613.426,0)        #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(     516881      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 1000 , useBaseline=True, SL=  5  ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### Weak FFF
            #lRange, bRange =        box(144.1861969,0.1836835,5895.833,3625.000,0)       #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(     612883      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,  cutDis= 3000 , useBaseline=True, SL=  4  ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            #################### Far, FFF
            #lRange, bRange =         box(140.6235781,0.3970823,3198.276,6086.771,0)       #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(     157843      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,  cutDis= 3000 , useBaseline=True, SL=  5  ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### Far, FFF
            #lRange, bRange =         box(106.0727873,1.5860770,5333.333,3833.333,0)      #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(     225400      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,  cutDis= 3000 , useBaseline=True, SL=  5  ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### AG GOOD, do not use mask
            #lRange, bRange =          box(142.1204632,-3.5906363,10517.940,5685.764,0)        #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(     248354      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,  useMask=False, cutDis= 1500 , useBaseline=True, SL=  5  ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### FFF, not enough on-cloud stars
            #lRange, bRange =           box(131.8877566,0.0700062,3872.042,4989.520,0)        #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(     267299      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  useMask=True, cutDis= 2000 , useBaseline=True, SL=  5  ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### Far, FFF
            #lRange, bRange =          box(111.4239398,2.3489552,4906.925,4991.319,0)       #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(     84386      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  useMask=True, cutDis= 3000 , useBaseline=True, SL=  5  ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### FFF, weak
            #lRange, bRange =         box(144.5804192,-0.2669718,4571.759,5107.060,0)       #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(     478786      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,  cutDis= 2500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            #################### FFF, unclear jump
            #lRange, bRange =          box(143.8093442,-3.5602895,5750.000,5520.833,0)         #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(     186740      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  , useMask=False, cutDis= 2000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 1   )




            #################### Far, FFF
            #lRange, bRange =          box(111.6652347,4.0336319,4311.811,4521.123,0)        #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(     209687      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  , useMask=False, cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            ####################FFF, weak
            #lRange, bRange =          box(135.7287973,-1.2425991,4572.997,3434.630,0)        #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(     496036      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 2500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### insufficient foreground stars, FFF
            #lRange, bRange =          box(113.1678539,-3.6948210,3882.137,5196.277,0)         #region2  749 pc, do they belong the same cloud?
            #useAV= False
            #doDisDBSCAN.calDisByID(     424151      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 1000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 1   )


            ####################Far, FFF
            #lRange, bRange =           box(120.7053606,-0.3943630,4633.918,3558.291,0)         #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(     105966      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            ####################Far, FFF
            #lRange, bRange =           box(109.6259750,1.8562648,4050.560,4079.631,0)         #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(     66721      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,  cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            #################### AV GOOD
            #lRange, bRange =           box(109.6259750,1.8562648,4050.560,4079.631,0)         #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(     66721      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,  cutDis= 1500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### FFF, two few on-cloud stars
            #lRange, bRange =            box(146.0636858,-4.3222124,5564.656,6660.470,0)       #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(     446969      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,  useMask=False, cutDis= 1500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### unclear jump, FFF
            #lRange, bRange =         box(117.9247002,4.9294781,4077.424,2291.360,0)       #region2  749 pc, do they belong the same cloud?
            #useAV= True
            #doDisDBSCAN.calDisByID(      512057      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### FFF, unclear jump
            #lRange, bRange =         box(139.5482955,3.1869228,5427.357,5425.347,0)        #
            #useAV= True
            #doDisDBSCAN.calDisByID(      94917      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 1   )



            #################### FFF, nodetection
            #lRange, bRange =         box(123.5918108,1.9833821,3672.357,4603.106,0)        #
            #useAV= False
            #doDisDBSCAN.calDisByID(      413699      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,   cutDis= 1500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            ####################FFF, nodetection, jump unclear
            lRange, bRange =         box(138.3404591,-2.0669541,3827.884,4028.823,0)        #
            useAV= True
            doDisDBSCAN.calDisByID(      112625      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### FFF, unclear jump
            lRange, bRange =         box(141.2820727,3.7949422,3968.541,4189.574,0)         #
            useAV= True
            doDisDBSCAN.calDisByID(      294369      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,   cutDis= 2000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            #################### FFF, insufficient on-cloud stars
            lRange, bRange =         box(114.2193552,-2.6361066,6017.232,5187.676,0)       #
            useAV= False
            doDisDBSCAN.calDisByID(      540708      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,   cutDis= 1000 , useBaseline=True, SL=  3,    lExpand=0.5, bExpand=0.5, NL= 1   )


            #################### FFF, unclear jump
            lRange, bRange =          box(145.3603100,-3.3048719,4392.361,3559.028,0)    #
            useAV= True
            doDisDBSCAN.calDisByID(      252482      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,   cutDis= 1500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            #################### FFF, no detection unclear jump
            lRange, bRange =         box(123.5178437,3.2158928,4842.547,2598.824,0)        #
            useAV= True
            doDisDBSCAN.calDisByID(      57909      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### FFF,unclear jump
            lRange, bRange =         box(117.4739209,-0.6087911,4137.385,3007.105,0)        #
            useAV= True
            doDisDBSCAN.calDisByID(      505158      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### FFF, no detection
            lRange, bRange =         box(138.5460378,1.4939173,4742.155,3858.025,0)        #
            useAV= True
            doDisDBSCAN.calDisByID(      422095      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### FFF, no detection
            lRange, bRange =         box(119.3940225,1.4164086,3982.728,4197.853,0)        #
            useAV= True
            doDisDBSCAN.calDisByID(      65975      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            ####################  AG GOOD
            lRange, bRange =         box(125.6333998,-2.5121796,5304.784,3858.025,0)        #
            useAV= False
            doDisDBSCAN.calDisByID(      236543      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,   cutDis= 2500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### FFF, far
            lRange, bRange =         box(131.9834973,-0.8129730,3884.817,3524.801,0)        #
            useAV= True
            doDisDBSCAN.calDisByID(      85005      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### FFF, no detection
            lRange, bRange =         box(149.9694350,3.6617247,2020.488,4056.857,0)        #
            useAV= True
            doDisDBSCAN.calDisByID(      330607      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### FFF, no detection
            lRange, bRange =         box(132.1596658,-0.3628601,3901.329,3021.989,0)        #
            useAV= True
            doDisDBSCAN.calDisByID(      560529      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### FFF, AG truncated
            lRange, bRange =         box(138.7753004,2.4189587,5104.167,3385.417,0)        #
            useAV= True
            doDisDBSCAN.calDisByID(      527583      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv= True   ,   cutDis= 2500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            #################### FFF, unclear jump
            lRange, bRange =         box(147.0795112,-4.2064376,4029.241,3860.629,0)        #
            useAV= True
            doDisDBSCAN.calDisByID(      355805      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,   cutDis= 700 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### FFF, far
            lRange, bRange =   box(117.1318929,-3.3423496,10323.307,4439.262,0)            #
            useAV= False
            doDisDBSCAN.calDisByID(      146753      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  ,  useMask=False,  cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### FFF, unclear jump
            lRange, bRange =         box(139.5257556,-1.7192524,3906.250,5659.722,0)         #
            useAV= True
            doDisDBSCAN.calDisByID(      116146      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv= True    ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            #################### AV GOOD
            lRange, bRange =         box(135.7161406,3.4395674,3906.250,3420.139,0)        #
            useAV= True
            doDisDBSCAN.calDisByID(      316030      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1000 , useBaseline=True, SL=  4  ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            #################### unclear jump, FFF
            lRange, bRange = box(142.6796662, -2.3218424, 4123.264, 3460.166, 0)
            useAV = False

            doDisDBSCAN.calDisByID(260265, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=2000,   useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)


            #################### AG GOOD
            lRange, bRange =         box(106.4580168,1.6852660,5353.009,2387.153,0)
            useAV= False
            doDisDBSCAN.calDisByID(      264979      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################FFF, nodetection
            lRange, bRange =          box(119.9315472,3.4908529,2871.750,3164.786,0)
            useAV= True
            doDisDBSCAN.calDisByID(      45301      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  True  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### FFF, unclear jump
            lRange, bRange =         box(149.6114494,3.5437226,4579.730,4345.301,0)
            useAV= True
            doDisDBSCAN.calDisByID(     324398       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  True  , useMask=False,  cutDis= 2000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################  Far
            lRange, bRange =          box(107.2527436,0.0981500,3089.434,3114.551,0)
            useAV= False
            doDisDBSCAN.calDisByID(      89499      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### unclear jump FFF
            lRange, bRange =         box(142.7489206,-0.9633617,3407.587,3173.158,0)
            useAV= True
            doDisDBSCAN.calDisByID(      140543      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  True  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### unclear jump FFF
            lRange, bRange =         box(128.7403614,1.3197981,4108.796,4021.991,0)
            useAV= True
            doDisDBSCAN.calDisByID(      75782      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################unclear jump FFF
            lRange, bRange =         box(146.2158827,-0.2688281,3508.056,3390.842,0)
            useAV= True
            doDisDBSCAN.calDisByID(      206164      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  True  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### AV GOOD
            lRange, bRange =         box(129.3505866,-1.6089215,3677.180,3345.631,0)
            useAV= True
            doDisDBSCAN.calDisByID(      303216      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL=  5  ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### AV GOOD
            lRange, bRange =         box(134.0382861,-4.0111832,3701.292,3339.603,0)
            useAV= True
            doDisDBSCAN.calDisByID(      270697      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### Far, FFF
            lRange, bRange =         box(111.1232753,3.4355762,3295.396,3345.631,0)
            useAV= True
            doDisDBSCAN.calDisByID(      39360      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### FFF, no detection
            lRange, bRange =         box(119.3350251,-1.6256208,4144.362,3198.276,0)
            useAV= True
            doDisDBSCAN.calDisByID(      186989      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### FFF, no detection
            lRange, bRange =    box(139.7397373,-0.5792985,4269.949,2721.046,0)
            useAV= True
            doDisDBSCAN.calDisByID(      152000      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### #unclear detection
            lRange, bRange =         box(107.0681115,2.4489860,4169.480,3024.129,0)
            useAV= True
            doDisDBSCAN.calDisByID(      277918      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### unclear jump
            lRange, bRange =         box(132.8939796,-1.6373980,3566.663,2796.398,0)
            useAV= True
            doDisDBSCAN.calDisByID(      307250      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### unclear jump
            lRange, bRange =         box(109.5524391,3.7051211,3488.521,3197.811,0)
            useAV= True
            doDisDBSCAN.calDisByID(      294302      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  True  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### unclear jump
            lRange, bRange =         box(144.1019012,-0.1508360,2972.220,4027.148,0)
            useAV= True
            doDisDBSCAN.calDisByID(      298021      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  True  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################  unclear jump
            lRange, bRange =         box(141.7249063,-1.0849492,3466.194,3315.490,0)
            useAV= True
            doDisDBSCAN.calDisByID(      148464      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### unclear jump
            lRange, bRange =         box(117.4854316,-0.1596701,4036.458,2965.856,0)
            useAV= True
            doDisDBSCAN.calDisByID(         121807     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### unclear jump
            lRange, bRange =         box(148.5310623,-1.4520640,4145.833,3458.333,0)
            useAV= True
            doDisDBSCAN.calDisByID(      227203      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### AG GOOD
            lRange, bRange =         box(108.3815358,0.2959505,3978.588,3327.546,0)
            useAV= False
            doDisDBSCAN.calDisByID(      307747      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### AG GOOD
            lRange, bRange =         box(105.8580949,-1.1969070,4310.137,4430.700,0)
            useAV= False
            doDisDBSCAN.calDisByID(      331873      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            ####################  unclear jump
            lRange, bRange =         box(105.1540601,2.7700257,2929.311,3213.007,0)
            useAV= True
            doDisDBSCAN.calDisByID(      236038      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### unclear jump
            lRange, bRange =         box(122.3563095,3.2550539,3304.404,2907.100,0)
            useAV= True
            doDisDBSCAN.calDisByID(      73236      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### unclear jump
            lRange, bRange =         box(122.5508685,1.5399114,4633.849,2586.082,0)
            useAV= True
            doDisDBSCAN.calDisByID(      76666      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### unclear jump
            lRange, bRange =         box(111.9238379,3.0286232,3046.641,3651.318,0)
            useAV= True
            doDisDBSCAN.calDisByID(      20525      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### unclear jump
            lRange, bRange =         box(114.6561785,-0.3492048,3275.333,2684.223,0)
            useAV= True
            doDisDBSCAN.calDisByID(      106840      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            #################### unclear jump
            lRange, bRange =         box(117.0693208,-2.2429357,3817.837,2813.143,0)
            useAV= True
            doDisDBSCAN.calDisByID(      129415      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### unclear jump
            lRange, bRange =         box(147.2969292,0.9177794,3592.785,3074.363,0)
            useAV= True
            doDisDBSCAN.calDisByID(     196396       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            #################### twoo few stars
            lRange, bRange =         box(146.9436881,-0.1132803,4079.861,3333.333,0)
            useAV= True
            doDisDBSCAN.calDisByID(     616826       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4  ,    lExpand=0.5, bExpand=0.5, NL= 2   )



            #################### Unclear jump
            lRange, bRange =         box(112.2122490,4.8955276,3613.704,2499.111,0)
            useAV= True
            doDisDBSCAN.calDisByID(     302890       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### Unclear jump
            lRange, bRange =         box(140.7427134,-1.0972318,3004.035,3446.100,0)
            useAV= True
            doDisDBSCAN.calDisByID(     146162       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  True  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### Unclear jump
            lRange, bRange =         box(144.4288605,0.4825094,3606.851,3305.443,0)
            useAV= True
            doDisDBSCAN.calDisByID(     150957       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### Unclear jump
            lRange, bRange =         box(130.6294775,1.8736745,3423.997,3436.053,0)
            useAV= True
            doDisDBSCAN.calDisByID(     130673       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### Unclear jump
            lRange, bRange =         box(114.0450049,3.4210679,3427.956,2955.552,0)
            useAV= True
            doDisDBSCAN.calDisByID(     330649       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500 , useBaseline=True, SL=  5  ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### AG GOOD
            lRange, bRange =         box(116.9641660,3.6574644,3935.051,3148.041,0)
            useAV= False
            doDisDBSCAN.calDisByID(     309918       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### Unclear jump
            lRange, bRange =         box(114.9638648,-1.0445918,3014.082,2888.495,0)
            useAV= True
            doDisDBSCAN.calDisByID(     421529       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### unclear jump
            lRange, bRange =         box(143.1124119,0.9009793,4751.365,2183.814,0)
            useAV=True
            doDisDBSCAN.calDisByID(     143715       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  True  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### unclear jump
            lRange, bRange =         box(116.4341827,0.4355616,3356.481,2792.245,0)
            useAV= True
            doDisDBSCAN.calDisByID(     405772       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            ####################  unclear jump
            lRange, bRange =         box(106.4280759,1.0347469,3888.166,2873.425,0)
            useAV= True
            doDisDBSCAN.calDisByID(     64256       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            ####################  unclear jump
            lRange, bRange =         box(138.7231660,-0.2986638,2437.500,4312.500,0)
            useAV= True
            doDisDBSCAN.calDisByID(     312792       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            ####################  unclear jump
            lRange, bRange =         box(110.9116959,3.1523015,3657.086,2793.049,0)
            useAV= False
            doDisDBSCAN.calDisByID(     210446       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### unclear jump
            lRange, bRange =         box(137.7904910,1.5678757,2662.439,-2592.110,0)
            useAV= True
            doDisDBSCAN.calDisByID(     147240       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            ####################  AG GOOD
            lRange, bRange =         box(114.6841598,-0.1312470,4064.368,3277.906,0)
            useAV= False
            doDisDBSCAN.calDisByID(      178245        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=   5  ,    lExpand=0.5, bExpand=0.5, NL=  1   )

            #################### unclear jump
            lRange, bRange =         box(127.9205510,-0.0378375,4311.343,2575.231,0)
            useAV= True
            doDisDBSCAN.calDisByID(     112690       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            ####################  unclear jump
            lRange, bRange =         box(111.7175330,2.5070752,3556.617,2913.612,0)
            useAV= False
            doDisDBSCAN.calDisByID(     476834       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=   4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### unclear jump
            lRange, bRange =         box(122.3416083,2.7083862,3660.301,2473.958,0)
            useAV= False
            doDisDBSCAN.calDisByID(     71392       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### AG GOOD
            lRange, bRange =         box(111.0382748,3.5244203,2235.444,3951.796,0)
            useAV= False
            doDisDBSCAN.calDisByID(     237244       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### no detection
            lRange, bRange =         box(113.7535870,3.5328420,4390.512,3034.176,0)
            useAV= True
            doDisDBSCAN.calDisByID(     249018       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### no detection
            lRange, bRange =         box(119.5175499,0.5535538,4026.808,3942.414,0)
            useAV= True
            doDisDBSCAN.calDisByID(     96329       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### un clear jump
            lRange, bRange =         box(113.4323727,-0.9482429,3238.510,2755.931,0)
            useAV= False
            doDisDBSCAN.calDisByID(     439533       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            ####################  AG GOOD
            lRange, bRange = box(145.9873551,1.3676858,5086.806,3871.528,0)    #region2  749 pc, do they belong the same cloud?
            useAV= False
            doDisDBSCAN.calDisByID( 160462   , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True,  cutDis= 2500 , useBaseline=True, SL=5 ,    lExpand=0.5, bExpand=0.5, NL= 2 )


            #################### Largest cloud
            lRange, bRange =         box(141.4359593,2.9489094,7384.201,5980.601,0)
            useAV= False
            doDisDBSCAN.calDisByID(      227384        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False, cutDis= 1500 , useBaseline=True, SL=  6 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### FFF
            lRange, bRange =         box(118.1028276,1.6850730,3114.551,3265.255,0)
            useAV= False
            doDisDBSCAN.calDisByID(     338264       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            ####################FFF
            lRange, bRange =         box(108.6869432,4.5675669,3432.704,2662.439,0)
            useAV= True
            doDisDBSCAN.calDisByID(     350539       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =         box(129.0669141,1.6458749,3906.250,3003.472,0)
            useAV= True
            doDisDBSCAN.calDisByID(     326903       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL= 4  ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            lRange, bRange =         box(131.2046272,-0.7257047,3322.298,2707.999,0)
            useAV= True
            doDisDBSCAN.calDisByID(     400399       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =         box(149.1863516,-2.1093134,3457.755,3168.403,0)
            useAV= True
            doDisDBSCAN.calDisByID(     307209       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False, cutDis= 2000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =         box(106.3150762,-0.1798487,3559.028,3081.597,0)
            useAV= True
            doDisDBSCAN.calDisByID(     83175       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False,  cutDis= 2000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =         box(122.6698714,-3.0212587,5987.221,5639.254,0)
            useAV= True
            doDisDBSCAN.calDisByID(     308688       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  , useMask=False,  cutDis= 2000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### AG GOOD
            lRange, bRange =         box(107.7561084,2.9155225,2883.137,2699.166,0)
            useAV= False
            doDisDBSCAN.calDisByID(     282690       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  True  ,   cutDis= 1500 , useBaseline=True, SL= 5,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =         box(136.4215005,-1.4268612,8598.990,2073.415,0)
            useAV= False
            doDisDBSCAN.calDisByID(     532327       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  , useMask=False,  cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =         box(131.9612055,-0.5246280,2724.328,3073.560,0)
            useAV= True
            doDisDBSCAN.calDisByID(     408144       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =         box(131.9612055,-0.5246280,2724.328,3073.560,0)
            useAV= True
            doDisDBSCAN.calDisByID(     408144       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 800 , useBaseline=True, SL= 3 ,    lExpand=0.5, bExpand=0.5, NL= 1   )


            #################### intital
            lRange, bRange =         box(146.2543506,-1.6086187,3662.947,2965.243,0)
            useAV= False
            doDisDBSCAN.calDisByID(     618036       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =         box(120.5753388,1.9530628,2553.597,3407.587,0)
            useAV= False
            doDisDBSCAN.calDisByID(     63453       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =         box(144.9740841,1.5794589,3976.913,2712.674,0)
            useAV= True
            doDisDBSCAN.calDisByID(     127798       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### AG GOOD
            lRange, bRange =         box(124.6732187,-1.5100658,3144.514,2572.784,0)
            useAV= False
            doDisDBSCAN.calDisByID(     251759       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### weak
            lRange, bRange =         box(105.0560352,2.4683309,2184.606,2980.324,0)
            useAV= True
            doDisDBSCAN.calDisByID(     266345       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL= 3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =         box(105.0117678,3.2904111,1892.679,3152.096,0)
            useAV= True
            doDisDBSCAN.calDisByID(     266345       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =         box(123.3967055,4.9966831,3684.981,1826.338,0)
            useAV= False
            doDisDBSCAN.calDisByID(     531701       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =         box(105.0850592,2.0301436,2365.426,2342.224,0)
            useAV= True
            doDisDBSCAN.calDisByID(     232272       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            ####################
            lRange, bRange =         box(121.6882942,2.8305192,3488.521,2517.549,0)
            useAV= True
            doDisDBSCAN.calDisByID(     49347       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### AG GOOD
            lRange, bRange =         box(121.7725262,-0.4722424,3300.141,3021.059,0)
            useAV= False
            doDisDBSCAN.calDisByID(     271145       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =         box(114.2574531,-0.9922566,2853.610,2790.816,0)
            useAV= True
            doDisDBSCAN.calDisByID(     197476       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange = box(125.2246801, 0.5804885, 2534.695, 2790.009, 0)
            useAV = True
            doDisDBSCAN.calDisByID(48294, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=3000,
                                   useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)

            ####################
            lRange, bRange =         box(148.2339018,1.3020552,3499.991,3229.158,0)
            useAV= False
            doDisDBSCAN.calDisByID(     454748       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(134.3623138,-1.7122963,3009.259,2835.648,0)
            useAV= True
            doDisDBSCAN.calDisByID(     114569       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =          box(138.3848589,-3.6507347,2395.451,2069.856,0)
            useAV= True
            doDisDBSCAN.calDisByID(     334123       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(129.9556027,2.1367605,2797.068,2604.167,0)
            useAV= True
            doDisDBSCAN.calDisByID(     444650       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################  AV GOOD
            lRange, bRange =        box(111.0855673,-0.1725342,3017.618,3437.902,0)
            useAV= True
            doDisDBSCAN.calDisByID(     309300       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  , useMask=True,   cutDis= 2000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            ####################
            lRange, bRange =          box(144.3733601,-1.5176945,3457.822,1883.801,0)
            useAV= True
            doDisDBSCAN.calDisByID(     589199       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =          box(137.3472756,4.9475624,2644.299,2155.906,0)
            useAV= True
            doDisDBSCAN.calDisByID(     596922       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =         box(147.9461907,3.2538674,2531.829,2531.829,0)
            useAV= True
            doDisDBSCAN.calDisByID(     480124       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            ####################
            lRange, bRange =         box(105.1485142,1.1642008,2848.958,2311.145,0)
            useAV= True
            doDisDBSCAN.calDisByID(     104193       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            ####################
            lRange, bRange =           box(138.9708307,0.9009631,3226.273,2416.088,0)
            useAV= True
            doDisDBSCAN.calDisByID(     351851       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =         box(108.8281296,1.6276420,2924.311,2897.691,0)
            useAV= True
            doDisDBSCAN.calDisByID(     199054       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(146.8241128,-0.7724909,2846.633,2771.281,0)
            useAV= True
            doDisDBSCAN.calDisByID(     603763       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(132.2771623,0.5921677,2560.574,2707.092,0)
            useAV= True
            doDisDBSCAN.calDisByID(     311059       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(142.9670811,3.2712099,3074.363,2748.843,0)
            useAV= True
            doDisDBSCAN.calDisByID(     551433       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(120.6881100,0.1710584,2652.392,2290.702,0)
            useAV= True
            doDisDBSCAN.calDisByID(     353994       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### AG GOOD
            lRange, bRange =           box(122.8705620,-2.2955322,2704.301,2260.561,0)
            useAV= False
            doDisDBSCAN.calDisByID(     270975       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =          box(148.1460128,1.0988682,2692.573,2772.948,0)
            useAV= True
            doDisDBSCAN.calDisByID(     609269       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(144.7983037,-4.5395466,2104.741,3029.199,0)
            useAV= True
            doDisDBSCAN.calDisByID(     405109       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(140.7544581,2.0714168,2728.023,2490.804,0)
            useAV= True
            doDisDBSCAN.calDisByID(     460870       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=   3  ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(143.2920751,-3.1187192,2670.811,2520.107,0)
            useAV= True
            doDisDBSCAN.calDisByID(     249180       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL= 3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(137.8327112,-0.7362822,2387.150,2647.567,0)
            useAV= True
            doDisDBSCAN.calDisByID(     133093       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =        box(140.9054996,3.0182443,2724.730,2278.646,0)
            useAV= True
            doDisDBSCAN.calDisByID(     576758       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(126.5856594,-1.4054298,2811.748,1834.962,0)
            useAV= False
            doDisDBSCAN.calDisByID(     150563       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2500  , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(108.0665450,3.3525067,2365.217,2030.319,0) 
            useAV= True
            doDisDBSCAN.calDisByID(     357534       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(116.3551418,3.1279619,2558.248,2552.434,0)
            useAV= True
            doDisDBSCAN.calDisByID(     298599       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(114.7674491,4.1262496,3142.361,2968.750,0)
            useAV= True
            doDisDBSCAN.calDisByID(     584718       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################   AV GOOD
            lRange, bRange =         box(121.7015519,-2.1891454,2779.653,2176.837,0)
            useAV=True
            doDisDBSCAN.calDisByID(     252641       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL=  4  ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################  
            lRange, bRange =          box(145.5692114,-0.4866630,2478.245,2394.521,0)
            useAV= True
            doDisDBSCAN.calDisByID(     543968       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(113.9796604,-4.0426408,2790.816,2284.981,0)
            useAV= True
            doDisDBSCAN.calDisByID(      478223        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =         box(112.4043135,-4.2617778,3071.836,2950.707,0)
            useAV= True
            doDisDBSCAN.calDisByID(     153938       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =         box(139.1374911,-3.2004256,3293.163,2455.918,0)
            useAV= True
            doDisDBSCAN.calDisByID(     201475       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(148.2904832,3.7123739,2874.541,2476.850,0)
            useAV= True
            doDisDBSCAN.calDisByID(     521951       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =          box(138.9292117,-1.0156494,2654.066,2520.107,0)
            useAV= True
            doDisDBSCAN.calDisByID(     418641       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =         box(136.3957051,0.0754191,2846.536,2596.202,0)
            useAV= True
            doDisDBSCAN.calDisByID(     83738       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(124.9666702,2.9267485,3139.669,2518.712,0)
            useAV= True
            doDisDBSCAN.calDisByID(     67768       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(124.9884100,2.9038530,2972.220,2148.929,0)
            useAV= True
            doDisDBSCAN.calDisByID(     67768       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(105.5643392,2.4411318,2821.181,2922.454,0)
            useAV= True
            doDisDBSCAN.calDisByID(      282636       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL= 3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =           box(114.6314746,2.1312795,5833.333,4295.750,0)
            useAV= True
            doDisDBSCAN.calDisByID(     54231       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  , useMask=False,  cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(108.0412346,4.3453113,2934.028,3177.083,0)
            useAV= True
            doDisDBSCAN.calDisByID(     536955       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(112.9174813,-1.4926970,2976.393,2338.950,0)
            useAV= True
            doDisDBSCAN.calDisByID(      195541        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =          box(140.0333173,2.4936456,2972.220,2637.322,0)
            useAV= True
            doDisDBSCAN.calDisByID(      440703       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(113.6642173,-1.7557853,3415.959,3014.082,0)
            useAV= True
            doDisDBSCAN.calDisByID(      201598        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =         box(125.5379143,2.0571888,2652.392,2447.434,0)
            useAV= True
            doDisDBSCAN.calDisByID(     89065       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =         box(109.1603301,2.3066859,2895.472,2139.626,0)
            useAV= True
            doDisDBSCAN.calDisByID(     458149       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(120.4554732,3.6373211,3245.161,2541.876,0)
            useAV= True
            doDisDBSCAN.calDisByID(     176986       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(107.7919825,-2.3471949,2864.583,4861.111,0)
            useAV= True
            doDisDBSCAN.calDisByID(     87053       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(145.8683922,-3.8336873,2809.124,3291.377,0)
            useAV= True
            doDisDBSCAN.calDisByID(     499604       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL= 3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(137.0126228,0.2303893,3034.176,2712.674,0)
            useAV= True
            doDisDBSCAN.calDisByID(      57353        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(112.7611988,4.8579404,3038.194,2797.068,0)
            useAV= True
            doDisDBSCAN.calDisByID(      300185        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(121.5647184,-2.8377384,2965.856,2423.322,0)
            useAV= True
            doDisDBSCAN.calDisByID(     243636       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(111.5463384,3.2772098,2755.931,2483.827,0)
            useAV= True
            doDisDBSCAN.calDisByID(      266343       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(149.6105091,4.8378447,3653.067,2857.350,0)
            useAV= True
            doDisDBSCAN.calDisByID(     348422       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3  ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =          box(146.8622057,-0.6873981,2679.184,2379.171,0)
            useAV= False
            doDisDBSCAN.calDisByID(      462937        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(132.5962648,-1.6057277,2622.251,2361.031,0)
            useAV= True
            doDisDBSCAN.calDisByID(      176427        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(134.6669368,2.1853405,2160.092,3325.537,0)
            useAV= True
            doDisDBSCAN.calDisByID(      158643        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(106.7597231,1.6573190,2712.674,2478.245,0)
            useAV= True
            doDisDBSCAN.calDisByID(     261639       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =          box(141.5487879,0.2578133,3236.690,3215.278,0)
            useAV= True
            doDisDBSCAN.calDisByID(     529800       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### AV GOOD
            lRange, bRange =          box(112.1463454,-1.5355212,3194.927,3448.110,0)
            useAV= True
            doDisDBSCAN.calDisByID(      353634        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  , useMask=False,    cutDis= 1500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 1   )

            ####################
            lRange, bRange =          box(110.7714507,2.7066191,2644.299,2888.495,0)
            useAV= True
            doDisDBSCAN.calDisByID(      40595       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            lRange, bRange =          box(126.6189535,0.8942016,2567.939,2446.810,0)
            useAV= True
            doDisDBSCAN.calDisByID(     403167       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### AV GOOD
            lRange, bRange =          box(132.9798639,0.9859154,3356.481,2604.167,0)
            useAV= True
            doDisDBSCAN.calDisByID(     241763       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 1   )

            ####################
            lRange, bRange =   box(106.4146451,1.2365953,2471.547,2447.434,0)
            useAV= True
            doDisDBSCAN.calDisByID(   400830            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 1   )
            ####################
            lRange, bRange =   box(122.5547740,2.4537760,2905.240,2026.133,0)
            useAV= True
            doDisDBSCAN.calDisByID(  52227            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(108.3703898,-2.0597429,4479.260,3131.296,0)
            useAV= True
            doDisDBSCAN.calDisByID(      199299           , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   useMask=False, cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =  box(136.3338128,1.8078323,3923.611,2288.773,0)
            useAV= True
            doDisDBSCAN.calDisByID(  315503            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(112.3650489,2.3528879,3125.000,3055.556,0)
            useAV= True
            doDisDBSCAN.calDisByID(          23297     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=   4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### AV GDOOD
            lRange, bRange =    box(132.0123736,-2.5034317,4340.279,3269.676,0)
            useAV= True
            doDisDBSCAN.calDisByID(   262071            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  , useMask=False,  cutDis= 1500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 1   )

            ####################
            lRange, bRange =    box(114.6136268,-2.7530428,2951.288,3188.508,0)
            useAV= True
            doDisDBSCAN.calDisByID(  351099      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False,  cutDis= 1500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 1   )
            ####################
            lRange, bRange =   box(134.7111519,-1.9275349,2432.274,2282.074,0)
            useAV= False
            doDisDBSCAN.calDisByID(  328860            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(137.3628099,2.6800087,3062.307,3110.532,0)
            useAV= True
            doDisDBSCAN.calDisByID(    137342      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(138.6104300,-0.4935901,2746.163,2386.148,0)
            useAV= True
            doDisDBSCAN.calDisByID(   583614            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL= 4 ,    lExpand=0.5, bExpand=0.5, NL= 1   )
            ####################
            lRange, bRange =   box(124.1996663,1.4473995,2451.655,2054.351,0)
            useAV= True
            doDisDBSCAN.calDisByID(   288425             , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2500 , useBaseline=True, SL= 4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(138.4051097,2.3605630,2579.793,2578.291,0)
            useAV= True
            doDisDBSCAN.calDisByID(   165114            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=   4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(140.2646890,-0.9084870,3279.321,2724.730,0)
            useAV= True
            doDisDBSCAN.calDisByID(       465463        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =    box(110.3939861,-1.8372248,2517.357,2552.079,0)
            useAV= True
            doDisDBSCAN.calDisByID(    114580            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(147.9079403,-1.8967175,2436.383,2478.245,0)
            useAV= True
            doDisDBSCAN.calDisByID(   339759            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4  ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(134.6365230,2.1084393,3038.194,2763.310,0)
            useAV= True
            doDisDBSCAN.calDisByID(      547890          , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(120.4333144,2.8326542,2574.528,2239.630,0)
            useAV= True
            doDisDBSCAN.calDisByID(     45718           , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL= 4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(115.0521298,-0.5666536,2106.302,2052.467,0)
            useAV= True
            doDisDBSCAN.calDisByID(        310811        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(115.1154263,3.4388937,2338.927,2363.040,0)
            useAV= True
            doDisDBSCAN.calDisByID(         260581       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### AV GOOD
            lRange, bRange =     box(109.1132135,-0.1513716,3122.589,2109.857,0)
            useAV= True
            doDisDBSCAN.calDisByID(          465607     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL= 1  )

            ####################
            lRange, bRange =    box(131.4323696,3.2720911,2735.000,2637.322,0)
            useAV= True
            doDisDBSCAN.calDisByID(  449313            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =     box(120.4662752,0.2156344,2488.478,2116.369,0)
            useAV= True
            doDisDBSCAN.calDisByID(   114058             , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(142.0013345,-4.9140122,5505.723,2371.078,0)
            useAV= True
            doDisDBSCAN.calDisByID(   294554            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  , useMask=False,  cutDis= 3000 , useBaseline=True, SL= 3  ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =   box(116.0763026,0.8706275,2880.123,2704.301,0)
            useAV= True
            doDisDBSCAN.calDisByID(   530085            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=   4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(113.2563588,0.3837630,3420.139,3072.917,0)
            useAV= True
            doDisDBSCAN.calDisByID(      526955          , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  , useMask=False,   cutDis= 2000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange = box(127.0275227,0.9548469,2693.138,2204.745,0)
            useAV= True
            doDisDBSCAN.calDisByID(          332559       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(109.7195636,-1.4837301,2921.985,2419.638,0)
            useAV= True
            doDisDBSCAN.calDisByID(      112639          , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(113.6043513,0.5961073,2570.342,2294.051,0)
            useAV= True
            doDisDBSCAN.calDisByID(   302330     ,  useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(142.6920142,0.3590988,2853.331,2722.721,0)
            useAV= True
            doDisDBSCAN.calDisByID(   133331            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################  AV GOOD
            lRange, bRange =     box(134.8962161,1.6825291,1821.008,3510.033,0)
            useAV= True
            doDisDBSCAN.calDisByID(   148021            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3,    lExpand=0.5, bExpand=0.5, NL= 1   )

            ####################
            lRange, bRange =    box(133.2032316,-0.3859139,2777.896,2769.821,0)
            useAV= True
            doDisDBSCAN.calDisByID(   161598            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 1   )
            ####################
            lRange, bRange =    box(134.5831769,1.8184924,2310.796,4189.574,0)
            useAV= True
            doDisDBSCAN.calDisByID(        216222       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  , useMask=False,   cutDis= 3500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(107.3738491,3.0368299,2616.390,2279.167,0)
            useAV= True
            doDisDBSCAN.calDisByID(   401308            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =    box(124.3859233,0.2279353,2654.066,2980.592,0)
            useAV= True
            doDisDBSCAN.calDisByID(          535520     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL= 3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(122.0991400,4.1498265,3008.602,1884.232,0)
            useAV= True
            doDisDBSCAN.calDisByID(        47672       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(108.8029202,0.5560662,2913.612,2602.157,0)
            useAV= True
            doDisDBSCAN.calDisByID(     457600           , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =       box(110.8890003,-4.1507600,3920.718,3153.935,0)
            useAV= True
            doDisDBSCAN.calDisByID(        459786       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  , useMask=False,  cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(122.5887782,0.4161591,2923.251,1942.105,0)
            useAV= True
            doDisDBSCAN.calDisByID(      143023         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL= 4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(111.3993701,4.8421342,2295.447,2225.676,0)
            useAV= True
            doDisDBSCAN.calDisByID(      333753         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(133.1021650,-1.5832942,2672.486,2330.890,0)
            useAV= True
            doDisDBSCAN.calDisByID(      269760            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(143.2311012,-1.4897749,2604.167,2712.674,0)
            useAV= True
            doDisDBSCAN.calDisByID(     145431           , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(137.2175774,-0.8351548,2610.576,2534.992,0)
            useAV= True
            doDisDBSCAN.calDisByID(        391430       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(111.9165767,3.2741979,3677.180,2461.500,0)
            useAV= True
            doDisDBSCAN.calDisByID(       518431         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(130.8681946,4.2492031,5562.500,5145.833,0)
            useAV= True
            doDisDBSCAN.calDisByID(      270455         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  , useMask=False,  cutDis= 3000 , useBaseline=True, SL=  6,    lExpand=0.5, bExpand=0.5, NL=  1   )
            ####################
            lRange, bRange =    box(116.3672710,3.0518151,2534.023,2524.332,0)
            useAV= True
            doDisDBSCAN.calDisByID(      354965         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3  ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### AV GOOD
            lRange, bRange =      box(108.7108229,-0.6949806,3082.528,3192.163,0)
            useAV= True
            doDisDBSCAN.calDisByID(      325988         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False,  cutDis= 2500 , useBaseline=True, SL= 6 ,    lExpand=0.5, bExpand=0.5, NL= 0.5  )

            ####################
            lRange, bRange =    box(121.8255605,-0.9624780,3243.152,2350.984,0)
            useAV= True
            doDisDBSCAN.calDisByID(         133052      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### AV GOOD
            lRange, bRange =    box(121.4695202,0.3812644,3211.806,2835.649,0)
            useAV= True
            doDisDBSCAN.calDisByID(         418923      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 1   )

            ####################
            lRange, bRange =    box(106.5175711,0.9531791,2486.615,2319.166,0)
            useAV= True
            doDisDBSCAN.calDisByID(   94029            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################AV GOOD
            lRange, bRange =     box(124.6312621,-1.9743158,2446.810,2437.119,0)
            useAV= True
            doDisDBSCAN.calDisByID(   305787            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL=4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =    box(115.2570387,3.6793458,2672.207,2351.263,0)
            useAV= True
            doDisDBSCAN.calDisByID(       71479        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(118.7615343,-1.3832421,2267.538,1904.151,0)
            useAV= True
            doDisDBSCAN.calDisByID(        161487       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =    box(116.6782306,4.3508955,2494.292,2087.298,0)
            useAV= True
            doDisDBSCAN.calDisByID(        177823        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(114.6317744,3.3605789,2277.229,2127.029,0)
            useAV= True
            doDisDBSCAN.calDisByID(   381940            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(110.6189639,2.3726265,2383.822,2068.887,0)
            useAV= True
            doDisDBSCAN.calDisByID(      213466         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(141.7785430,1.5799965,2762.909,2451.454,0)
            useAV= True
            doDisDBSCAN.calDisByID(         123804       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(120.2766436,2.1313927,2348.937,2314.052,0)
            useAV= True
            doDisDBSCAN.calDisByID(   107306     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(136.1886099,0.8055306,2802.445,2418.708,0)
            useAV= True
            doDisDBSCAN.calDisByID(        33759       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(138.2995034,-0.6316176,2465.221,2691.975,0)
            useAV= True
            doDisDBSCAN.calDisByID(          165675     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL= 3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(140.5726093,0.3208530,3362.934,2602.436,0)
            useAV= True
            doDisDBSCAN.calDisByID(          138104     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =      box(144.9085312,2.5694193,2561.970,2300.749,0)
            useAV= False
            doDisDBSCAN.calDisByID(        558095       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(140.6558105,-0.2265989,3327.546,2488.426,0)
            useAV= True
            doDisDBSCAN.calDisByID(       612846        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL= 3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### AV GOOD
            lRange, bRange =    box(130.0402813,-4.4565050,2511.735,2609.413,0)
            useAV= True
            doDisDBSCAN.calDisByID(          149173     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(138.3412118,-3.5733296,2813.143,2551.923,0)
            useAV = True
            doDisDBSCAN.calDisByID(  400048     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, useMask=False, cutDis=2500,   useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=1)

            ####################
            lRange, bRange =    box(137.0693877,1.8309352,2700.617,2604.167,0)
            useAV= True
            doDisDBSCAN.calDisByID(      576249         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =     box(146.8938225,-1.3491265,3156.413,2478.245,0)
            useAV= True
            doDisDBSCAN.calDisByID(   595842             , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(130.9135757,-0.9195770,2511.735,2302.424,0)
            useAV= True
            doDisDBSCAN.calDisByID(     195596          , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(119.5883663,3.2779036,2239.630,2197.768,0)
            useAV= False
            doDisDBSCAN.calDisByID(      128566          , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(119.5883663,3.2779036,2239.630,2197.768,0)
            useAV= True
            doDisDBSCAN.calDisByID(      128566          , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =    box(124.4914718,2.4387596,2902.449,2323.355,0)
            useAV= True
            doDisDBSCAN.calDisByID(      5580         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(119.6678176,-0.5110542,2458.333,2375.000,0)
            useAV= True
            doDisDBSCAN.calDisByID(    106817            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL= 4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(131.1794558,-0.7955527,2679.184,2662.439,0)
            useAV= True
            doDisDBSCAN.calDisByID(     260933           , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =    box(143.8713025,0.5011034,3003.472,2986.111,0)
            useAV= True
            doDisDBSCAN.calDisByID(      596174         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(124.3994517,2.0823998,2452.058,2495.548,0)
            useAV= True
            doDisDBSCAN.calDisByID(          42178     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL= 4  ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =  box(115.1207900,3.8908750,2525.689,2260.561,0)
            useAV= True
            doDisDBSCAN.calDisByID(        345863       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### AG GOOD
            lRange, bRange =            box(143.8073971,-3.3211090,5714.699,4533.179,0)
            useAV= False
            doDisDBSCAN.calDisByID(     186740      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False  , useMask=False, cutDis= 2500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 1   )

            lRange, bRange =    box(108.7217953,3.6185217,3290.373,3148.041,0)
            useAV= True
            doDisDBSCAN.calDisByID(   348047            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################  AV GOOD
            lRange, bRange =    box(122.5996191,-2.3075627,2476.850,2267.538,0)
            useAV= True
            doDisDBSCAN.calDisByID(   250260            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =    box(128.4314770,1.7095680,2030.125,2122.183,0)
            useAV= False
            doDisDBSCAN.calDisByID(  113604            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(145.8631823,-1.9342404,2343.123,2180.325,0)
            useAV= True
            doDisDBSCAN.calDisByID(    227167            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(137.7375852,0.9918984,2465.221,2122.183,0)
            useAV= True
            doDisDBSCAN.calDisByID(      513929         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =   box(146.0327693,1.8084771,2277.229,2103.610,0)
            useAV= True
            doDisDBSCAN.calDisByID(   125216            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(114.3052591,3.2020661,2742.814,2150.045,0)
            useAV= True
            doDisDBSCAN.calDisByID(  277984            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(134.4868783,0.1127507,2114.043,2190.791,0)
            useAV= True
            doDisDBSCAN.calDisByID(  258214            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =  box(131.2217699,-0.5849484,2372.194,1837.288,0)
            useAV= True
            doDisDBSCAN.calDisByID(  280221            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### AV GOOD
            lRange, bRange =  box( 131.2217699,-0.5849484,2372.194,1837.288,0)
            useAV= True
            doDisDBSCAN.calDisByID(  280221            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =  box(116.1195235,-0.5829809,2216.664,1970.368,0)
            useAV= True
            doDisDBSCAN.calDisByID(  126794            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =  box(121.4280388,1.2718033,2069.856,2017.528,0)
            useAV= True
            doDisDBSCAN.calDisByID(  208652            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =    box(113.0772247,3.7082410,2407.079,2127.998,0)
            useAV= False
            doDisDBSCAN.calDisByID(    333454     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =    box(113.0772247,3.7082410,2407.079,2127.998,0)
            useAV= True
            doDisDBSCAN.calDisByID(    333454     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =    box(140.9094016,3.8978258,2511.735,2323.355,0)
            useAV= True
            doDisDBSCAN.calDisByID(  417207            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(131.5728743,-3.8019671,2628.949,2268.934,0)
            useAV= True
            doDisDBSCAN.calDisByID(  559492            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(145.4322103,0.2803350,2131.874,2102.803,0)
            useAV= True
            doDisDBSCAN.calDisByID(  165866            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =  box(137.2702439,4.9184988,2592.110,2290.702,0)
            useAV= True
            doDisDBSCAN.calDisByID(  495083            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(145.3165719,2.7611469,2244.282,2040.785,0)
            useAV= True
            doDisDBSCAN.calDisByID(  299887            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2500 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(114.5806225,3.6059490,2393.125,2372.194,0)
            useAV= True
            doDisDBSCAN.calDisByID(     474580        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =     box(147.2889011,-2.0483345,2670.811,2185.209,0)
            useAV= True
            doDisDBSCAN.calDisByID(  618012            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=   2 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(143.3570496,-0.7067565,2406.433,2394.320,0)
            useAV= True
            doDisDBSCAN.calDisByID(  509336            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =      box(124.8649199,0.5670564,2832.679,2707.092,0)
            useAV= True
            doDisDBSCAN.calDisByID(  329672            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL= 4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(128.6753848,1.8263605,2034.505,2528.480,0)
            useAV= True
            doDisDBSCAN.calDisByID(         17926       , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL= 4,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =   box(139.9738922,3.5363748,2419.638,2227.072,0)
            useAV= True
            doDisDBSCAN.calDisByID(  584580            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### AV GOOD
            lRange, bRange =    box(124.1429992,0.5965697,2447.779,1755.889,0)
            useAV= True
            doDisDBSCAN.calDisByID(  293831            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =     box(118.1872051,-0.3034908,2986.174,2107.066,0)
            useAV= True
            doDisDBSCAN.calDisByID(  564582            , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL= 3  ,    lExpand=0.5, bExpand=0.5, NL= 1   )
            ####################
            lRange, bRange =  box(145.1304596,-3.8096847,2611.545,2349.906,0)
            useAV= True
            doDisDBSCAN.calDisByID(    397124    , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =     box(147.0143527,-1.3864550,2627.610,2461.500,0)
            useAV= True
            doDisDBSCAN.calDisByID(        310570     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################
            lRange, bRange =   box(105.2274027,1.4234433,2323.355,2002.411,0)
            useAV= True
            doDisDBSCAN.calDisByID(    342750        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(139.6020665,-0.0358123,2903.565,2089.763,0)
            useAV= True
            doDisDBSCAN.calDisByID(    526779        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL= 4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(143.3460058,0.2120669,2402.893,2118.230,0)
            useAV= True
            doDisDBSCAN.calDisByID(    152174        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(128.3746131,0.4674426,2679.184,2285.679,0)
            useAV= False
            doDisDBSCAN.calDisByID(    286475        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(115.3848118,3.1112439,2558.248,2186.140,0)
            useAV= True
            doDisDBSCAN.calDisByID(    248198        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(120.0931019,-0.2321225,2411.265,2230.421,0)
            useAV= True
            doDisDBSCAN.calDisByID(    108483        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(143.8169828,4.7318129,3167.577,1932.640,0)
            useAV= True
            doDisDBSCAN.calDisByID(    537118        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL= 3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(126.3226243,-0.2382207,2511.735,2343.123,0)
            useAV= False
            doDisDBSCAN.calDisByID(        79416      , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### AV GOOD
            lRange, bRange =   box(128.3482152,0.9011521,2165.790,1446.832,0)
            useAV= True
            doDisDBSCAN.calDisByID(    296887        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 1   )

            #################### AV GDOOD
            lRange, bRange =     box(122.7525255,0.7395408,2193.582,2227.072,0)
            useAV= True
            doDisDBSCAN.calDisByID(    278640        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL= 1  ,    lExpand=0.5, bExpand=0.5, NL=  1    )


            ####################
            lRange, bRange =    box(142.4949200,4.0032407,2383.822,2214.242,0)
            useAV= True
            doDisDBSCAN.calDisByID(    324590    , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2000 , useBaseline=True, SL=  2  ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(128.8253129,1.3132500,2592.110,2511.735,0)
            useAV= True
            doDisDBSCAN.calDisByID(    17205        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL= 2 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(108.6179019,4.6337694,2552.434,2343.123,0)
            useAV= True
            doDisDBSCAN.calDisByID(    548833        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  1.5 ,    lExpand=0.5, bExpand=0.5, NL= 1   )
            ####################
            lRange, bRange =   box(149.5714798,-2.1908241,2401.218,2391.172,0)
            useAV= True
            doDisDBSCAN.calDisByID(    183865        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(140.8900045,0.3056710,2226.839,1971.014,0)
            useAV= True
            doDisDBSCAN.calDisByID(    279472        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  2 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(110.3830812,3.3909948,2301.455,2008.726,0)
            useAV= True
            doDisDBSCAN.calDisByID(    39879        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(115.8152661,3.7914968,2267.538,2136.719,0)
            useAV= True
            doDisDBSCAN.calDisByID(    54791        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3  ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(149.0337759,-3.0022991,3159.049,2131.874,0)
            useAV= True
            doDisDBSCAN.calDisByID(     244775        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1500  , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 1   )
            ####################
            lRange, bRange =   box(124.3089473,-0.8429343,2818.055,2559.560,0)
            useAV= True
            doDisDBSCAN.calDisByID(     349136        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False,  cutDis= 2000 , useBaseline=True, SL=  6 ,    lExpand=0.5, bExpand=0.5, NL= 1   )
            ####################
            lRange, bRange =   box(110.2803589,4.0106815,2055.159,1752.336,0)
            useAV= True
            doDisDBSCAN.calDisByID(     503250         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(113.4095047,1.1435205,2401.265,1889.615,0)
            useAV= True
            doDisDBSCAN.calDisByID(    293924        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(137.0141413,0.4275531,2010.745,1855.699,0)
            useAV= True
            doDisDBSCAN.calDisByID(    303594        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2500 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(148.2299073,-2.3851554,2950.707,2374.132,0)
            useAV= True
            doDisDBSCAN.calDisByID(    205812        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(127.9325645,1.1323643,2595.459,2393.125,0)
            useAV= True
            doDisDBSCAN.calDisByID(    70419        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL= 4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =  box(135.0355972,-0.0099357,2520.721,1954.330,0)
            useAV= True
            doDisDBSCAN.calDisByID(     182540         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(135.2321316,-0.0452188,2093.112,1768.486,0)
            useAV= True
            doDisDBSCAN.calDisByID(     129630         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(114.5698073,-1.5604702,2160.137,1715.997,0)
            useAV= True
            doDisDBSCAN.calDisByID(     421478         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  5,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(125.3105061,0.7670292,2354.751,1870.235,0)
            useAV= True
            doDisDBSCAN.calDisByID(        89598     , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =    box(117.1530713,5.0015267,2608.315,1655.432,0)
            useAV= True
            doDisDBSCAN.calDisByID(     40683        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =     box(120.7283209,-0.9163541,2197.768,2009.388,0)
            useAV= True
            doDisDBSCAN.calDisByID(     314778        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 1000 , useBaseline=True, SL= 3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(145.2057240,-0.0981502,2441.964,2100.089,0)
            useAV= True
            doDisDBSCAN.calDisByID(     193202        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(142.7741198,0.0846830,2453.593,2011.714,0)
            useAV= True
            doDisDBSCAN.calDisByID(     296761         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(139.9834070,2.5788121,2150.604,1511.311,0)
            useAV= True
            doDisDBSCAN.calDisByID(     561991        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(130.6942393,3.2330432,2208.589,1841.164,0)
            useAV= True
            doDisDBSCAN.calDisByID(    545135        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(112.1393786,2.2589222,2228.777,2059.196,0)
            useAV= True
            doDisDBSCAN.calDisByID(     100064         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(146.2827927,-2.3573739,2287.996,1776.561,0)
            useAV= True
            doDisDBSCAN.calDisByID(    619326        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL= 2 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(118.3758282,-0.5251898,2560.574,2428.010,0)
            useAV= True
            doDisDBSCAN.calDisByID(    488606        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  3 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(134.2042824,-0.9459448,2253.003,1990.556,0)
            useAV= True
            doDisDBSCAN.calDisByID(     172989         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(134.6879063,5.0813641,2907.101,1380.873,0)
            useAV= True
            doDisDBSCAN.calDisByID(     258967         , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(129.1325854,0.6420529,2261.724,1930.315,0)
            useAV= True
            doDisDBSCAN.calDisByID(     175849        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=   4  ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(122.1918231,-0.1215179,2813.143,2079.716,0)
            useAV= True
            doDisDBSCAN.calDisByID(    96754        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ####################
            lRange, bRange =   box(118.6534530,2.5256047,2821.515,2277.306,0)
            useAV= True
            doDisDBSCAN.calDisByID(     49336        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  4 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### Largest cloud RX7
            lRange, bRange =       box(148.1588470,2.6170538,14687.500,14427.083,0)
            useAV= False
            doDisDBSCAN.calDisByID(      227384        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False, cutDis= 500  , useBaseline=True, SL= 10 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### Largest cloud RX8
            lRange, bRange =      box(143.9353241,4.0589674,9563.079,8015.046,0)
            useAV= False
            doDisDBSCAN.calDisByID(      227384        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False, cutDis= 1500  , useBaseline=True, SL= 10 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### Largest cloud   RX9
            lRange, bRange =     box(134.6871099,0.2076384,5726.755,4963.188,0)
            useAV= False
            doDisDBSCAN.calDisByID(      227384        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False, cutDis= 1500  , useBaseline=True, SL= 10 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### Largest cloud   RX10
            lRange, bRange =    box(149.4269647,-1.5103541,4398.892,6453.813,0)
            useAV= False
            doDisDBSCAN.calDisByID(      227384        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False, cutDis= 1500  , useBaseline=True, SL= 10 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### Largest cloud   RX11
            lRange, bRange =    box(124.8339230,-0.8728112,4036.458,3732.639,0)
            useAV= False
            doDisDBSCAN.calDisByID(      227384        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False, cutDis= 1500  , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### Largest cloud   RX12
            lRange, bRange =    box(112.5211329,-2.8129080,14998.071,5774.981,0)
            useAV= False
            doDisDBSCAN.calDisByID(      227384        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False, cutDis= 800  , useBaseline=True, SL=5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### Largest cloud   RX13
            lRange, bRange =    box(105.6950933,-1.3000795,4360.651,4197.853,0)
            useAV= False
            doDisDBSCAN.calDisByID(      227384        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False, cutDis= 1500  , useBaseline=True, SL=5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            #################### Largest cloud   RX 14
            lRange, bRange =     box(107.1700387,-0.1954013,3410.259,3296.801,0)
            useAV= False
            doDisDBSCAN.calDisByID(      227384        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False, cutDis= 1500  , useBaseline=True, SL= 6 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### Largest cloud   RX15
            lRange, bRange =      box(121.2507312,-3.8218493,6540.557,2692.579,0)
            useAV= False
            doDisDBSCAN.calDisByID(      227384        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False, cutDis= 1500  , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            #################### Largest cloud  RX17
            lRange, bRange =       box(110.7551618,3.8044724,4667.928,4477.947,0)
            useAV= False
            doDisDBSCAN.calDisByID(      227384        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False, cutDis= 1500  , useBaseline=True, SL= 6 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### Largest cloud   RX18
            lRange, bRange =    box(134.0244930,-4.0537287,3906.250,3342.014,0)
            useAV= True
            doDisDBSCAN.calDisByID(      227384        , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  useMask=False, cutDis= 2000  , useBaseline=True, SL= 8 ,    lExpand=0.5, bExpand=0.5, NL= 2   )

            ####################   AG Candidates
            lRange, bRange  =  box(114.7876086,4.3273689,4955.473,6558.521,0)
            useAV= False
            doDisDBSCAN.calDisByID(    222339        ,  useForegroundFITS=True,  useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )


            #################### AV GOOD
            lRange, bRange =   box(143.1263942,0.7612341,11611.111,12277.199,0)
            useAV= True
            doDisDBSCAN.calDisByID(    454540        ,    useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 500 , useBaseline=True, SL= 3 ,    lExpand=0.5, bExpand=0.5, NL= 1   )

            ###################
            lRange, bRange =    box(116.2567739,4.0607890,3790.509,3573.495,0)
            useAV= False
            doDisDBSCAN.calDisByID(         222339          ,   useForegroundFITS=True ,  useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 3500 , useBaseline=True, SL=  5 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ###################
            lRange, bRange =    box(133.5449203,0.7156154,7754.630,6640.625,0)
            useAV= False
            doDisDBSCAN.calDisByID(        89568  ,   useForegroundFITS=True,    useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,  lowerDisCut=1000,  cutDis= 3000 , useBaseline=True, SL=  10 ,    lExpand=0.5, bExpand=0.5, NL= 2   )
            ###################
            lRange, bRange =   box(145.8254079,1.3522415,3368.056,3315.972,0)
            useAV= False
            doDisDBSCAN.calDisByID(     160462      ,    useForegroundFITS=True,   useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False  ,   cutDis= 2500 , useBaseline=True, SL= 5 ,    lExpand=0.5, bExpand=0.5, NL= 1  )
            ###################
            lRange, bRange =  box(148.5220859,-0.1596938,11887.539,5160.108,0) 
            useAV= False
            doDisDBSCAN.calDisByID(     181901    ,     cutDis= 3500 ,  SL=  5 ,  NL= 2  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=True,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ################### AG GOOD
            lRange, bRange =   box(109.8285909,1.7517648,7656.250,10329.861,0)
            useAV= False
            doDisDBSCAN.calDisByID(     458317     ,     cutDis= 500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=True,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ###################
            lRange, bRange =   box(144.0426486,-1.5442145,5867.991,3805.974,0)
            useAV= False
            doDisDBSCAN.calDisByID(     181719     ,     cutDis= 3500 ,  SL=  6 ,  NL= 2  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=True,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ################### AG GOOD
            lRange, bRange =    box(143.5637502,-1.7532601,4834.587,6667.149,0)
            useAV= False
            doDisDBSCAN.calDisByID(     258072     ,     cutDis= 1500 ,  SL=  5,  NL= 1  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=True,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =   box(138.7868938,1.9423974,5280.671,5353.009,0)
            useAV= True
            doDisDBSCAN.calDisByID(     158599    ,     cutDis= 3500 ,  SL=  5 ,  NL= 2  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=True,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =   box(136.7715826,1.4882143,5702.643,5184.221,0)
            useAV= False
            doDisDBSCAN.calDisByID(     139159     ,     cutDis= 3500 ,  SL=  4  ,  NL= 1  ,  lowerDisCut=1500,   lExpand=0.5, bExpand=0.5,   useForegroundFITS=True,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ################### AG GOOD
            lRange, bRange =   box(144.5320592,-2.5559178,5786.480,4921.203,0)
            useAV= False
            doDisDBSCAN.calDisByID(     247603     ,     cutDis= 1100 ,  SL=  4 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=True,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ################### AV GOOD
            lRange, bRange =   box(131.5835131,-0.1102546,2413.194,3784.722,0)
            useAV= True
            doDisDBSCAN.calDisByID(     267299     ,     cutDis= 1500 ,  SL= 4,  NL= 2  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=True,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ################### AG GOOD
            lRange, bRange =   box(129.7549294,-0.0526735,4591.451,5405.253,0)
            useAV= False
            doDisDBSCAN.calDisByID(     523676     ,     cutDis= 1500 ,  SL=  5 ,  NL= 2  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=False,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =   box(143.8264570,-3.4325486,5399.306,3541.667,0)
            useAV= False
            doDisDBSCAN.calDisByID(     186740     ,  lowerDisCut=100,   cutDis= 2500 ,  SL=  6 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=True,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  True       )
            ###################
            lRange, bRange =   box(111.6603823,4.0211684,4615.163,4716.436,0)
            useAV= True
            doDisDBSCAN.calDisByID(     209687    ,     cutDis= 3500 ,  SL=  8 ,  NL= 2  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=True,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =   box(114.1041911,-2.5178716,7447.917,6562.500,0)
            useAV= True
            doDisDBSCAN.calDisByID(     540708     ,     cutDis= 500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=False,  useForegroundFITS=False,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =   box(117.4933255,-0.6279153,4617.573,4062.982,0)
            useAV= True
            doDisDBSCAN.calDisByID(     505158     ,     cutDis= 500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=False,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =   box(138.5404139,1.3787931,4751.365,2448.941,0)
            useAV= False
            doDisDBSCAN.calDisByID(    422095    ,     cutDis= 1500 ,  SL=  4 ,  NL= 2  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=False,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ################### AV GOOD
            lRange, bRange =   box(134.0338977,-4.0083215,3948.447,3375.772,0)
            useAV= True
            doDisDBSCAN.calDisByID(       270697      ,     cutDis= 1500 ,  SL= 5   ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=False,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =    box(124.0153097,0.9017730,2197.768,1901.244,0)
            useAV= True
            doDisDBSCAN.calDisByID(     311159     ,     cutDis= 3500 ,  SL=  4 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =    box(124.0153097,0.9017730,2197.768,1901.244,0)
            useAV= True
            doDisDBSCAN.calDisByID(     311159     ,     cutDis= 2000 ,  SL=  3  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(132.3786706,-2.3608631,2441.964,2065.204,0)
            useAV= True
            doDisDBSCAN.calDisByID(     283247     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,   useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(138.2130924,-3.5370802,3090.829,2895.472,0)
            useAV= True
            doDisDBSCAN.calDisByID(    421203    ,     cutDis= 1500 ,  SL=  4,  NL= 1  ,    lExpand=0.5, bExpand=0.5, useMask=False,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =  box(132.0536543,-2.4878083,4412.281,2821.515,0)
            useAV = True
            doDisDBSCAN.calDisByID(283247, cutDis= 3000, SL= 5 , NL=1, lExpand=0.5, bExpand=0.5, useForegroundFITS=True, useMask=False,
                                   useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False)
            ###################
            lRange, bRange =     box(132.4222895,-1.1583397,2034.505,1917.291,0)
            useAV= True
            doDisDBSCAN.calDisByID(    287817    ,     cutDis= 3000 ,  SL= 4 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True,  useForegroundFITS=True ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ################### AV GOOD
            lRange, bRange =    box(130.7095659,1.3864350,5416.975,2402.893,0)
            useAV= True
            doDisDBSCAN.calDisByID(     335449     ,     cutDis= 1500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=False ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(142.2654297,-0.5757086,1833.759,1309.428,0)
            useAV= True
            doDisDBSCAN.calDisByID(     550120     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(135.3430283,1.0386993,2133.218,1988.536,0)
            useAV= True
            doDisDBSCAN.calDisByID(    124499    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(110.3597286,1.7129652,2172.250,1837.126,0)
            useAV= False
            doDisDBSCAN.calDisByID(    76661    ,     cutDis= 3500 ,  SL=  4 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS= True ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(117.3039870,-1.8053979,2386.148,2076.367,0)
            useAV= True
            doDisDBSCAN.calDisByID(    516615    ,     cutDis= 1500 ,  SL=   3  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(134.1402711,1.8810063,1901.055,1682.350,0)
            useAV= True
            doDisDBSCAN.calDisByID(     145824     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(127.1830273,3.7871027,2928.315,1849.599,0)
            useAV= True
            doDisDBSCAN.calDisByID(     573915     ,     cutDis= 3500 ,  SL=  4 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ################### AV GOOD
            lRange, bRange =     box(105.2099600,4.9838862,2117.338,1865.389,0)
            useAV= True
            doDisDBSCAN.calDisByID(     355655     ,     cutDis= 1500 ,  SL=  3  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS= False  ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(137.7215650,2.0920862,2035.643,1554.491,0)
            useAV= True
            doDisDBSCAN.calDisByID(     180631    ,     cutDis= 3500 ,  SL=  4 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange = box(139.1409895,-1.4517495,2312.871,1701.611,0)
            useAV= True
            doDisDBSCAN.calDisByID(     603525     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(122.9798474,-3.4779906,2139.626,2098.927,0)
            useAV= True
            doDisDBSCAN.calDisByID(     180086     ,     cutDis= 4000 ,  SL= 3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(120.4833269,3.8430123,3340.607,2570.342,0)
            useAV= True
            doDisDBSCAN.calDisByID(     185469     ,     cutDis= 3500 ,  SL=   3  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=True ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(149.9905296,-2.6513171,1883.655,1976.680,0)
            useAV= True
            doDisDBSCAN.calDisByID(    231253    ,     cutDis= 3500 ,  SL= 3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(119.7125590,-0.0305118,2397.348,1783.291,0)
            useAV= True
            doDisDBSCAN.calDisByID(    125587    ,     cutDis= 3500 ,  SL=  3  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(112.8639845,3.2584642,1931.338,1887.597,0)
            useAV= True
            doDisDBSCAN.calDisByID(     343332    ,     cutDis= 2000 ,  SL=  4  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(136.4687981,0.7577725,2223.932,1850.854,0)
            useAV= True
            doDisDBSCAN.calDisByID(     95246     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(131.8371757,4.6966123,2248.965,1643.319,0)
            useAV= True
            doDisDBSCAN.calDisByID(     288987    ,     cutDis= 3500 ,  SL=  3  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(147.4923417,-2.1077544,2058.227,1889.615,0)
            useAV= True
            doDisDBSCAN.calDisByID(     568048     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(134.5618618,2.2400642,2114.043,1848.916,0)
            useAV= True
            doDisDBSCAN.calDisByID(     212385     ,     cutDis= 3500 ,  SL=  3  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(123.8317725,-1.1104701,2461.345,2131.874,0)
            useAV= True
            doDisDBSCAN.calDisByID(     413188     ,     cutDis= 3500 ,  SL= 3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(127.4779533,3.6079864,2378.170,1602.943,0)
            useAV= True
            doDisDBSCAN.calDisByID(     195947     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS= True,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(120.6052836,0.7042758,2030.035,1727.212,0)
            useAV= True
            doDisDBSCAN.calDisByID(     428353     ,     cutDis= 2000 ,  SL=   4  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(119.1755718,-0.8207475,2486.617,2227.072,0)
            useAV= False
            doDisDBSCAN.calDisByID(     182414     ,     cutDis= 3500 ,  SL=  4 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(138.0403388,0.7732298,1850.585,1557.856,0)
            useAV= True
            doDisDBSCAN.calDisByID(     174879     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(121.2925001,1.6552244,2030.126,1826.629,0)
            useAV= True
            doDisDBSCAN.calDisByID(    538971    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(131.8868268,-0.8697651,2497.781,2169.860,0)
            useAV= True
            doDisDBSCAN.calDisByID(     460231     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(147.3970489,-2.1964202,1860.545,3120.288,0)
            useAV= True
            doDisDBSCAN.calDisByID(     177990     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=False ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ###################
            lRange, bRange =    box(135.7068587,0.2571372,2111.685,1881.540,0)
            useAV= True
            doDisDBSCAN.calDisByID(    186104     ,     cutDis= 3500 ,  SL=   3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(142.6157556,-1.6518973,2164.174,1893.652,0)
            useAV= True
            doDisDBSCAN.calDisByID(     265790     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(124.9555141,-3.1571011,2116.369,1907.058,0)
            useAV= True
            doDisDBSCAN.calDisByID(     165377     ,     cutDis= 3500 ,  SL=3  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =    box(142.4954162,-4.7907115,3203.572,2688.531,0)
            useAV= True
            doDisDBSCAN.calDisByID(    246162    ,     cutDis= 3500 ,  SL= 3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(150.0386884,1.9158297,1518.152,2313.567,0)
            useAV= True
            doDisDBSCAN.calDisByID(     601009     ,     cutDis= 3500 ,  SL=   3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(118.7122324,2.5923942,2581.504,2279.166,0)
            useAV= True
            doDisDBSCAN.calDisByID(     385539     ,     cutDis= 3500 ,  SL=  3  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(147.5399096,-2.7799143,2111.685,1796.750,0)
            useAV= True
            doDisDBSCAN.calDisByID(     236960     ,     cutDis= 3500 ,  SL=  3  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(118.3985271,3.7095605,2434.986,2309.399,0)
            useAV= True
            doDisDBSCAN.calDisByID(     125934     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(123.4357541,2.6759571,2122.427,2008.859,0)
            useAV= False
            doDisDBSCAN.calDisByID(     489857     ,     cutDis= 3500 ,  SL=  4 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(123.4357541,2.6759571,2122.427,2008.859,0)
            useAV= True
            doDisDBSCAN.calDisByID(     489857     ,     cutDis=  1200 ,  SL=  2 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(131.1693633,2.9058214,1901.055,1877.502,0)
            useAV= True
            doDisDBSCAN.calDisByID(     566170     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(141.0732182,-4.4476549,2114.153,1609.448,0)
            useAV= True
            doDisDBSCAN.calDisByID(     291044     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(141.0732182,-4.4476549,2114.153,1609.448,0)
            useAV= True
            doDisDBSCAN.calDisByID(     291044     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(117.8082480,-1.2545180,2156.100,2049.506,0)
            useAV= True
            doDisDBSCAN.calDisByID(     540970     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(128.5552317,-2.8947619,3625.271,3282.000,0)
            useAV= True
            doDisDBSCAN.calDisByID(     592287     ,     cutDis= 3500 ,  SL=  3  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(128.5552317,-2.8947619,3625.271,3282.000,0)
            useAV= True
            doDisDBSCAN.calDisByID(     592287     ,     cutDis= 3500 ,  SL=  3  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(143.2389624,-2.7236664,2065.926,2035.643,0)
            useAV= True
            doDisDBSCAN.calDisByID(     175419     ,     cutDis= 3500 ,  SL=  3  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(143.2389624,-2.7236664,2065.926,2035.643,0)
            useAV= True
            doDisDBSCAN.calDisByID(     175419     ,     cutDis= 3500 ,  SL=  3  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(142.5347316,1.5143380,2393.125,2079.158,0)
            useAV= True
            doDisDBSCAN.calDisByID(     536017    ,     cutDis= 3500 ,  SL=   2 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(106.4659453,0.3848808,2107.646,1792.710,0)
            useAV= True
            doDisDBSCAN.calDisByID(     85052     ,     cutDis= 3500 ,  SL=  4 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(113.5033452,-1.1771599,3375.772,3546.570,0)
            useAV= True
            doDisDBSCAN.calDisByID(    424583     ,     cutDis= 3500 ,  SL= 4 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=False ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(128.3433555,-3.6279349,2223.932,2112.493,0)
            useAV= True
            doDisDBSCAN.calDisByID(     590616 ,     cutDis= 3500 ,  SL=  2 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(107.1386764,5.0863670,1827.985,1555.880,0)
            useAV= True
            doDisDBSCAN.calDisByID(    609956    ,     cutDis= 1000 ,  SL= 2   ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(114.3767078,3.2429928,2123.798,1760.411,0)
            useAV= True
            doDisDBSCAN.calDisByID(    45269    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(143.4344756,-0.8534912,2337.309,2261.724,0)
            useAV= True
            doDisDBSCAN.calDisByID(    196857    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(128.7376390,-0.9484517,2485.571,1792.712,0)
            useAV= True
            doDisDBSCAN.calDisByID(     587493    ,     cutDis= 1500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(134.9751058,-2.9542447,2131.874,1836.318,0)
            useAV= True
            doDisDBSCAN.calDisByID(    307102    ,     cutDis= 2000 ,  SL=  5  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(138.3708474,-1.0247575,2151.254,2044.661,0)
            useAV= True
            doDisDBSCAN.calDisByID(     424646     ,     cutDis= 2000 ,  SL=  2 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(144.7590957,-0.8402892,2351.263,1988.457,0)
            useAV= True
            doDisDBSCAN.calDisByID(    175633    ,     cutDis= 3500 ,  SL= 4 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(124.2637665,1.3455286,1928.376,2136.718,0)
            useAV= True
            doDisDBSCAN.calDisByID(    83834    ,     cutDis= 3500 ,  SL=  4  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(120.8651537,4.1526339,2407.079,2093.112,0)
            useAV= True
            doDisDBSCAN.calDisByID(    47682    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =    box(114.2247541,2.8278665,2123.125,1483.833,0)
            useAV= True
            doDisDBSCAN.calDisByID(     285354     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(149.3279669,1.3576390,2197.768,1883.801,0)
            useAV= True
            doDisDBSCAN.calDisByID(    194343    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(148.5741263,-3.1096621,2797.068,2302.758,0)
            useAV= True
            doDisDBSCAN.calDisByID(     540504     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(122.5494172,1.7953462,2476.850,2372.194,0)
            useAV= True 
            doDisDBSCAN.calDisByID(    440382    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ###################
            lRange, bRange =     box(143.6496877,1.5701299,2412.892,2059.195,0)
            useAV= False
            doDisDBSCAN.calDisByID(    196496    ,     cutDis= 3500 ,  SL=  4 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(123.4232639,1.5150298,2662.439,2126.602,0)
            useAV= True
            doDisDBSCAN.calDisByID(     573038     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(135.1206968,-2.3162417,2721.048,2372.196,0)
            useAV= True
            doDisDBSCAN.calDisByID(    519984    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =      box(124.9395284,-0.2226757,2337.309,1883.801,0)
            useAV= True
            doDisDBSCAN.calDisByID(     144894     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(122.6651826,0.2190796,1982.369,2175.839,0)
            useAV= True
            doDisDBSCAN.calDisByID(     298228     ,     cutDis= 1500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(143.5396763,3.2171773,2156.772,1968.349,0)
            useAV= True 
            doDisDBSCAN.calDisByID(    539470    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =      box(122.8081109,2.8952438,1928.377,1753.951,0)
            useAV= True
            doDisDBSCAN.calDisByID(     56023     ,     cutDis= 2000 ,  SL= 4  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(107.6541434,-0.3932600,2243.817,2001.016,0)
            useAV= True
            doDisDBSCAN.calDisByID(    465617    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(122.0956273,0.9913132,1917.878,1788.673,0)
            useAV= True 
            doDisDBSCAN.calDisByID(     268775     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(113.9302464,-0.1219861,2197.768,1834.962,0)
            useAV= True 
            doDisDBSCAN.calDisByID(     174677     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(138.7190124,-0.3129850,2143.314,1793.385,0)
            useAV= True
            doDisDBSCAN.calDisByID(     179685     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(146.2529353,-1.6438863,1981.674,1860.544,0)
            useAV= True
            doDisDBSCAN.calDisByID(     470557     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(128.3704712,-3.8364852,1765.191,2483.827,0)
            useAV= True
            doDisDBSCAN.calDisByID(    331340     ,     cutDis= 1500 ,  SL=  2 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ################### AV GOOD 
            lRange, bRange =     box(148.5458397,-4.6290356,2572.016,3968.541,0)
            useAV= True
            doDisDBSCAN.calDisByID(     247523     ,     cutDis= 2000 ,  SL=  5  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask= False ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ###################
            lRange, bRange =     box(139.2615479,1.2915613,2218.699,2065.204,0)
            useAV= True
            doDisDBSCAN.calDisByID(     292151     ,     cutDis= 2000 ,  SL=  4  ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ################### AV GOOD
            lRange, bRange =     box(132.9202600,-4.5433709,2843.170,1749.993,0)
            useAV= True
            doDisDBSCAN.calDisByID(    275704     ,     cutDis= 2000 ,  SL= 3 ,  NL=  0.5   ,    lExpand=0.5, bExpand=0.5,  useMask=True  ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(115.5063999,-1.5336392,2240.467,2109.857,0)
            useAV= True
            doDisDBSCAN.calDisByID(    529414    ,     cutDis= 1500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=True ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ###################
            lRange, bRange =     box(108.5301165,-0.7704131,3074.363,3014.082,0)
            useAV= True
            doDisDBSCAN.calDisByID(     546878     ,     cutDis= 1000 ,  SL=  2,  NL= 0.1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(122.9319875,1.7922384,1792.712,1773.331,0)
            useAV= True 
            doDisDBSCAN.calDisByID(    436310    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(114.8623047,3.4396919,2575.691,2011.714,0)
            useAV= True 
            doDisDBSCAN.calDisByID(     47866     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange =     box(149.5340321,-1.4451328,2261.078,2113.031,0)
            useAV= True
            doDisDBSCAN.calDisByID(     192997     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ###################
            lRange, bRange = box(143.6346962,-1.1321016,2116.396,1975.079,0)
            useAV= True
            doDisDBSCAN.calDisByID(    592581    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(144.6943147,-4.1178905,2064.041,2087.298,0)
            useAV= True
            doDisDBSCAN.calDisByID(    279866    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(141.4123788,-1.2794374,1938.067,1843.855,0) 
            useAV= True
            doDisDBSCAN.calDisByID(    193952    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ###################
            lRange, bRange =     box(107.9277698,0.9643634,1791.390,1450.637,0)
            useAV= True 
            doDisDBSCAN.calDisByID(    83307    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################
            lRange, bRange =     box(142.1593824,3.3392803,2325.678,1986.517,0)
            useAV= True 
            doDisDBSCAN.calDisByID(     181388     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


 
            ################### AV GOOD
            lRange, bRange =    box(123.0272224,-1.2390893,1767.518,1651.234,0)
            useAV= True
            doDisDBSCAN.calDisByID(     339915     ,     cutDis= 2000 ,  SL=  1  ,  NL=  0.1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS= True  ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ################### 
            lRange, bRange =     box(118.8986377,0.0664718,2141.952,1911.709,0)
            useAV= True
            doDisDBSCAN.calDisByID(     117329     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ################### 
            lRange, bRange =     box(109.5388763,3.9016272,1780.599,1643.319,0)
            useAV= True
            doDisDBSCAN.calDisByID(     555096     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ################### 
            lRange, bRange =     box(125.5799786,1.0115412,2296.609,1957.448,0)
            useAV= True
            doDisDBSCAN.calDisByID(     81260     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ################### 
            lRange, bRange =     box(114.4536331,-0.9006720,2345.867,2018.819,0)
            useAV= True 
            doDisDBSCAN.calDisByID(    126117    ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ################### 
            lRange, bRange =     box(127.2601397,0.8859695,2151.254,1996.209,0)
            useAV= True 
            doDisDBSCAN.calDisByID(     86717     ,     cutDis= 3500 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ################### 
            lRange, bRange =     box(141.8530994,3.5422883,1986.519,1836.318,0)
            useAV= True
            doDisDBSCAN.calDisByID(     341117     ,     cutDis= 1500 ,  SL= 4  ,  NL= 0.5  ,    lExpand=0.5, bExpand=0.5,  useMask=True ,  useForegroundFITS=True  ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################  
            lRange, bRange =     box(111.1766954,-1.3435776,2714.069,2309.401,0)
            useAV= True
            doDisDBSCAN.calDisByID(     214215     ,     cutDis= 3500 ,  SL=  4 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask= False  ,  useForegroundFITS= False  ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )



            ###################  
            lRange, bRange =     box(113.3251699,1.6342449,2139.949,1631.206,0)
            useAV= True
            doDisDBSCAN.calDisByID(     517730     ,     cutDis= 2000 ,  SL=  2 ,  NL= 0.5  ,    lExpand=0.5, bExpand=0.5,  useMask= True ,  useForegroundFITS= False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ###################  
            lRange, bRange =     box(148.4151191,5.0433696,1749.842,1691.235,0)
            useAV= True
            doDisDBSCAN.calDisByID(     487657     ,     cutDis= 3500 ,  SL=  2,  NL=  0.5  ,    lExpand=0.5, bExpand=0.5,  useMask= True ,  useForegroundFITS= False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ###################
            lRange, bRange =     box(120.2537026,-4.2901456,2203.878,1564.585,0)
            useAV= True
            doDisDBSCAN.calDisByID(     277351     ,     cutDis= 2000 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask= True ,  useForegroundFITS= False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ###################
            lRange, bRange =     box(120.2537026,-4.2901456,2203.878,1564.585,0)
            useAV= True
            doDisDBSCAN.calDisByID(     277351     ,     cutDis= 2000 ,  SL=  3 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask= True ,  useForegroundFITS= False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ###################  
            lRange, bRange =     box(121.0553613,1.9242837,2523.363,2139.626,0)
            useAV= True
            doDisDBSCAN.calDisByID(    539039    ,     cutDis= 3500 ,  SL=  2 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask= True ,  useForegroundFITS= False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )
            ############ 66+

            ###################   AV GOOD
            lRange, bRange =     box(111.3605110,-0.5373801,1942.105,1939.722,0)
            useAV= True
            doDisDBSCAN.calDisByID(     349203     ,     cutDis= 2500 ,  SL=  4 ,  NL= 0.5  ,    lExpand=0.5, bExpand=0.5,  useMask= True ,  useForegroundFITS= False   ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ################### AV GOOD
            lRange, bRange =  box(111.7227096,0.0367231,2093.112,1984.271,0)
            useAV= True
            doDisDBSCAN.calDisByID(      176665    , lowerDisCut = 800,    cutDis= 3000 ,  SL= 3 ,  NL= 0.5  ,    lExpand=0.5, bExpand=0.5,  useMask= True ,  useForegroundFITS= True  ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################   AV GOOD
            lRange, bRange =  box(143.6427735,-3.2409614,6636.271,4487.834,0)
            useAV= True
            doDisDBSCAN.calDisByID(      186740    ,   lowerDisCut = 1200,   cutDis= 3000 ,  SL= 3.5 ,  NL= 0.5  ,    lExpand=0.5, bExpand=0.5,  useMask= False  ,  useForegroundFITS= True  ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ###################   AV GOOD
            lRange, bRange =   box(150.0671411,-1.4291097,1283.970,2102.803,0)
            useAV= True
            doDisDBSCAN.calDisByID(       301901   ,     cutDis= 2000 ,  SL= 5 ,  NL= 1  ,    lExpand=0.5, bExpand=0.5,  useMask= True ,  useForegroundFITS= False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ###################   AV GOOD
            lRange, bRange =    box(116.4657572,-4.5258504,1947.755,1682.531,0)
            useAV= True

            doDisDBSCAN.calDisByID(       459778     ,     cutDis= 2000 ,  SL= 3 ,  NL= 0.5  ,    lExpand=0.5, bExpand=0.5,  useMask= True  ,  useForegroundFITS= False   ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ################### AV GOOD
            lRange, bRange =     box(123.6329542,1.0038744,2246.607,1625.651,0)
            useAV= True
            doDisDBSCAN.calDisByID(       365278  ,     cutDis= 2000 ,  SL= 3  ,  NL= 0.5  ,    lExpand=0.5, bExpand=0.5,  useMask= True   ,  useForegroundFITS= True  ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ################### AV GOOD
            lRange, bRange =   box(136.8672987,-2.7192093,2020.435,1889.615,0)
            useAV= True
            doDisDBSCAN.calDisByID(        413003 ,     cutDis= 2000 ,  SL= 2 ,  NL= 0.5  ,    lExpand=0.5, bExpand=0.5,  useMask= True ,  useForegroundFITS= False   ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

            ################### barely no detection
            lRange, bRange =   box(147.9672285,-5.0126169,2080.506,1633.692,0)
            useAV= True
            doDisDBSCAN.calDisByID(      269315    ,     cutDis= 2000 ,  SL=   3   ,  NL= 0.5  ,    lExpand=0.5, bExpand=0.5,  useMask= True   ,  useForegroundFITS= False   ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )


            ###################
            lRange, bRange =   box(136.7033449,-2.8819219,2176.960,2001.995,0)
            useAV= False
            doDisDBSCAN.calDisByID(        276696  ,     cutDis= 1500 ,  SL= 2.5 ,  NL= 0.5  ,    lExpand=0.5, bExpand=0.5,  useMask= True ,  useForegroundFITS= False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )

        else: #test distances      useForegroundFITS=True,  lowerDisCut  True

            ###################
            lRange, bRange =   box(126.8728955,0.7492943,2072.655,1588.138,0)
            useAV= True
            doDisDBSCAN.calDisByID(      115715    ,     cutDis= 3500 ,  SL= 3 ,  NL= 0.5  ,    lExpand=0.5, bExpand=0.5,  useMask= True ,  useForegroundFITS= False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )





            return
            #################### intital
            lRange, bRange = XXX
            useAV= True
            doDisDBSCAN.calDisByID(     XXX   ,     cutDis= 3500 ,  SL= 3 ,  NL= 0.5  ,    lExpand=0.5, bExpand=0.5,  useMask= True ,  useForegroundFITS= False ,    useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=  False       )



    def runGoodSources(self, useAV=True):
        """
        A function that deal with distance
        :return:
        """
        print "Tesing distance"

        doDisDBSCAN = dbscanDis("disQ2", useAV)

        print doDisDBSCAN.maskPath
        print doDisDBSCAN.CO12FITS



        return


        #################### AG GOOD
        lRange, bRange = box(146.2001185, -4.3333877, 7052.951, 6510.417,  0)  # region2  749 pc, do they belong the same cloud?
        useAV = True
        doDisDBSCAN.calDisByID(325105, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True, cutDis= 800 ,
                               useBaseline=True,  SL=4 ,  lExpand=0.5, bExpand=0.5, NL=1 )




        #################### AV GOOD
        lRange, bRange =  box(128.3145165,0.9473845,1974.406,1829.051,0)
        useAV = True
        doDisDBSCAN.calDisByID(296887, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1500,
                               useBaseline=True, SL=3 , lExpand=0.5, bExpand=0.5, NL=1)


        #################### AG GOOD
        lRange, bRange = box(122.8705620, -2.2955322, 2704.301, 2260.561, 0)
        useAV = False
        doDisDBSCAN.calDisByID(270975, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=2000,
                               useBaseline=True, SL=3, lExpand=0.5, bExpand=0.5, NL=1)



        #################### GOOD
        lRange, bRange = box(115.5554637, -3.0557037, 13644.000, 11628.000,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(402638, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=800, useBaseline=True, SL= 6 ,
                               lExpand=0.5, bExpand=0.5, NL=2)



        #################### AG GOOD
        lRange, bRange = box(141.3882957, -3.5589113, 5791.667, 6208.333,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(252392, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis= 1500,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)



        #################### AV GOOD
        lRange, bRange = box(135.7161406, 3.4395674, 3906.250, 3420.139, 0)  #
        useAV = True
        doDisDBSCAN.calDisByID(316030, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1000,
                               useBaseline=True, SL= 4   , lExpand=0.5, bExpand=0.5, NL=1 )






        #################### AV GOOD
        lRange, bRange = box(131.2217699, -0.5849484, 2372.194, 1837.288, 0)
        useAV = True
        doDisDBSCAN.calDisByID(280221, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1500,
                               useBaseline=True, SL= 4 , lExpand=0.5, bExpand=0.5, NL=1)

        

        ####################
        lRange, bRange = box(145.1316260,3.4706396,24552.000,12672.000,0) # box(145.2605568, 2.3391936, 23483.520, 20373.120, 0)
        useAV = False
        doDisDBSCAN.calDisByID(360300, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=500, useBaseline=True, SL= 6 ,
                               lExpand=0.5, bExpand=0.5, NL=1)
        
        

        #################### AG GOOD, do not use mask
        lRange, bRange = box(142.1204632, -3.5906363, 10517.940, 5685.764,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(248354, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True, useMask=False,
                               cutDis=1500, useBaseline=True, SL=4 , lExpand=0.5, bExpand=0.5, NL=2)



        ####################  AV GOOD
        lRange, bRange = box(134.8962161, 1.6825291, 1821.008, 3510.033, 0)
        useAV = True
        doDisDBSCAN.calDisByID(148021, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=3500,
                               useBaseline=True, SL= 2 , lExpand=0.5, bExpand=0.5, NL=1)



        #################### AG GOOD
        lRange, bRange = box(131.6478138, -1.5676965, 4042.219, 3305.443,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = True
        doDisDBSCAN.calDisByID(267103, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True, cutDis=1500,
                               useBaseline=True, SL= 4 , lExpand=0.5, bExpand=0.5, NL=1)



        ####################
        lRange, bRange = box(122.5996191, -2.3075627, 2476.850, 2267.538, 0)
        useAV = True
        doDisDBSCAN.calDisByID(250260, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=2000,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)


        #################### AG GOOD
        lRange, bRange = box(144.2079849, -3.8978160, 7042.904, 5505.723,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(339429, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True, cutDis=1000,
                               useBaseline=True, SL=4, lExpand=0.5, bExpand=0.5, NL=2)



        
        ####################   AV GOOD
        lRange, bRange = box(121.7015519, -2.1891454, 2779.653, 2176.837, 0)
        useAV = True
        doDisDBSCAN.calDisByID(252641, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=2000,
                               useBaseline=True, SL=4, lExpand=0.5, bExpand=0.5, NL=2)

        
        
        #################### AV GOOD
        lRange, bRange = box(109.1132135, -0.1513716, 3122.589, 2109.857, 0)
        useAV = True
        doDisDBSCAN.calDisByID(465607, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1500,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=1)

        
        #################### AG GOOD
        lRange, bRange = box(130.3226687, 0.6311805, 5058.836, 3284.144,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(293789, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1500,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)

        ####################  AV GOOD
        lRange, bRange = box(133.2848100, -2.6336915, 4448.785, 3592.785,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = True
        doDisDBSCAN.calDisByID(277600, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1500,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)


        #################### AG GOOD
        lRange, bRange = box(144.2706059, 4.1616455, 9461.806, 6822.917,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(281664, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1500, useBaseline=True, SL=5,
                               lExpand=0.5, bExpand=0.5, NL=2)


        #################### AG GOOD
        lRange, bRange = box(146.2876565, -2.3809948, 6293.396, 3098.473,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(275042, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=2000,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)


        ###################
        lRange, bRange = box(119.5661502, 2.8275915, 3183.969, 3097.806, 0)
        useAV = False
        doDisDBSCAN.calDisByID(238938, cutDis=2000, SL=5, NL=1, lExpand=0.5, bExpand=0.5, useForegroundFITS=True,
                               useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False)

        ################### AG GOOD
        lRange, bRange = box(129.7549294, -0.0526735, 4591.451, 5405.253, 0)
        useAV = False
        doDisDBSCAN.calDisByID(523676, cutDis=1500, SL=5, NL=2, lExpand=0.5, bExpand=0.5, useForegroundFITS=False,
                               useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False)

        ################### AV GOOD
        lRange, bRange = box(131.5835131, -0.1102546, 2413.194, 3784.722, 0)
        useAV = True
        doDisDBSCAN.calDisByID(267299, cutDis=1500, SL=4, NL=2, lExpand=0.5, bExpand=0.5, useForegroundFITS=True,
                               useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False)


        ################### AG GOOD
        lRange, bRange = box(144.5320592, -2.5559178, 5786.480, 4921.203, 0)
        useAV = False
        doDisDBSCAN.calDisByID(247603, cutDis=1100, SL=4, NL=1, lExpand=0.5, bExpand=0.5, useForegroundFITS=True,
                               useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False)

        ################### AG GOOD
        lRange, bRange = box(143.5637502, -1.7532601, 4834.587, 6667.149, 0)
        useAV = False
        doDisDBSCAN.calDisByID(258072, cutDis=1500, SL=5, NL=1, lExpand=0.5, bExpand=0.5, useForegroundFITS=True,
                               useBaseline=True, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False)

        #################### AV GOOD
        lRange, bRange = box(143.1263942, 0.7612341, 11611.111, 12277.199, 0)
        useAV = True
        doDisDBSCAN.calDisByID(454540, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=500,
                               useBaseline=True, SL=3, lExpand=0.5, bExpand=0.5, NL=1)



        #################### AV GOOD
        lRange, bRange = box(124.1429992, 0.5965697, 2447.779, 1755.889, 0)
        useAV = True
        doDisDBSCAN.calDisByID(293831, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=2000,
                               useBaseline=True, SL=3, lExpand=0.5, bExpand=0.5, NL=2)



        ####################AG GOOD
        lRange, bRange = box(140.1897786, -3.3067479, 12208.333, 11166.667, 0)
        useAV = False
        doDisDBSCAN.calDisByID(259585, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1000, useBaseline=True, SL=10,
                               lExpand=0.5, bExpand=0.5, NL=1)



        #################### AV GOOD
        lRange, bRange = box(112.1463454, -1.5355212, 3194.927, 3448.110, 0)
        useAV = True
        doDisDBSCAN.calDisByID(353634, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, useMask=False,
                               cutDis=1500, useBaseline=True, SL=3, lExpand=0.5, bExpand=0.5, NL=1)





        #################### Good
        lRange, bRange = box(136.3514000, -1.8665333, 8418.240, 6432.480,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(260426, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1500, useBaseline=True, SL=5,
                               lExpand=0.5, bExpand=0.5, NL=2)

        #################### Good
        lRange, bRange = box(148.8437794, 3.2313697, 9994.338, 9470.260,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(539067, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=500, useBaseline=True, SL=5,
                               lExpand=0.5, bExpand=0.5, NL=2)

        #################### GooD
        lRange, bRange = box(112.3148208, -2.5667783, 17020.833, 12875.000,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(412961, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1000, useBaseline=True, SL=5,
                               lExpand=0.5, bExpand=0.5, NL=2)

        #################### GOOD
        lRange, bRange = box(106.5445307, 4.1271120, 12378.472, 7760.417,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(264174, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1500, useBaseline=True, SL=5,
                               lExpand=0.5, bExpand=0.5, NL=2)



        #################### GOOD
        lRange, bRange = box(122.4631200, -0.3014400, 12312.000, 13478.400,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(234734, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1500, useBaseline=True, SL=5,
                               lExpand=0.5, bExpand=0.5, NL=2)



        #################### AG GOOD
        lRange, bRange = box(105.7726078, 0.7864653, 7173.515, 7685.909,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(403058, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1000, useBaseline=True, SL=5,
                               lExpand=0.5, bExpand=0.5, NL=2)

        #################### # AG  GooD
        lRange, bRange = box(148.1421364, 0.0953276, 7263.937, 6761.590,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(295121, useAV=useAV, lRange=lRange, bRange=bRange, cutDis=1500, useHighAv=False,
                               useBaseline=True, SL=7, lExpand=0.5, bExpand=0.5, NL=2)

        #################### AV GOOD
        lRange, bRange = box(147.8746820, -4.4293344, 7414.641, 5401.235,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(215739, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True, cutDis=1500,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)







        #################### AG GOOD
        lRange, bRange = box(143.3018896, 2.9311760, 5895.544, 4605.517,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(285283, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=2000,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)




        #################### AG GOOD
        lRange, bRange = box(129.9984357, -2.2159297, 7660.791, 3642.015,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(263681, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True, cutDis=2000,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)

        #################### AG GOOD
        lRange, bRange = box(115.9509875, -0.5314069, 5309.606, 3964.120,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(371984, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=800,
                               useBaseline=True, SL=3, lExpand=0.5, bExpand=0.5, NL=2)



        #################### AG GOOD
        lRange, bRange = box(106.8876165, 0.9620629, 6510.417, 2742.814,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(277035, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1000,
                               useBaseline=True, SL=4, lExpand=0.5, bExpand=0.5, NL=2)


        #################### AG GOOD
        lRange, bRange = box(130.9951651, -0.9276344, 5987.976, 3154.739,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(312594, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1500,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)

        #################### AG GOOD
        lRange, bRange = box(149.5081494, -1.1220855, 3566.663, 4738.806,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(296464, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=2000,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)

        #################### AG GOOD
        lRange, bRange = box(120.9062689, -1.4743259, 4895.557, 4005.984,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(237970, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=2000,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)

        #################### AG GOOD
        lRange, bRange = box(121.9820486, -0.4433275, 5280.671, 5613.426,
                             0)  # region2  749 pc, do they belong the same cloud?
        useAV = False
        doDisDBSCAN.calDisByID(516881, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1000,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)


        ####################
        lRange, bRange = box(143.8073971, -3.3211090, 5714.699, 4533.179, 0)
        useAV = False
        doDisDBSCAN.calDisByID(186740, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, useMask=False,
                               cutDis=2500, useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=1)

        ####################  AG GOOD
        lRange, bRange = box(125.6333998, -2.5121796, 5304.784, 3858.025, 0)  #
        useAV = False
        doDisDBSCAN.calDisByID(236543, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=2500,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)



        #################### AG GOOD
        lRange, bRange = box(106.4580168, 1.6852660, 5353.009, 2387.153, 0)
        useAV = False
        doDisDBSCAN.calDisByID(264979, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1500,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)

        #################### AV GOOD
        lRange, bRange = box(129.3505866, -1.6089215, 3677.180, 3345.631, 0)
        useAV = True
        doDisDBSCAN.calDisByID(303216, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=2000,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)

        #################### AV GOOD
        lRange, bRange = box(134.0382861, -4.0111832, 3701.292, 3339.603, 0)
        useAV = True
        doDisDBSCAN.calDisByID(270697 , useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1500,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)

        #################### AG GOOD
        lRange, bRange = box(108.3815358, 0.2959505, 3978.588, 3327.546, 0)
        useAV = False
        doDisDBSCAN.calDisByID(307747, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1500,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)

        #################### AG GOOD
        lRange, bRange = box(105.8580949, -1.1969070, 4310.137, 4430.700, 0)
        useAV = False
        doDisDBSCAN.calDisByID(331873, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1500,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)

        #################### AG GOOD
        lRange, bRange = box(116.9641660, 3.6574644, 3935.051, 3148.041, 0)
        useAV = False
        doDisDBSCAN.calDisByID(309918, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1000,
                               useBaseline=True, SL=4, lExpand=0.5, bExpand=0.5, NL=2)

        ####################  AG GOOD
        lRange, bRange = box(114.6841598, -0.1312470, 4064.368, 3277.906, 0)
        useAV = False
        doDisDBSCAN.calDisByID(178245, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=3000,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=1)

        #################### AG GOOD
        lRange, bRange = box(111.0382748, 3.5244203, 2235.444, 3951.796, 0)
        useAV = False
        doDisDBSCAN.calDisByID(237244, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1500,
                               useBaseline=True, SL=4, lExpand=0.5, bExpand=0.5, NL=2)

        #################### AG GOOD
        lRange, bRange = box(107.7561084, 2.9155225, 2883.137, 2699.166, 0)
        useAV = False
        doDisDBSCAN.calDisByID(282690, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=True, cutDis=1500,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)

        #################### AG GOOD
        lRange, bRange = box(124.6732187, -1.5100658, 3144.514, 2572.784, 0)
        useAV = False
        doDisDBSCAN.calDisByID(251759, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=2000,
                               useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=2)

        #################### AG GOOD
        lRange, bRange = box(121.7725262, -0.4722424, 3300.141, 3021.059, 0)
        useAV = False
        doDisDBSCAN.calDisByID(271145, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1500,
                               useBaseline=True, SL=4, lExpand=0.5, bExpand=0.5, NL=2)

        ####################  AV GOOD
        lRange, bRange = box(111.0855673, -0.1725342, 3017.618, 3437.902, 0)
        useAV = True
        doDisDBSCAN.calDisByID(309300, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, useMask=True,
                               cutDis=2000, useBaseline=True, SL=4, lExpand=0.5, bExpand=0.5, NL=2)




        #################### AV GDOOD
        lRange, bRange = box(132.0123736, -2.5034317, 4340.279, 3269.676, 0)
        useAV = True
        doDisDBSCAN.calDisByID(262071, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, useMask=False,
                               cutDis=1500, useBaseline=True, SL=5, lExpand=0.5, bExpand=0.5, NL=1)





        #################### AV GOOD
        lRange, bRange = box(121.4695202, 0.3812644, 3211.806, 2835.649, 0)
        useAV = True
        doDisDBSCAN.calDisByID(418923, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=1500,
                               useBaseline=True, SL=3, lExpand=0.5, bExpand=0.5, NL=1)

        ####################AV GOOD
        lRange, bRange = box(124.6312621, -1.9743158, 2446.810, 2437.119, 0)
        useAV = True
        doDisDBSCAN.calDisByID(305787, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=2000,
                               useBaseline=True, SL=4, lExpand=0.5, bExpand=0.5, NL=2)

        #################### AV GOOD
        lRange, bRange = box(130.0402813, -4.4565050, 2511.735, 2609.413, 0)
        useAV = True
        doDisDBSCAN.calDisByID(149173, useAV=useAV, lRange=lRange, bRange=bRange, useHighAv=False, cutDis=3000,
                               useBaseline=True, SL=3, lExpand=0.5, bExpand=0.5, NL=2)





    def ZZZ(self):
        """
        pass
        :return:
        """