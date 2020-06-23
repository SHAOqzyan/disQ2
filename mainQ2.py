
from distanceQ2 import disQ2
doQ2=disQ2()

if 0: #crop fits
    doQ2.prepareData()


if 0:
    pass
    doQ2.processDBSCAN( doQ2.cropRawCO12FITSLocal )

if 1:

    doQ2.mergeByVaxis(doQ2.cropRawCO12FITSPer ,doQ2.cropRawCO12FITSLocal )