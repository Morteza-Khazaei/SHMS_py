#! /usr/bin/python
import gdal, osr
import numpy as np
import numexpr as ne
import matplotlib.pyplot as plt
from scipy.stats import linregress
import gdalnumeric



class DISPATCH():


    def __init__(self, hdfFile, MOD11A1, h5File, ndv=None):

        # Product subdataset information
        gq_lst = 0  # 1km LST
        gq_red = 1  # 250m red band
        gq_nir = 2  # 250m NIR band

        # Open and find subdatasets
        gq_ds = gdal.Open(hdfFile, gdal.GA_ReadOnly)      # MODIS
        mod11_ds = gdal.Open(MOD11A1, gdal.GA_ReadOnly)   # MODIS
        self.sm_ds = gdal.Open(h5File, gdal.GA_ReadOnly)  # SMAP

        # Read in datasets
        self.array_sm = gdalnumeric.LoadFile(h5File)
        ds_red = gdal.Open(gq_ds.GetSubDatasets()[gq_red][0], gdal.GA_ReadOnly)     # MODIS Dataset
        ds_nir = gdal.Open(gq_ds.GetSubDatasets()[gq_nir][0], gdal.GA_ReadOnly)     # MODIS Dataset
        ds_lst = gdal.Open(mod11_ds.GetSubDatasets()[gq_lst][0], gdal.GA_ReadOnly)  # MODIS Dataset

        # Wrap To WGS84
        geotref = ds_red.GetProjectionRef()
        geoTransform = ds_red.GetGeoTransform()
        # pixel = geoTransform[1] / 100000
        pixel = 0.001853250866278332

        # MOD11A1 Geotransform
        geotref_lst = ds_lst.GetProjectionRef()
        geoTransform_lst = ds_lst.GetGeoTransform()
        pixel_lst = geoTransform_lst[1] / 100000

        self.w_ds_red = self.reproject_dataset(dataset=ds_red, pixel_spacing=pixel, epsg_from=geotref, epsg_to=4326)
        self.w_ds_nir = self.reproject_dataset(dataset=ds_nir, pixel_spacing=pixel, epsg_from=geotref, epsg_to=4326)
        self.w_ds_lst = self.reproject_dataset(dataset=ds_lst, pixel_spacing=pixel_lst, epsg_from=geotref_lst, epsg_to=4326)

        # Write out image
        red = self.w_ds_red.GetRasterBand(1).ReadAsArray().astype(float)
        nir = self.w_ds_nir.GetRasterBand(1).ReadAsArray().astype(float)
        lst = self.w_ds_lst.GetRasterBand(1).ReadAsArray().astype(float)

        # masking invalid values
        eqn = '(red <= 0) | (red >= 10000) | (nir <= 0) | (nir >= 10000)'

        invalid = ne.evaluate(eqn)

        # Apply valid range mask to data
        red[invalid] = ndv
        nir[invalid] = ndv

        # Set Variabels
        self.red = red
        self.nir = nir
        self.lst = ne.evaluate('(lst / 50.0)').astype(float)
        self.sm = self.Cilpper(srcArray=self.array_sm, smapDataset=self.sm_ds, modisWrap=self.w_ds_lst)
        self.clippedSM = self.smapClipRsizer(modisWrap=self.w_ds_lst, clipArray=self.smapClipped,
                                             transform=self.clippedGeotranse)
        self.ndvi = self.NDVI()
        self.fvc = self.FVC()
        self.fv1km = self.resmapleFvc(fvc250m=self.fvc, lst1km=self.lst)
        self.Tv1km = self.calculateTv1km(LST=self.lst, FVC=self.fv1km)
        self.Ts1km = self.calculateTs1km(Tmodis=self.lst, fv1km=self.fv1km, Tv1km=self.Tv1km)
        self.see = self.SEEMODIS1km(Ts1km=self.Ts1km) #, Tsmax=self.Tsmax, Tsmin=self.Tsmin)
        self.SEE40km = self.calculateMeanSEE(SEE=self.see, clip=self.sm)
        self.SEEavg = self.SEEMean(modisWrap=self.w_ds_lst, seeMean=self.SEE40km, transform=self.clippedGeotranse)
        self.smp = self.soilMoistureParameter(rescaleSMAP=self.clippedSM, rescaleMeanSEE=self.SEEavg)
        self.rescalclippedSM = self.smapClipRsizer(modisWrap=self.w_ds_lst, clipArray=self.smapClipped,
                                                   transform=self.clippedGeotranse)
        self.disaggregated = self.Disaggregation(rescaleSMAP=self.rescalclippedSM,
                                                 rescaleMeanSEE=self.SEEavg,
                                                 SEEmodis1km=self.see, SMp=self.smp)







    def reproject_dataset(self,dataset, pixel_spacing=None, epsg_from=None, epsg_to=None):
        """
        A sample function to reproject and resample a GDAL dataset from within
        Python. The idea here is to reproject from one system to another, as well
        as to change the pixel size. The procedure is slightly long-winded, but
        goes like this:

        1. Set up the two Spatial Reference systems.
        2. Open the original dataset, and get the geotransform
        3. Calculate bounds of new geotransform by projecting the UL corners
        4. Calculate the number of pixels with the new projection & spacing
        5. Create an in-memory raster dataset
        6. Perform the projection
        """
        # Define the UK OSNG, see <http://spatialreference.org/ref/epsg/27700/>
        out = osr.SpatialReference()
        out.ImportFromEPSG(epsg_to)
        sinu = osr.SpatialReference()
        # wgs84.ImportFromEPSG ( epsg_from )
        sinu.ImportFromWkt(epsg_from)
        tx = osr.CoordinateTransformation(sinu, out)
        # Up to here, all  the projection have been defined, as well as a
        # transformation from the from to the  to :)
        # We now open the dataset
        # g = gdal.Open ( dataset )

        # Get the Geotransform vector
        geo_t = dataset.GetGeoTransform()
        x_size = dataset.RasterXSize  # Raster xsize
        y_size = dataset.RasterYSize  # Raster ysize

        # Work out the boundaries of the new dataset in the target projection
        (ulx, uly, ulz) = tx.TransformPoint(geo_t[0], geo_t[3])
        (lrx, lry, lrz) = tx.TransformPoint(geo_t[0] + geo_t[1] * x_size, geo_t[3] + geo_t[5] * y_size)
        # print ulx, uly, ulz
        # print lrx, lry, lrz
        # See how using 27700 and WGS84 introduces a z-value!
        # Now, we create an in-memory raster
        mem_drv = gdal.GetDriverByName('MEM')
        # The size of the raster is given the new projection and pixel spacing
        # Using the values we calculated above. Also, setting it to store one band
        # and to use Float32 data type.
        dest = mem_drv.Create('', int((lrx - ulx) / pixel_spacing), int((uly - lry) / pixel_spacing), 1,
                              gdal.GDT_Float32)
        # Calculate the new geotransform
        new_geo = (ulx, pixel_spacing, geo_t[2], uly, geo_t[4], -pixel_spacing)
        # Set the geotransform
        dest.SetGeoTransform(new_geo)
        dest.SetProjection(out.ExportToWkt())
        # Perform the projection/resampling
        res = gdal.ReprojectImage(dataset, dest, sinu.ExportToWkt(), out.ExportToWkt(), gdal.GRA_NearestNeighbour)

        return dest



    def NDVI(self, type=float, ndv=None):

        # Calculate NDVI
        red = self.red
        nir = self.nir
        ndvi = ne.evaluate('(nir - red) / (nir + red)').astype(type)

        # masking invalid values
        eqn = '(ndvi <= -1) | (ndvi >= 1)'
        invalid = ne.evaluate(eqn)

        # Apply valid range mask to data
        ndvi[invalid] = ndv


        print "NDVI is ok!"
        return ndvi



    def FVC(self, type=float, ndv=0):

        # NDVI Variables (Denote pure vegetation (set to 0.9) and vegetation-free soil (set to 0.15))
        NDVIv = 0.9
        NDVIs = 0.15
        NDVImodis = self.ndvi

        # Calculate FVC (represents the fraction of vegetation within each pixel)
        fv1km = ne.evaluate('( NDVImodis - NDVIs )/( NDVIv - NDVIs )').astype(type)

        # masking invalid values
        eqn = '(fv1km <= -1) | (fv1km >= 1)'
        invalid = ne.evaluate(eqn)

        # Apply valid range mask to data
        fv1km[invalid] = ndv


        print "FVC is ok!"
        return fv1km



    def SEEMODIS1km(self, Ts1km=None, ndv=0, type=float):

        # Set Variables
        Tsmax = self.Tsmax
        Tsmin = self.Tsmin


        # Calculate SEE-MODIS 1KM
        see = ne.evaluate('((Tsmax - Ts1km) / (Tsmax - Tsmin))').astype(type)

        # # masking invalid values
        # eqn = '(see < 0) | (see > 1)'
        # invalid = ne.evaluate(eqn)
        #
        # # Apply valid range mask to data
        # see[invalid] = ndv


        print "SEE-MODIS-1KM is ok!"
        return see



    def calculateTs1km(self, Tmodis=None, fv1km=None, Tv1km=None, ndv=None, type=float):

        # Set Variables
        # Tmodis = Tmodis / 50.0

        # Calculate SEE-MODIS 1KM
        Ts1km = ne.evaluate('(Tmodis - fv1km * Tv1km) / (1 - fv1km)').astype(type)

        # masking invalid values
        eqn = '(Ts1km < -500) | (Ts1km > 500)'
        invalid = ne.evaluate(eqn)

        # Apply valid range mask to data
        Ts1km[invalid] = ndv

        print "Ts-1KM is ok!"
        return Ts1km



    def Disaggregation(self, rescaleSMAP=None, rescaleMeanSEE=None, SEEmodis1km=None, SMp=None, type=float):

        """  where SMH is the disaggregated soil moisture at 250 m resolution, SM is the SMAP soil moisture, NSMI
        is the variable derived from the MODIS 250 m resolution surface reflectance, NSMImean is the mean
        value at 36 km resolution scale and SM / NSMI (the variance conversion factor) is the partial derivative
        computed statistically. For clarity, bold font represents variables derived at SMAP-equal resolution. """

        D_SMp_SEE = (SMp / np.pi) * np.sqrt(-SEEmodis1km**2 + SEEmodis1km)

        SMh = ne.evaluate('rescaleSMAP + D_SMp_SEE * (SEEmodis1km - rescaleMeanSEE)').astype(type)


        return SMh



    def Cilpper(self, smapDataset=None, srcArray=None, modisWrap=None, ndv=None):

        # Set Variables
        # smapDataset = None
        # modisWrap = None
        # srcArray = smapDataset.GetRasterBand(1).ReadAsArray()

        def world2Pixel(geoMatrix, x, y):
            """
            Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
            the pixel location of a geospatial coordinate
            """
            ulX = geoMatrix[0]
            ulY = geoMatrix[3]
            xDist = geoMatrix[1]
            yDist = geoMatrix[5]
            rtnX = geoMatrix[2]
            rtnY = geoMatrix[4]
            pixel = int((x - ulX) / xDist)
            line = int((ulY - y) / xDist)
            return (pixel, line)

        def Corrector(geoMatrix, pixel, line):
            ulX = geoMatrix[0]
            ulY = geoMatrix[3]
            xDist = geoMatrix[1]
            yDist = geoMatrix[5]
            x = xDist * pixel + ulX
            y = yDist * line + ulY
            return (x, y)

        # Start Clip Raster
        ulXx, xres, xskew, ulYy, yskew, yres = modisWrap.GetGeoTransform()
        lrXx = ulXx + (modisWrap.RasterXSize * xres)
        lrYy = ulYy + (modisWrap.RasterYSize * yres)
        # print ulXx, lrXx
        # print ulYy, lrYy

        geoTrans = smapDataset.GetGeoTransform()
        print geoTrans
        ulX, ulY = world2Pixel(geoTrans, ulXx, ulYy)
        lrX, lrY = world2Pixel(geoTrans, lrXx, lrYy)
        # print ulX, lrX
        # print ulY, lrY

        # Clip section : in this section clipper select ulY:lrY and ulX:lrX
        clip = srcArray[ulY:lrY, ulX:lrX]
        # masking invalid values
        eqn = '(clip <= -1) | (clip >= 1)'
        invalid = ne.evaluate(eqn)

        # Apply valid range mask to data
        clip[invalid] = ndv
        print "Smap Soil Moisture has {} Column and {} Row That They Were Clipped {} Column and {} Row From It ".format(srcArray.shape[1], srcArray.shape[0], clip.shape[1], clip.shape[0])

        # Create a new geomatry for the image
        correct = Corrector(geoTrans, ulX, ulY)
        geoTransWrap = list(modisWrap.GetGeoTransform())
        geoTransWrap[0] = correct[0]  # -115.7676365011655
        geoTransWrap[3] = correct[1]  # 30.14823810492634
        geoTransWrap[1] = geoTrans[1]
        geoTransWrap[5] = geoTrans[5]

        # # Create a new geomatrix for the image
        # geoTransWrap = list(modisWrap.GetGeoTransform())
        # geoTransWrap[1] = geoTrans[1]
        # geoTransWrap[5] = geoTrans[5]
        # print geoTransWrap
        self.smapClipped = clip
        self.clippedGeotranse = np.asarray(geoTransWrap)

        return clip



    def smapClipRsizer(self, modisWrap=None, clipArray=None, transform=None):
        '''

        :param modisWrap:   get MODIS wrapped array
        :param clipArray:   get SMAP clipped array
        :param transform:   get SMAP corrected geotransform
        :return:            return rescaled and resampled array from SMAP clipped using MODIS wrap
        '''


        # give geotransform from clipped smap raster
        xOrigin = transform[0]
        yOrigin = transform[3]
        pixelWidth = transform[1]
        pixelHeight = -transform[5]

        # get MODIS wrapped transform for calculate center of each pixel in clipped raster
        geoTrans = list(modisWrap.GetGeoTransform())
        lrXx_m = geoTrans[0] + (modisWrap.RasterXSize * geoTrans[1])
        lrYy_m = geoTrans[3] + (modisWrap.RasterYSize * geoTrans[5])

        # Calculate X and Y offset
        ulXx, xres, xskew, ulYy, yskew, yres = self.clippedGeotranse
        lrXx_s = ulXx + (self.smapClipped.shape[1] * xres)
        lrYy_s = ulYy + (self.smapClipped.shape[0] * yres)
        xOffset = int(round(((lrXx_m - lrXx_s) / geoTrans[1])))
        yOffset = int(round(((lrYy_m - lrYy_s) / geoTrans[5])))

        # get MODIS raster x and y size for resample SMAP pixels
        cols = modisWrap.RasterXSize - xOffset
        rows = modisWrap.RasterYSize - yOffset
        print modisWrap.RasterXSize, modisWrap.RasterYSize, cols, rows
        coords = []
        px = ((geoTrans[5] / 2))
        for i in xrange(cols):

            # loop from MODIS column raster for get x center of each pixel
            if i != 0:
                px = px + geoTrans[5]

            py = (geoTrans[1] / 2)
            for j in xrange(rows):

                # loop from MODIS column raster for get y center of each pixel
                if j != 0:
                    py = py + geoTrans[1]

                # print i, j

                # subtraction shift of center each pixel at x and y coordinate
                x = geoTrans[0] - px
                y = geoTrans[3] - py
                # print x, y

                # calculate each pixel of SMAP clipped that contain n pixel of MODIS raster
                col = int((x - xOrigin) / pixelWidth)
                row = int((yOrigin - y) / pixelHeight)

                # get and add calculated row and column of SMAP that contain n pixel of MODIS raster
                coords.append(clipArray[row][col])
                # print row, col, clipArray[row][col]
                # print
                # print

        # reshape coord list for given (X * Y) for create final raster
        a = np.reshape(coords, (cols, rows))
        print a.shape

        # transpose the coord reshaped list in (rows, cols)
        transpose = a.transpose()
        print transpose.shape

        # create zero array for push gap and create full like array for given Nan value in column xOffset * rows
        zero = np.zeros((xOffset,rows), dtype=np.double)
        fZero = np.full_like(zero, np.nan, dtype=np.double)
        print zero.shape
        print

        # create zero array for push gap and create full like array for given Nan value in rows yOffset * modisWrap.RasterXSize
        zero1 = np.zeros((yOffset,modisWrap.RasterXSize), dtype=np.int)
        fZero1 = np.full_like(zero1, np.nan, dtype=np.double)
        print zero1.shape

        # insert Nan array into final array for end column
        transpose_y = np.insert(transpose, cols, fZero, axis=1)
        print transpose_y.shape
        print
        print

        # insert Nan array into final array for end row
        final = np.insert(transpose_y, rows, fZero1, axis=0)
        print final.shape


        return final



    def calculateMeanSEE(self, SEE=None, clip=None):
        '''
        <SEEMODIS>40 km its average within a SMAP pixel,
        SEE-MODIS the MODIS-derived soil evaporative efficiency (ratio of actual to potential evaporation)

        :param SEE:     MODIS-derived soil evaporative efficiency
        :return:        Soil Evaporative Efficiency at 37KM resolution with (27L, 53L)
        '''

        # Set variables
        meanSEE = []
        SEECols = SEE.shape[1]
        SEERows = SEE.shape[0]
        clipCols = clip.shape[1]
        clipRows = clip.shape[0]
        # print SEECols, SEERows
        # print clipCols, clipRows

        # Set X-Axis for kernel filter
        yy = 0
        yyy = int(round(SEECols / clipCols)) #114
        addy = int(round(SEECols / clipCols))
        # print yyy
        # print

        # Loop through 53 column for SMAP cilpped by MODIS scene
        for row in xrange(clipCols): #53

            # Set Y-Axis for kernel filter
            xx = 0
            xxx = int(round(SEERows / clipRows)) #113
            addx = int(round(SEERows / clipRows))
            minSEERows = SEERows - xxx
            minSEERows_1 = SEERows - (2*xxx)
            # print xxx
            # print minSEERows, minSEERows_1

            # Loop through 114 column within SEE-MODIS array
            for i in xrange(yyy):

                # Loop through 113 column within SEE-MODIS array
                for j in xrange(xxx):

                    # masking invalid values for each column in SEE-MODIS array
                    if (xx >= minSEERows and xxx >= SEERows):
                        break

                    N9 = []
                    # print xx, xxx
                    # print yy, yyy
                    # print
                    # print

                    # Loop through each column and row within each kernel 144*113
                    for x in xrange(xx, xxx):
                        for y in xrange(yy, yyy):
                            # print x, y

                            # masking invalid values for each kernel 144*113
                            if (xx <= minSEERows_1 and xxx <= minSEERows):

                                # Append each selected SEE to list
                                N9.append(SEE[x][y])
                                # print SEE[x][y]

                    # Add 113*113 to kernel window for move kernel window in length SEE-MODIS array columns
                    # addx = xxx - 0
                    xx += addx  #113
                    xxx += addx #113

                    # print N9
                    # print
                    # print

                    # Calculate average for each input list than append it to mean list
                    meanSEE.append(np.nanmean(N9))

            # Add 114*114 to kernel window for move kernel window in length SEE-MODIS array rows
            # addy = yyy - 0
            yy += addy  #114
            yyy += addy #114

        # print meanSEE
        # print
        # print

        # Reshape mean list to create 27L and 53L equal SMAP clipped raster
        mean = np.array(meanSEE).reshape(clipCols, clipRows) # 27 * 53
        print mean.shape
        # print mean

        # transpose the coord reshaped list in (rows, cols)
        transpose_mean = mean.transpose()
        print transpose_mean.shape

        print "<SEEMODIS>40 km its average within a SMAP pixel!"
        return transpose_mean



    def SEEMean(self, modisWrap=None, seeMean=None, transform=None):
        '''

        :param modisWrap:       MODIS wrapped transform
        :param seeMean:
        :param transform:
        :return:
        '''

        # give geotransform from clipped smap raster
        xOrigin = transform[0]
        yOrigin = transform[3]
        pixelWidth = transform[1]
        pixelHeight = -transform[5]

        # get MODIS wrapped transform for calculate center of each pixel in clipped raster
        geoTrans = list(modisWrap.GetGeoTransform())
        lrXx_m = geoTrans[0] + (modisWrap.RasterXSize * geoTrans[1])
        lrYy_m = geoTrans[3] + (modisWrap.RasterYSize * geoTrans[5])

        # Calculate X and Y offset
        ulXx, xres, xskew, ulYy, yskew, yres = self.clippedGeotranse
        lrXx_s = ulXx + (self.smapClipped.shape[1] * xres)
        lrYy_s = ulYy + (self.smapClipped.shape[0] * yres)
        xOffset = int(round(((lrXx_m - lrXx_s) / geoTrans[1])))
        yOffset = int(round(((lrYy_m - lrYy_s) / geoTrans[5])))

        # get MODIS raster x and y size for resample SMAP pixels
        cols = modisWrap.RasterXSize - xOffset
        rows = modisWrap.RasterYSize - yOffset
        print modisWrap.RasterXSize, modisWrap.RasterYSize, cols, rows
        coords = []
        px = ((geoTrans[5] / 2))
        for i in xrange(cols):

            # loop from MODIS column raster for get x center of each pixel
            if i != 0:
                px = px + geoTrans[5]

            py = (geoTrans[1] / 2)
            for j in xrange(rows):

                # loop from MODIS column raster for get y center of each pixel
                if j != 0:
                    py = py + geoTrans[1]

                # print i, j

                # subtraction shift of center each pixel at x and y coordinate
                x = geoTrans[0] - px
                y = geoTrans[3] - py
                # print x, y

                # calculate each pixel of SMAP clipped that contain n pixel of MODIS raster
                col = int((x - xOrigin) / pixelWidth)
                row = int((yOrigin - y) / pixelHeight)

                # get and add calculated row and column of SMAP that contain n pixel of MODIS raster
                coords.append(seeMean[row][col])
                # print row, col, seeMean[row][col]
                # print
                # print

        # reshape coord list for given (X * Y) for create final raster
        a = np.reshape(coords, (cols, rows))
        print a.shape

        # transpose the coord reshaped list in (rows, cols)
        transpose = a.transpose()
        print transpose.shape

        # create zero array for push gap and create full like array for given Nan value in column xOffset * rows
        zero = np.zeros((xOffset, rows), dtype=np.double)
        fZero = np.full_like(zero, np.nan, dtype=np.double)
        print zero.shape
        print

        # create zero array for push gap and create full like array for given Nan value in rows yOffset * modisWrap.RasterXSize
        zero1 = np.zeros((yOffset, modisWrap.RasterXSize), dtype=np.int)
        fZero1 = np.full_like(zero1, np.nan, dtype=np.double)
        print zero1.shape

        # insert Nan array into final array for end column
        transpose_y = np.insert(transpose, cols, fZero, axis=1)
        print transpose_y.shape
        print
        print

        # insert Nan array into final array for end row
        final = np.insert(transpose_y, rows, fZero1, axis=0)
        print final.shape

        return final



    def soilMoistureParameter(self, rescaleSMAP=None, rescaleMeanSEE=None, type=float, ndv=None):
        '''
        A value of SMp is obtained for each SMap pixel and each input data set.
        The soil moisture parameter SMp used to compute *SMmod/*SEE.
        SMp was set to the soil moisture at field capacity. In DisPATCh,
        SMp is retrieved at 40-km resolution from SMAP and aggregated MODIS data.

        :param rescaleSMAP:         SMAP re-scale to MODIS resolution
        :param rescaleMeanSEE:      MeanSEE re-scale to MODIS resolution
        :param ndv:                 No data value
        :return:                    Soil moisture parameter SMp

        '''

        # masking invalid values
        eqn = '(rescaleMeanSEE <= 0) | (rescaleMeanSEE >= 1)'
        invalid = ne.evaluate(eqn)

        # Apply valid range mask to data
        rescaleMeanSEE[invalid] = ndv

        # soil moisture parameter SMp calculator
        pi = np.pi
        acos = np.arccos(1 - 2 * rescaleMeanSEE)
        SMp = ne.evaluate('pi * rescaleSMAP / acos').astype(type)


        print "soil moisture parameter is ok!"
        return SMp



    def calculateTv1km(self, FVC=None, LST=None):

        # masking invalid values
        # LST = LST / 50.0

        # Set Variables
        x = FVC.flatten()
        y = LST.flatten()


        getOffsetlist = zip(x, y)
        getOffset = np.array(getOffsetlist)
        print getOffset.shape

        #################### Ts max ###################

        lstAxisMax = np.nanmax(getOffset, axis=0) - 0.5

        lstMaxIndex = np.where(y > lstAxisMax[1])

        lstMax = []
        for i in lstMaxIndex[0]:
            lstMax.append([x[i], y[i]])


        minlstMax = np.nanmin(lstMax, axis=0)

        Tssmax = []
        for i in lstMax:
            if (i[0] == minlstMax[0]):
                Tssmax.append(i)
        Tsmax = np.mean(Tssmax, axis=0)
        self.Tsmax = Tsmax[1]
        print 'Tsmax'
        print Tsmax


        ##################### Ts min ###################

        lstAxisMin = np.nanmin(getOffset, axis=0) + 0.1

        lstMinIndex = np.where(y < lstAxisMin[1])

        lstMin = []
        for i in lstMinIndex[0]:
            lstMin.append([x[i], y[i]])

        minlstMin = np.nanmin(lstMin, axis=0)

        Tssmin = []
        for i in lstMin:
            if (i[0] == minlstMin[0]):
                Tssmin.append(i)
        Tsmin = np.average(Tssmin, axis=0)
        self.Tsmin = Tsmin[1]
        print 'Tsmin'
        print Tsmin

        ##################### Tv max ##################

        fvcAxisMax = np.nanmax(getOffset, axis=0) - 0.1

        fvcMaxIndex = np.where(x > fvcAxisMax[0])

        fvcMax = []
        for i in fvcMaxIndex[0]:
            fvcMax.append([x[i], y[i]])

        maxArrfvcMax = np.nanmax(fvcMax, axis=0)

        Tvmmax = []
        for i in fvcMax:
            if (i[1] > (maxArrfvcMax[1] - 10)):
                Tvmmax.append(i)

        Tvmax = np.mean(Tvmmax, axis=0)
        print 'Tvmax'
        print Tvmax

        ############################ Tv min ############################

        minArrfvcMax = np.nanmin(fvcMax, axis=0)
        Tvmmin = []
        for i in fvcMax:
            if (i[1] < (minArrfvcMax[1] + 10)):
                Tvmmin.append(i)

        Tvmin = np.mean(Tvmmin, axis=0)
        print 'Tvmin'
        print Tvmin

        ########################## intersect between line A and B #############

        # line segment intersection using vectors

        def perp(a):
            b = np.empty_like(a)
            b[0] = -a[1]
            b[1] = a[0]
            return b

        # line segment a given by endpoints a1, a2
        # line segment b given by endpoints b1, b2
        # return
        def seg_intersect(a1, a2, b1, b2):
            da = a2 - a1
            db = b2 - b1
            dp = a1 - b1
            dap = perp(da)
            denom = np.dot(dap, db)
            num = np.dot(dap, dp)
            return (num / denom.astype(float)) * db + b1

        intersectAB = seg_intersect(Tsmax, Tvmin, Tsmin, Tvmax)
        print 'intersect point between line A and B'
        print intersectAB

        def FindTriangleContainingPoint(pt):

            # Find Endmembers for Tsmax , Tsmin , Tvmax , Tvmin , diameter

            Tsmaxx = Tsmax.tolist()
            Ts_max = Tsmaxx[1]
            Tsminn = Tsmin.tolist()
            Ts_min = Tsminn[1]
            Tvmaxx = Tvmax.tolist()
            Tv_max = Tvmaxx[1]
            Tvminn = Tvmin.tolist()
            Tv_min = Tvminn[1]
            diameter = intersectAB.tolist()
            point = pt.tolist()

            Tsmaxx.append(1)
            Tsminn.append(1)
            Tvmaxx.append(1)
            Tvminn.append(1)
            diameter.append(1)
            point.append(1)

            if (point):

                # def ZoneA(Tsmax, Tsmin, diameter, point):

                a1 = []
                a1.append(point)
                a1.append(Tsminn)
                a1.append(diameter)
                arr_a1 = np.asarray(a1).reshape(3, 3)
                det_arr_a1 = np.linalg.det(arr_a1)

                a2 = []
                a2.append(Tsmaxx)
                a2.append(point)
                a2.append(diameter)
                arr_a2 = np.asarray(a2).reshape(3, 3)
                det_arr_a2 = np.linalg.det(arr_a2)

                a3 = []
                a3.append(Tsmaxx)
                a3.append(Tsminn)
                a3.append(point)
                arr_a3 = np.asarray(a3).reshape(3, 3)
                det_arr_a3 = np.linalg.det(arr_a3)

                if (det_arr_a1 > 0 and det_arr_a2 > 0 and det_arr_a3 > 0):
                    Tv1km_A = (Tv_min + Tv_max) / 2
                    # print 'A'
                    # print Tv1km_A
                    return Tv1km_A

                # def ZoneB(Tsmax, diameter, Tvmax):

                b1 = []
                b1.append(point)
                b1.append(diameter)
                b1.append(Tvmaxx)
                arr_b1 = np.asarray(b1).reshape(3, 3)
                det_arr_b1 = np.linalg.det(arr_b1)

                b2 = []
                b2.append(Tsmaxx)
                b2.append(point)
                b2.append(Tvmaxx)
                arr_b2 = np.asarray(b2).reshape(3, 3)
                det_arr_b2 = np.linalg.det(arr_b2)

                b3 = []
                b3.append(Tsmaxx)
                b3.append(diameter)
                b3.append(point)
                arr_b3 = np.asarray(b3).reshape(3, 3)
                det_arr_b3 = np.linalg.det(arr_b3)

                if (det_arr_b1 > 0 and det_arr_b2 > 0 and det_arr_b3 > 0):
                    Tv1km_B = (Ts_max + Tv_max) / 2
                    # print 'B'
                    # print Tv1km_B
                    return Tv1km_B

                # def ZoneC(Tvmax, diameter, Tvmin):

                c1 = []
                c1.append(point)
                c1.append(Tsminn)
                c1.append(Tvminn)
                arr_c1 = np.asarray(c1).reshape(3, 3)
                det_arr_c1 = np.linalg.det(arr_c1)

                c2 = []
                c2.append(diameter)
                c2.append(point)
                c2.append(Tvminn)
                arr_c2 = np.asarray(c2).reshape(3, 3)
                det_arr_c2 = np.linalg.det(arr_c2)

                c3 = []
                c3.append(diameter)
                c3.append(Tsminn)
                c3.append(point)
                arr_c3 = np.asarray(c3).reshape(3, 3)
                det_arr_c3 = np.linalg.det(arr_c3)

                if (det_arr_c1 > 0 and det_arr_c2 > 0 and det_arr_c3 > 0):
                    Tv1km_C = (Tv_min + Ts_min) / 2
                    # print 'C'
                    # print Tv1km_C
                    return Tv1km_C

                # def ZoneD():

                d1 = []
                d1.append(point)
                d1.append(Tvmaxx)
                d1.append(diameter)
                arr_d1 = np.asarray(d1).reshape(3, 3)
                det_arr_d1 = np.linalg.det(arr_d1)

                d2 = []
                d2.append(Tvminn)
                d2.append(point)
                d2.append(diameter)
                arr_d2 = np.asarray(d2).reshape(3, 3)
                det_arr_d2 = np.linalg.det(arr_d2)

                d3 = []
                d3.append(Tvminn)
                d3.append(Tvmaxx)
                d3.append(point)
                arr_d3 = np.asarray(d3).reshape(3, 3)
                det_arr_d3 = np.linalg.det(arr_d3)

                if (det_arr_d1 > 0 and det_arr_d2 > 0 and det_arr_d3 > 0):
                    Tv1km_D = (Ts_min + Ts_max) / 2
                    # print 'D'
                    # print Tv1km_D
                    return Tv1km_D


            else:
                return None

        Tvs = []
        for i in getOffset:
            Tv = FindTriangleContainingPoint(i)
            Tvs.append(Tv)

        x_pixels = self.w_ds_lst.RasterXSize
        y_pixels = self.w_ds_lst.RasterYSize
        arrTvs = np.array(Tvs).reshape(y_pixels, x_pixels).astype(float)
        print arrTvs.shape


        print "Vegetation Temperature 1 km is ok!"
        return arrTvs



    def resmapleFvc(self, fvc250m=None, lst1km=None):

        # Set variables
        meanSEE = []
        FVCCols = fvc250m.shape[1]
        FVCRows = fvc250m.shape[0]
        lstCols = lst1km.shape[1]
        lstRows = lst1km.shape[0]
        # print SEECols, SEERows
        # print clipCols, clipRows

        # Set X-Axis for kernel filter
        yy = 0
        yyy = int(round(FVCCols / lstCols))  # 114
        addy = int(round(FVCCols / lstCols))
        # print yyy
        # print


        meanFVC = []

        # yy = 0
        # yyy = 5

        for row in xrange(lstCols):

            # Set Y-Axis for kernel filter
            xx = 0
            xxx = int(round(FVCRows / lstRows))  # 113
            addx = int(round(FVCRows / lstRows))
            minFVCRows = FVCRows - xxx
            minFVCRows_1 = FVCRows - (2 * xxx)
            # print xxx
            # print minSEERows, minSEERows_1

            # xx = 0
            # xxx = 5

            for i in xrange(yyy):

                for j in xrange(xxx):

                    if (xx >= minFVCRows and xxx >= FVCRows):
                        break

                    N9 = []
                    # print xx, xxx
                    # print yy, yyy
                    # print
                    # print


                    for x in xrange(xx, xxx):

                        for y in xrange(yy, yyy):

                            # print x,'', y

                            if (xx <= minFVCRows_1 and xxx <= minFVCRows):
                                N9.append(fvc250m[x][y])
                                # print FVC[x][y]

                    # Add 5*5 to kernel window for move kernel window in length SEE-MODIS array columns
                    # addx = xxx - 0
                    xx += addx  # 5
                    xxx += addx  # 5

                    # print N9
                    # print
                    # print

                    meanFVC.append(np.nanmean(N9))

            # Add 5*5 to kernel window for move kernel window in length SEE-MODIS array rows
            # addy = yyy - 0
            yy += addy  # 5
            yyy += addy  # 5


        a = np.reshape(meanFVC, (lstCols, lstRows - 1))
        print a.shape

        zero1 = np.zeros((1, lstCols), dtype=np.double)
        fZero1 = np.full_like(zero1, np.nan, dtype=np.double)
        print zero1.shape

        transpose = a.transpose()
        print transpose.shape

        final = np.insert(transpose, lstRows - 1, fZero1, axis=0).astype(float)
        print final.shape

        print "Resmaple FVC to Fv1km is ok!"
        return final



    def save_tif(self, inputRaster, file_path, fileName):
        '''Save the given Raster as a .tif file.

        Parameters
            inputRaster 	-	  A mask generated with masker.
            file_path       -	  Path of .tif file.
        '''

        # Get the size of the geotiffs
        x_pixels = self.w_ds_red.RasterXSize
        y_pixels = self.w_ds_red.RasterYSize

        # Get transformation/projection info
        geotransform = self.w_ds_red.GetGeoTransform()
        ref = self.w_ds_red.GetProjectionRef()

        # Create new tiff to store a raster stack as a single tif
        driver = gdal.GetDriverByName('GTiff')
        dataset = driver.Create(file_path + fileName + '.tif', x_pixels, y_pixels, 1, gdal.GDT_Float32)

        # Set transform/proj from those found above
        dataset.SetGeoTransform(geotransform)
        dataset.SetProjection(ref)

        # Write Raster To TIFF File
        dataset.GetRasterBand(1).WriteArray(inputRaster)


        print "Save File is ok!"
        dataset.FlushCache()



    def save_tif_1km(self, inputRaster, file_path, fileName):
        '''Save the given Raster as a .tif file.

        Parameters
            inputRaster 	-	  A mask generated with masker.
            file_path       -	  Path of .tif file.
        '''

        # Get the size of the geotiffs
        x_pixels = self.w_ds_lst.RasterXSize
        y_pixels = self.w_ds_lst.RasterYSize

        # Get transformation/projection info
        geotransform = self.w_ds_lst.GetGeoTransform()
        ref = self.w_ds_lst.GetProjectionRef()

        # Create new tiff to store a raster stack as a single tif
        driver = gdal.GetDriverByName('GTiff')
        dataset = driver.Create(file_path + fileName + '.tif', x_pixels, y_pixels, 1, gdal.GDT_Float32)

        # Set transform/proj from those found above
        dataset.SetGeoTransform(geotransform)
        dataset.SetProjection(ref)

        # Write Raster To TIFF File
        dataset.GetRasterBand(1).WriteArray(inputRaster)


        print "Save File is ok!"
        dataset.FlushCache()



    def save_tif_clipSmap(self, inputRaster, file_path, fileName):
        '''Save the given Raster as a .tif file.

        Parameters
            inputRaster 	-	  A mask generated with masker.
            file_path       -	  Path of .tif file.
        '''

        # Get the size of the geotiffs
        x_pixels = self.smapClipped.shape[1]
        y_pixels = self.smapClipped.shape[0]

        # Get transformation/projection info
        # geotransform = self.sm_ds.GetGeoTransform()
        ref = self.sm_ds.GetProjectionRef()

        # Create new tiff to store a raster stack as a single tif
        driver = gdal.GetDriverByName('GTiff')
        dataset = driver.Create(file_path + fileName + '.tif', x_pixels, y_pixels, 1, gdal.GDT_Float32)

        # Set transform/proj from those found above
        dataset.SetGeoTransform(self.clippedGeotranse)
        dataset.SetProjection(ref)

        # Write Raster To TIFF File
        dataset.GetRasterBand(1).WriteArray(inputRaster)

        print "Save File is ok!"
        dataset.FlushCache()
