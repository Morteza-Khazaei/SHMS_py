# Satellite-based Hydrological Monitoring System (SHMS)


The DISPATCH (DISaggagregation based on Physical And Theoretical CHange) (Merlin et al., 2008, 2012) is an algorithm that downscales SMAP soil moisture data from 36 km (low resolution) to 1 km resolution (high resolution). This algorithm uses Terra and Aqua satellite data to estimate NDVI and LST twice a day using the Moderate Resolution Imaging Spectroradiometer (MODIS) sensor. These estimations have a resolution of 1 km and can be conducted only if there is no cloud cover. This downscaling process provides the final user with the possibility of estimating soil moisture using remote sensing techniques at high resolution.


website: http://185.105.184.53/en

To use the SHMS_py, you must do the following:

    from SHMS_py import DisPATCh

	source = DisPATCh.DISPATCH(hdfFile=df_MOD09GQ + '\\' + MOD09GQ, MOD11A1=df_MOD11A1 + '\\' + MOD11A1, h5File=fullPath, ndv=None)
	SMhr = source.disaggregated
	source.save_tif_1km(inputRaster=SMhr, file_path=file_dir, fileName=soilMositure_fullname)



