# Satellite-based Hydrological Monitoring System (SHMS)

Introduction of system
Soil moisture is one of the key parameters in environmental science studies such as water resources management, watershed, agriculture and etc. Field measurement of this parameter on a large scale is very costly and time consuming. In recent years, SMAP and SMOS missions have made it possible to estimation of soil moisture at a global scale. SMAP and SMOS coarse-resolution soil moisture products with 2-3 days of temporal resolution are not usable for water resources management, watershed management, and agriculture at the catchment scale. Therefore, downscaling algorithms have been developed to generating of soil moisture data with better spatial and temporal resolution using SMAP and SMOS soil moisture products. The SHMS system has been developed with the aim of generating soil moisture data with 1-km spatial resolution and 2-3 day temporal resolution using SMAP soil moisture product and MODIS optical products. The downscaling algorithm implemented in this system is the DisPATCh algorithm. This algorithm uses MODIS NDVI and LST products to downscaling of SMAP 36-km soil moisture product and generating new soil moisture data with 1-km spatial resolution. The SHMS system automatically receives new MODIS and SMAP products for every day, then executes the DisPATCh algorithm and sharing the generated soil moisture data through the user interface.

website: http://185.105.184.53/en

To use the SHMS_py, you must do the following:

    from SHMS_py import DisPATCh

	source = DisPATCh.DISPATCH(hdfFile=df_MOD09GQ + '\\' + MOD09GQ, MOD11A1=df_MOD11A1 + '\\' + MOD11A1, h5File=fullPath, ndv=None)
	SMhr = source.disaggregated
	source.save_tif_1km(inputRaster=SMhr, file_path=file_dir, fileName=soilMositure_fullname)



