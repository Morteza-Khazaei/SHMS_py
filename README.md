============================================================================================
SHMS_py
============================================================================================

To use the SHMS_py, you must do the following:::

    from SHMS_py import DisPATCh

	source = DisPATCh.DISPATCH(hdfFile=df_MOD09GQ + '\\' + MOD09GQ, MOD11A1=df_MOD11A1 + '\\' + MOD11A1, h5File=fullPath, ndv=None)
	SMhr = source.disaggregated
	source.save_tif_1km(inputRaster=SMhr, file_path=file_dir, fileName=soilMositure_fullname)

============================================================================================
SHMS website: http://185.105.184.53/en
