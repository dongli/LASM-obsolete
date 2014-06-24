load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

	file_prefix = "barotropic-output"
    run_type = "Release"
    time_stamp = "%5.5d"
    start_time = 4700
    end_time = 4700
    time_step = 10

    min_lon = 120.0
    max_lon = 270.0
    min_lat = 10.0
    max_lat = 80.0

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '../build/"+run_type+"/"+file_prefix+"."+time_stamp+".nc\n' $i; done")

    num_time = dimsizes(fs)

    wks = gsn_open_wks("pdf", file_prefix)

    gsn_define_colormap(wks, "GMT_panoply")

    res = True
    res@gsnFrame = False
    res@cnLinesOn = False
    res@cnFillOn = True
    res@cnLevelSelectionMode = "ManualLevels"
    res@gsnSpreadColors = True
    res@mpCenterLonF = 180.0
    res@mpOutlineOn = False
    res@pmTickMarkDisplayMode = "Always"
    res@lbBoxMinorExtentF = 0.2
    if (isdefined("min_lon")) then
        res@mpLimitMode = "LatLon"
        res@mpMinLonF = min_lon
        res@mpMaxLonF = max_lon
        res@mpMinLatF = min_lat
        res@mpMaxLatF = max_lat
    end if

    do l = 0, num_time-1, 1
        system("echo file: "+fs(l))
        f = addfile(fs(l), "r")

        ; res@cnMinLevelValF = -2.0e-5
        ; res@cnMaxLevelValF =  2.0e-5
        ; res@cnLevelSpacingF = 2.0e-6
        ; plot = gsn_csm_contour_map(wks, f->div(0,:,:), res)
        ; frame(wks)

        res@cnMinLevelValF = -1.0e-4
        res@cnMaxLevelValF =  1.0e-4
        res@cnLevelSpacingF = 1.0e-6
        plot = gsn_csm_contour_map(wks, f->vor_xy(0,:,:), res)
        frame(wks)

        ; res@cnMinLevelValF = 62000.0
        ; res@cnMaxLevelValF = 84000.0
        ; res@cnLevelSpacingF = (res@cnMaxLevelValF-res@cnMinLevelValF)/20.0
        ; plot = gsn_csm_contour_map(wks, f->gd(0,:,:), res)
        ; frame(wks)
	end do

end