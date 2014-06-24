load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    root = "../build/Release"
    test_case = "deform"
    ic = "sc"
    resolution = "240x120"
    time_stamp = "%5.5d"
    start_time = 0
    end_time = 600
    time_step = 10

    ; min_lon = 120.0
    ; max_lon = 270.0
    ; min_lat = 10.0
    ; max_lat = 80.0

    file_prefix = "lasm."+test_case+"."+resolution

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '"+root+"/"+file_prefix+"."+time_stamp+".nc\n' $i; done")

    num_time = dimsizes(fs)

    if (ic .ne. "") then
        wks = gsn_open_wks("pdf", file_prefix+"."+ic+".cells")
    else
        wks = gsn_open_wks("pdf", file_prefix+".cells")
    end if

    if (isStrSubset(test_case, "rotation")) then
        gsn_define_colormap(wks, "WhBlGrYeRe")
    end if
    if (isStrSubset(test_case, "barotropic")) then
        if (ic .eq. "gd") then
            gsn_define_colormap(wks, "GMT_panoply")
        end if
        if (ic .eq. "step") then
            cmap = (/(/1.0, 1.0, 1.0 /),(/0.0, 0.0, 0.0 /), \
                     (/0.84,0.84,1.0 /),(/0.38,0.4, 0.98/), \
                     (/0.04,0.15,0.89/),(/0.06,0.4, 0.46/), \
                     (/0.15,0.71,0.14/),(/0.54,0.85,0.17/), \
                     (/0.98,0.98,0.22/),(/0.99,0.53,0.14/), \
                     (/0.99,0.09,0.11/)/)
            gsn_define_colormap(wks, cmap)
        end if
    end if
    if (isStrSubset(test_case, "deform")) then
        if (ic .eq. "ch" .or. ic .eq. "sc") then
            cmap = (/(/1.0, 1.0, 1.0 /),(/0.0, 0.0, 0.0 /), \
                     (/0.84,0.84,1.0 /),(/0.38,0.4, 0.98/), \
                     (/0.04,0.15,0.89/),(/0.06,0.4, 0.46/), \
                     (/0.15,0.71,0.14/),(/0.54,0.85,0.17/), \
                     (/0.98,0.98,0.22/),(/0.99,0.53,0.14/), \
                     (/0.99,0.09,0.11/)/)
            gsn_define_colormap(wks, cmap)
        end if
        if (ic .eq. "gh") then
            gsn_define_colormap(wks, "WhBlGrYeRe")
        end if
    end if

    res = True
    res@gsnFrame = False
    res@cnLinesOn = False
    res@cnFillOn = True
    res@gsnSpreadColors = True
    res@cnLevelSelectionMode = "ManualLevels"
    if (test_case .eq. "barotropic") then
        if (ic .eq. "step") then
            res@cnMinLevelValF = 0.2
            res@cnMaxLevelValF = 0.9
            res@cnLevelSpacingF = 0.1
        else if (ic .eq. "gd") then
            ; res@cnMinLevelValF = 102000.0
            ; res@cnMaxLevelValF = 124000.0
            res@cnMinLevelValF = 61000
            res@cnMaxLevelValF = 86000
            res@cnLevelSpacingF = (res@cnMaxLevelValF-res@cnMinLevelValF)/20.0
        end if
        end if
    end if
    if (isStrSubset(test_case, "rotation")) then
        res@cnMinLevelValF = 0.0
        res@cnMaxLevelValF = 1000.0
        res@cnLevelSpacingF = 50.0
    end if
    if (isStrSubset(test_case, "deform")) then
        if (ic .eq. "ch") then
            res@cnMinLevelValF = 0.2
            res@cnMaxLevelValF = 0.9
            res@cnLevelSpacingF = 0.1
        end if
        if (ic .eq. "sc") then
            res@cnMinLevelValF = 0.2
            res@cnMaxLevelValF = 0.9
            res@cnLevelSpacingF = 0.1
        end if
        if (ic .eq. "gh") then
            res@cnMinLevelValF = 0.0
            res@cnMaxLevelValF = 1.0
            res@cnLevelSpacingF = 0.1
        end if
    end if
    res@mpOutlineOn = False
    res@lbBoxMinorExtentF = 0.2
    res@pmTickMarkDisplayMode = "Always"
    if (.not. isdefined("min_lon")) then
        res@mpCenterLonF = 180.0
    else
        res@mpCenterLonF = (min_lon+max_lon)*0.5
        res@mpLimitMode = "LatLon"
        res@mpMinLonF = min_lon
        res@mpMaxLonF = max_lon
        res@mpMinLatF = min_lat
        res@mpMaxLatF = max_lat
    end if

    do l = 0, num_time-1
        system("echo file: "+fs(l))
        f = addfile(fs(l), "r")
        num_lat = dimsizes(f->lat)
        if (isStrSubset(test_case, "barotropic")) then
            if (ic .eq. "step") then
                q = f->q2/f->q0
            end if
            if (ic .eq. "gd") then
                q = f->q1
            end if
        end if
        if (isStrSubset(test_case, "rotation")) then
            q = f->q1/f->q0
        end if
        if (isStrSubset(test_case, "deform")) then
            if (ic .eq. "ch") then
                q = f->q1/f->q0
            end if
            if (ic .eq. "sc") then
                q = f->q3/f->q0
            end if
            if (ic .eq. "gh") then
                q = f->q4/f->q0
            end if
        end if
        copy_VarCoords(f->q0, q)
        res@gsnLeftString = ""
        res@gsnRightString = ""
        ; res@gsnRightString = "Min:"+sprintf("%8.5f", min(q))+" "+ \
        ;                      "Max:"+sprintf("%8.5f", max(q))
        plot = gsn_csm_contour_map(wks, q, res)

        frame(wks)
    end do

end
