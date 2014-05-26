load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    test_case = "deform"
    ic = "sc"
    resolution = "240x120"
    time_stamp = "%4.4d"
    run_type = "Release"
    start_time = 0
    end_time = 300
    time_step = 300

    file_prefix = "tracers."+test_case+"."+resolution

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '../build/"+run_type+"/"+file_prefix+"."+time_stamp+".nc\n' $i; done")

    num_time = dimsizes(fs)

    if (ic .ne. "") then
        wks = gsn_open_wks("pdf", file_prefix+"."+ic+".cells")
    else
        wks = gsn_open_wks("pdf", file_prefix+".cells")
    end if

    if (isStrSubset(test_case, "rotation") .or. \\
        isStrSubset(test_case, "barotropic")) then
        gsn_define_colormap(wks, "WhBlGrYeRe")
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
            ; gsn_define_colormap(wks, "WhBlGrYeRe")
        end if
        if (ic .eq. "gh") then
            gsn_define_colormap(wks, "WhBlGrYeRe")
        end if
    end if

    res = True
    res@cnLinesOn = False
    res@cnFillOn = True
    res@gsnSpreadColors = True
    res@cnLevelSelectionMode = "ManualLevels"
    if (test_case .eq. "barotropic") then
        res@cnMinLevelValF = 75000
        res@cnMaxLevelValF = 85000
        res@cnLevelSpacingF = 500.0
    end if
    if (isStrSubset(test_case, "rotation")) then
        res@cnMinLevelValF = 0.0
        res@cnMaxLevelValF = 1000.0
        res@cnLevelSpacingF = 50.0
    end if
    if (isStrSubset(test_case, "deform")) then
        if (ic .eq. "ch") then
            ; res@cnMinLevelValF = 0.0
            ; res@cnMaxLevelValF = 1.0
            ; res@cnLevelSpacingF = 0.1
            res@cnMinLevelValF = 0.2
            res@cnMaxLevelValF = 0.9
            res@cnLevelSpacingF = 0.1
        end if
        if (ic .eq. "sc") then
            res@cnMinLevelValF = 0.2
            res@cnMaxLevelValF = 0.9
            res@cnLevelSpacingF = 0.1
            ; res@cnMinLevelValF = 0.0
            ; res@cnMaxLevelValF = 1.1
            ; res@cnLevelSpacingF = 0.1
        end if
        if (ic .eq. "gh") then
            res@cnMinLevelValF = 0.0
            res@cnMaxLevelValF = 1.0
            res@cnLevelSpacingF = 0.1
        end if
    end if
    res@mpCenterLonF = 180.0
    res@mpOutlineOn = False
    res@pmTickMarkDisplayMode = "Always"
    res@lbBoxMinorExtentF = 0.2

    do l = 0, num_time-1
        system("echo file: "+fs(l))
        f = addfile(fs(l), "r")
        num_lat = dimsizes(f->lat)
        if (isStrSubset(test_case, "barotropic")) then
            q = f->q1
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
        res@gsnRightString = "Min:"+sprintf("%8.5f", min(q))+" "+ \
                             "Max:"+sprintf("%8.5f", max(q))
        plot = gsn_csm_contour_map(wks, q, res)
    end do

end
