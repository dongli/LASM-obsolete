load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

function plot_cell(wks, test_case, ic, time, data)
begin

    ; Define colormap for each test case.
    if (isStrSubset(test_case, "rotation")) then
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
        end if
        if (ic .eq. "gh") then
            gsn_define_colormap(wks, "WhBlGrYeRe")
        end if
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
    if (isStrSubset(test_case, "terminator_chemistry")) then
        gsn_define_colormap(wks,"wh-bl-gr-ye-re")
    end if

    ; Set plot resources.
    res = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@gsnLeftStringOrthogonalPosF = 0
    res@gsnRightStringOrthogonalPosF = 0
    res@gsnLeftString = "Test tracer: "+ic
    res@cnLinesOn = False
    res@cnFillOn = True
    res@gsnSpreadColors = True
    res@cnConstFEnableFill = True
    res@cnConstFLabelOn = False
    res@cnLevelSelectionMode = "ManualLevels"
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
    if (isStrSubset(test_case, "barotropic")) then
        if (ic .eq. "step") then
            res@cnMinLevelValF = 0.2
            res@cnMaxLevelValF = 0.9
            res@cnLevelSpacingF = 0.1
        else if (ic .eq. "gd") then
            if (isStrSubset(test_case, "case1")) then
                res@cnMinLevelValF = 61000
                res@cnMaxLevelValF = 86000
                res@cnLevelSpacingF = (res@cnMaxLevelValF-res@cnMinLevelValF)/20.0
            else if (isStrSubset(test_case, "case2")) then
                res@cnMinLevelValF = 102000
                res@cnMaxLevelValF = 124000
                res@cnLevelSpacingF = (res@cnMaxLevelValF-res@cnMinLevelValF)/20.0
            end if
            end if
        end if
        end if
    end if
    if (isStrSubset(test_case, "terminator_chemistry")) then
        res@gsnSpreadColorStart = 10
        res@gsnSpreadColorEnd = -4
        res@cnLevelSelectionMode = "ExplicitLevels"
        if (ic .eq. "cl") then
            res@gsnLeftString = "Cl"
            res@cnLevels = (/-1e-7, -1e-8, 1e-8, 1e-7, 1e-6, 2e-6, 3e-6, 3.9e-6, 3.99e-6, 4.01e-6, 4.1e-6/)
        end if
        if (ic .eq. "cl2") then
            res@gsnLeftString = "2*Cl2"
            res@cnLevels = (/-1e-7, -1e-8, 1e-8, 1e-7, 1e-6, 2e-6, 3e-6, 3.9e-6, 3.99e-6, 4.01e-6, 4.1e-6/)
        end if
        if (ic .eq. "cly") then
            res@gsnLeftString = "Cly = Cl + 2*Cl2"
            res@cnLevels = (/3.e-6,  3.7e-6, 3.9e-6, 3.97e-6, 3.99e-6, 4.01e-6, 4.03e-6, 4.1e-6, 4.3e-6, 5e-6 /)
        end if
        res@lbLabelStrings  = sprintf("%4.2fe-6", 1e6*res@cnLevels)
        res@lbOrientation = "Vertical"

        ;************************************************
        ;  Define terminator (boundary of photolysis region)
        ;************************************************

        pi = 4.0*atan(1.0)
        rtd = 180.0/pi  ; radians to degrees
        k1_lat_center =  pi* 20.0/180.0
        k1_lon_center =  pi*300.0/180.0
        longs = fspan(-pi,pi,361)
        latit = rtd*atan( -cos(longs-k1_lon_center)/tan(k1_lat_center) )
        latit!0 = "longs"
        longs = longs*rtd

        ;************************************************
        ;  Characterize and plot terminator 
        ;************************************************

        resxy = True
        resxy@gsnDraw = False
        resxy@gsnFrame = False
        resxy@xyLineThicknesses = 3.0
        resxy@xyLineColors = "black"

        term_plot  = gsn_csm_xy(wks, longs, latit, resxy)
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

    res@gsnRightString = cd_string(time, "%Y-%N-%D %H-%M-%S")
    ; res@gsnRightString = "Min:"+sprintf("%8.5f", min(q))+" "+ \
    ;                      "Max:"+sprintf("%8.5f", max(q))
    ; q = smth9_Wrap(q, 0.5, 0.25, True)
    plot = gsn_csm_contour_map(wks, data(0,:,:), res)

    if (isStrSubset(test_case, "terminator_chemistry")) then
        overlay(plot, term_plot)
    end if

    return plot
end

begin

    ; terminator_chemistry, deform, barotropic
    root = "."
    test_case = "barotropic.case1"
    ic = "gd"
    resolution = "240x120"
    time_stamp = "%5.5d"
    start_time = 3150
    end_time = 3690
    time_step = 90

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

        if (isStrSubset(test_case, "terminator_chemistry")) then
            q = f->q1/f->q0
            copy_VarCoords(f->q0, q)
            plot0 = plot_cell(wks, test_case, "cl", f->time, q)
            q = 2*f->q2/f->q0
            copy_VarCoords(f->q0, q)
            plot1 = plot_cell(wks, test_case, "cl2", f->time, q)
            q = (f->q1+2*f->q2)/f->q0
            copy_VarCoords(f->q0, q)
            plot2 = plot_cell(wks, test_case, "cly", f->time, q)

            gsn_panel(wks, (/plot0,plot1,plot2/), (/3,1/), False)
        else
            copy_VarCoords(f->q0, q)
            plot = plot_cell(wks, test_case, ic, f->time, q)
            draw(plot)
            frame(wks)
        end if
    end do

end
