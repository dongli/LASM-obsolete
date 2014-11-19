load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    root = "../build/Debug"
    test_case = "rotation"
    resolution = "240x120"
    time_stamp = "%5.5d"
    start_time = 0
    end_time = 40
    time_step = 1
    point_size = 0.001
    line_thickness = 1

    PI = atan(1.0)*4
    RAD = PI/180

    file_prefix = "lasm."+test_case+"."+resolution

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '"+root+"/"+file_prefix+"."+time_stamp+".nc\n' $i; done")

    f = addfile(fs(0), "r")

    num_time = dimsizes(fs)
    dims = dimsizes(f->c)
    num_tracer = dims(0)
    delete(dims)
    if (isfilevar(f, "s2")) then
        dims = dimsizes(f->s2)
        num_skel2 = dims(1)
        delete(dims)
        s2 = new((/2,num_skel2+1/), double)
    end if

    wks = gsn_open_wks("pdf", file_prefix+".tracers")

    gsn_define_colormap(wks, "GMT_panoply")
    colormap = gsn_retrieve_colormap(wks)
    levels = flt2dble(fspan(61000, 84000, 15))
    colors = span_color_indexes(colormap, dimsizes(levels))

    res_map = True
    res_map@gsnFrame = False
    res_map@mpCenterLonF = 180.0
    res_map@mpOutlineOn = False
    res_map@mpFillOn = False
    res_map@pmTickMarkDisplayMode = "Always"

    res_tracer = True
    res_tracer@gsMarkerIndex = 16
    res_tracer@gsMarkerSizeF = point_size

    res_s1 = True
    res_s1@gsMarkerIndex = 16
    res_s1@gsMarkerSizeF = point_size
    res_s1@gsMarkerColor = "green"

    res_s2 = True
    res_s2@gsMarkerIndex = 16

    ; check = ispan(15000, 18000, 10)
    ; check_colors = "green"
    ; min_lon = 120.0
    ; max_lon = 270.0
    ; min_lat = 10.0
    ; max_lat = 80.0
    ; min_lon = 0
    ; max_lon = 360
    ; min_lat = -90
    ; max_lat = 90

    do l = 0, num_time-1
        system("echo file: "+fs(l))
        f = addfile(fs(l), "r")
        if (isdefined("id")) then
            delete(id)
        end if
        if (isdefined("c")) then
            delete(c)
        end if
        id = f->id
        c = f->c/RAD
        res_map@gsnPolar = "NH"
        res_map@mpMinLatF = 70
        ; res_map@gsnPolar = "SH"
        ; res_map@mpMaxLatF = -60.0
        ; if (isdefined("check")) then
        ;     i = ind(id .eq. check(0))
        ;     if (ismissing(i)) then
        ;         continue
        ;     end if
        ;     res_map@mpProjection = "Satellite"
        ;     res_map@mpCenterLonF = doubletofloat(c(i,0))
        ;     res_map@mpCenterLatF = doubletofloat(c(i,1))
        ;     res_map@mpLimitMode = "Angles"
        ;     res_map@mpLeftAngleF = 2.0
        ;     res_map@mpRightAngleF = res_map@mpLeftAngleF
        ;     res_map@mpTopAngleF = res_map@mpLeftAngleF
        ;     res_map@mpBottomAngleF = res_map@mpLeftAngleF
        ; end if
        if (isdefined("min_lon")) then
            res_map@mpLimitMode = "LatLon"
            res_map@mpMinLonF = min_lon
            res_map@mpMaxLonF = max_lon
            res_map@mpMinLatF = min_lat
            res_map@mpMaxLatF = max_lat
            ; res_map@mpGridAndLimbOn = True
            ; res_map@mpGridLineColor = "grey"
            ; res_map@mpGridLonSpacingF = 1.5
            ; res_map@mpGridLatSpacingF = 1.5
        end if
        map = gsn_csm_map(wks, res_map)
        ; tracer centroid
        do i = 0, num_tracer-1, 10
            if (isdefined("min_lon")) then
                if (.not. (c(i,0) .ge. min_lon .and. c(i,0) .le. max_lon .and. \
                           c(i,1) .ge. min_lat .and. c(i,1) .le. max_lat)) then
                    continue
                end if
            end if
            if (isdefined("check") .and. any(check .eq. id(i))) then
                res_tracer@gsMarkerColor = "blue"
            else
                res_tracer@gsMarkerColor = "black"
            end if
            gsn_polymarker(wks, map, c(i,0), c(i,1), res_tracer)
        end do
        ; draw the shape within a region
        if (isdefined("min_lon")) then
            res_s2@gsLineColor = "gray"
            res_s2@gsLineThicknessF = line_thickness
            do i = 0, num_tracer-1, 2
                if (c(i,0) .ge. min_lon .and. c(i,0) .le. max_lon .and. \
                    c(i,1) .ge. min_lat .and. c(i,1) .le. max_lat) then
                    do j = 0, num_skel2-1
                        s2(:,j) = f->s2(i,j,:)/RAD
                    end do
                    s2(:,num_skel2) = s2(:,0)
                    ; res_s2@gsLineColor = get_color_index(colormap, levels, f->rho(i,1))
                    gsn_polyline(wks, map, s2(0,:), s2(1,:), res_s2)
                end if
            end do
        end if
        ; tracer shape
        if (.not. isdefined("check")) then
            frame(wks)
            continue
        end if
        res_s2@gsLineColor = "red"
        res_s2@gsLineThicknessF = 2*line_thickness
        do k = 0, dimsizes(check)-1
            i = ind(id .eq. check(k))
            if (ismissing(i)) then
                continue
            end if
            gsn_polymarker(wks, map, f->s1(i,:,0)/RAD, f->s1(i,:,1)/RAD, res_s1)
            do j = 0, num_skel2-1
                s2(:,j) = f->s2(i,j,:)/RAD
            end do
            s2(:,num_skel2) = s2(:,0)
            gsn_polyline(wks, map, s2(0,:), s2(1,:), res_s2)
        end do

        frame(wks)
    end do

end
