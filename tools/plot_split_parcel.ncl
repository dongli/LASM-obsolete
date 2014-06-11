load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    load "/Users/dongli/OneDrive/LASM/build/Debug/split_parcel.txt"

    PI = atan(1.0)*4
    RAD = PI/180
    point_size = 0.01
    line_thickness = 3
    satellite_angle = 3

    wks = gsn_open_wks("pdf", "split_parcel")

    res_map = True
    res_map@gsnMaximize = False
    res_map@gsnFrame = False
    res_map@mpCenterLonF = 180.0
    res_map@mpProjection = "Satellite"
    if (isdefined("centroid3")) then
        res_map@mpCenterLonF = centroid0(0)
        res_map@mpCenterLatF = centroid0(1)
    else
        res_map@mpCenterLonF = (centroid0(0)+centroid1(0))*0.5
        res_map@mpCenterLatF = (centroid0(1)+centroid1(1))*0.5
    end if
    res_map@mpLimitMode = "Angles"
    res_map@mpLeftAngleF = satellite_angle
    res_map@mpRightAngleF = satellite_angle
    res_map@mpTopAngleF = satellite_angle
    res_map@mpBottomAngleF = satellite_angle

    res_map@mpFillOn = True
    res_map@mpOutlineOn = True
    res_map@mpGridAndLimbOn = True
    res_map@mpGridLineColor = "gray"
    res_map@mpGridLonSpacingF = 1.5
    res_map@mpGridLatSpacingF = 1.5
    map = gsn_csm_map(wks, res_map)

    res_shape = True
    res_shape@gsMarkerIndex = 16
    res_shape@gsMarkerSizeF = point_size
    res_shape@gsLineThicknessF = line_thickness

    res_shape@gsLineColor = "red"
    gsn_polyline(wks, map, shape0(:,0), shape0(:,1), res_shape)
    if (isdefined("centroid3")) then
        res_shape@gsLineColor = "blue"
        gsn_polyline(wks, map, shape1(:,0), shape1(:,1), res_shape)
        res_shape@gsLineColor = "gray"
        gsn_polyline(wks, map, shape2(:,0), shape2(:,1), res_shape)
        res_shape@gsLineColor = "gray"
        gsn_polyline(wks, map, shape3(:,0), shape3(:,1), res_shape)
    else
        res_shape@gsLineColor = "gray"
        gsn_polyline(wks, map, shape1(:,0), shape1(:,1), res_shape)
        res_shape@gsLineColor = "gray"
        gsn_polyline(wks, map, shape2(:,0), shape2(:,1), res_shape)
    end if

    frame(wks)

end
