load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    load "../build/Debug/disorder_tracers.txt"

    PI = atan(1.0)*4
    RAD = PI/180
    point_size = 0.01
    line_thickness = 1
    satellite_angle = 4

    wks = gsn_open_wks("pdf", "merge_parcel")

    res_map = True
    res_map@gsnMaximize = False
    res_map@gsnFrame = False
    res_map@mpCenterLonF = 180.0
    res_map@mpProjection = "Satellite"
    res_map@mpCenterLonF = centroid0(0)
    res_map@mpCenterLatF = centroid0(1)
    res_map@mpLimitMode = "Angles"
    res_map@mpLeftAngleF = satellite_angle
    res_map@mpRightAngleF = satellite_angle
    res_map@mpTopAngleF = satellite_angle
    res_map@mpBottomAngleF = satellite_angle

    res_map@mpFillOn = True
    res_map@mpOutlineOn = True
    res_map@mpGridAndLimbOn = True
    res_map@mpGridLineColor = "gray"
    res_map@mpGridLonSpacingF = 2.25
    res_map@mpGridLatSpacingF = 2.25
    map = gsn_csm_map(wks, res_map)

    res_shape = True
    res_shape@gsMarkerIndex = 16
    res_shape@gsMarkerSizeF = point_size
    res_shape@gsLineThicknessF = line_thickness
    ; res_shape@gsLineColor = "red"

    res_cell = True
    res_cell@gsMarkerIndex = 16
    res_cell@gsMarkerSizeF = 0.005
    res_cell@gsMarkerColor = "blue"

    if (isdefined("weights")) then
	    gsn_define_colormap(wks, "MPL_Reds")
	    colormap = gsn_retrieve_colormap(wks)
	    levels = fspan(0, 1, 101)
	    colors = span_color_indexes(colormap(3:120,:), dimsizes(levels))
    end if

    if (isdefined("shape1")) then
	    if (isdefined("weights")) then
	        res_shape@gsLineColor = colors(closest_val(weights(0), levels))+3
        end if
        gsn_polyline(wks, map, shape1(:,0), shape1(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells1(:,0), ngb_cells1(:,1), res_cell)
    end if
    if (isdefined("shape2")) then
	    if (isdefined("weights")) then
	        res_shape@gsLineColor = colors(closest_val(weights(1), levels))+3
        end if
        gsn_polyline(wks, map, shape2(:,0), shape2(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells2(:,0), ngb_cells2(:,1), res_cell)
    end if
    if (isdefined("shape3")) then
	    if (isdefined("weights")) then
	        res_shape@gsLineColor = colors(closest_val(weights(2), levels))+3
	    end if
        gsn_polyline(wks, map, shape3(:,0), shape3(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells3(:,0), ngb_cells3(:,1), res_cell)
    end if
    if (isdefined("shape4")) then
	    if (isdefined("weights")) then
	        res_shape@gsLineColor = colors(closest_val(weights(3), levels))+3
	    end if
        gsn_polyline(wks, map, shape4(:,0), shape4(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells4(:,0), ngb_cells4(:,1), res_cell)
    end if
    if (isdefined("shape5")) then
	    if (isdefined("weights")) then
	        res_shape@gsLineColor = colors(closest_val(weights(4), levels))+3
        end if
        gsn_polyline(wks, map, shape5(:,0), shape5(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells5(:,0), ngb_cells5(:,1), res_cell)
    end if
    if (isdefined("shape6")) then
	    if (isdefined("weights")) then
	        res_shape@gsLineColor = colors(closest_val(weights(5), levels))+3
        end if
        gsn_polyline(wks, map, shape6(:,0), shape6(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells6(:,0), ngb_cells6(:,1), res_cell)
    end if
    if (isdefined("shape7")) then
        if (isdefined("weights")) then
            res_shape@gsLineColor = colors(closest_val(weights(6), levels))+3
        end if
        gsn_polyline(wks, map, shape7(:,0), shape7(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells7(:,0), ngb_cells7(:,1), res_cell)
    end if
    if (isdefined("shape8")) then
	    if (isdefined("weights")) then
	        res_shape@gsLineColor = colors(closest_val(weights(7), levels))+3
        end if
        gsn_polyline(wks, map, shape8(:,0), shape8(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells8(:,0), ngb_cells8(:,1), res_cell)
    end if
    if (isdefined("shape9")) then
	    if (isdefined("weights")) then
	        res_shape@gsLineColor = colors(closest_val(weights(8), levels))+3
        end if
        gsn_polyline(wks, map, shape9(:,0), shape9(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells9(:,0), ngb_cells9(:,1), res_cell)
    end if
    if (isdefined("shape10")) then
	    if (isdefined("weights")) then
	        res_shape@gsLineColor = colors(closest_val(weights(9), levels))+3
        end if
        gsn_polyline(wks, map, shape10(:,0), shape10(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells10(:,0), ngb_cells10(:,1), res_cell)
    end if
    if (isdefined("shape11")) then
	    if (isdefined("weights")) then
	        res_shape@gsLineColor = colors(closest_val(weights(10), levels))+3
        end if
        gsn_polyline(wks, map, shape11(:,0), shape11(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells11(:,0), ngb_cells11(:,1), res_cell)
    end if
    if (isdefined("shape12")) then
	    if (isdefined("weights")) then
	        res_shape@gsLineColor = colors(closest_val(weights(11), levels))+3
        end if
        gsn_polyline(wks, map, shape12(:,0), shape12(:,1), res_shape)
    end if
    if (isdefined("shape13")) then
	    if (isdefined("weights")) then
	        res_shape@gsLineColor = colors(closest_val(weights(12), levels))+3
        end if
        gsn_polyline(wks, map, shape13(:,0), shape13(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells13(:,0), ngb_cells13(:,1), res_cell)
    end if
    if (isdefined("shape14")) then
	    if (isdefined("weights")) then
	        res_shape@gsLineColor = colors(closest_val(weights(13), levels))+3
        end if
        gsn_polyline(wks, map, shape14(:,0), shape14(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells14(:,0), ngb_cells14(:,1), res_cell)
    end if
    if (isdefined("shape15")) then
	    if (isdefined("weights")) then
        	res_shape@gsLineColor = colors(closest_val(weights(14), levels))+3
    	end if
        gsn_polyline(wks, map, shape15(:,0), shape15(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells15(:,0), ngb_cells15(:,1), res_cell)
    end if
    if (isdefined("shape16")) then
	    if (isdefined("weights")) then
        	res_shape@gsLineColor = colors(closest_val(weights(15), levels))+3
    	end if
        gsn_polyline(wks, map, shape16(:,0), shape16(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells16(:,0), ngb_cells16(:,1), res_cell)
    end if
    if (isdefined("shape17")) then
	    if (isdefined("weights")) then
        	res_shape@gsLineColor = colors(closest_val(weights(16), levels))+3
    	end if
        gsn_polyline(wks, map, shape17(:,0), shape17(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells17(:,0), ngb_cells17(:,1), res_cell)
    end if
    if (isdefined("shape18")) then
	    if (isdefined("weights")) then
        	res_shape@gsLineColor = colors(closest_val(weights(17), levels))+3
    	end if
        gsn_polyline(wks, map, shape18(:,0), shape18(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells18(:,0), ngb_cells18(:,1), res_cell)
    end if
    if (isdefined("shape19")) then
	    if (isdefined("weights")) then
        	res_shape@gsLineColor = colors(closest_val(weights(18), levels))+3
    	end if
        gsn_polyline(wks, map, shape19(:,0), shape19(:,1), res_shape)
        gsn_polymarker(wks, map, ngb_cells19(:,0), ngb_cells19(:,1), res_cell)
    end if

    delete(res_shape@gsLineColor)
    res_shape@gsLineColor = "blue"
    gsn_polyline(wks, map, shape0(:,0), shape0(:,1), res_shape)
    gsn_polymarker(wks, map, ngb_cells0(:,0), ngb_cells0(:,1), res_cell)

    if (isdefined("vertices")) then
        res_point = True
        res_point@gsMarkerIndex = 16
        res_point@gsMarkerColor = "green"
        gsn_polymarker(wks, map, vertices(:,0), vertices(:,1), res_point)
    end if

    frame(wks)

end
