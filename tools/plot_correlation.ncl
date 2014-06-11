load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

function correlation_function(q1)
begin

	a = -0.8
	b =  0.9
	return a*q1^2+b

end

function distance(x0, y0, x1, y1, Rx, Ry)
begin

	return sqrt(((x1-x0)/Rx)^2+((y1-y0)/Ry)^2)

end

begin

    check_parcel = True
    root = "../build/Release"
	file_prefix = "lasm.deform.case4.240x120"

    f = addfile(root+"/"+file_prefix+".00300.nc", "r")

    if (check_parcel) then
        q1 = f->rho(:,1)/f->rho(:,0)
        q2 = f->rho(:,2)/f->rho(:,0)
        dA = new(dimsizes(f->id), double)
        do i = 0, dimsizes(f->id)-1
            dA(i) = f->h(i,0,0)*f->h(i,1,1)-f->h(i,0,1)*f->h(i,1,0)
        end do
    else
        q1 = ndtooned(f->q1/f->q0)
        q2 = ndtooned(f->q2/f->q0)
        dA = ndtooned(f->volume)
    end if
	A0 = sum(dA)

    ; reference curve
    a = -0.8
    b =  0.9
    refer_q1 = fspan(0.1, 1.0, 100)
    refer_q2 = correlation_function(refer_q1)

    ; calcuate the mixing diagnostics
    q1_min = 0.1
    q1_max = 1.0
    q2_min = correlation_function(q1_min)
    q2_max = correlation_function(q1_max)
    R_q1 = q1_max-q1_min
    R_q2 = q2_max-q2_min
    eps = 1.0d-7
    l_r = 0.0d0 ; 'real' mixing diagnostic
    l_u = 0.0d0 ; unmixing diagnostic
    l_o = 0.0d0 ; overshooting diagnostic
    n = dimsizes(q1)
    do i = 0, n-1
    	q2_upper = correlation_function(q1(i))
    	q2_lower = R_q2/R_q1*(q1(i)-q1_min)+q2_min
		sqrt_arg = 12.0d0*(125.0d0*q2(i)-52.0)^3+29648025.0d0*q1(i)^2
		if (sqrt_arg .lt. 0) then
			system("echo [Error]: Can not find closest point!")
			exit
		end if
		c = (65340.0d0*q1(i)+12.0d0*sqrt(sqrt_arg))^(1.0d0/3.0d0)/60.0d0
		Q1 = c-(-13.0d0/75.0d0+5.0d0/12.0d0*q2(i))/c
		Q1 = min((/1.0d0,max((/0.1d0,Q1/))/))
		Q2 = correlation_function(Q1)
    	if (q2(i) .lt. q2_upper+eps .and. q2(i) .gt. q2_lower-eps) then
			; 'real' mixing region
			l_r = l_r+distance(q1(i), q2(i), Q1, Q2, R_q1, R_q2)*dA(i)
    	else if (q1(i) .lt. q1_max+eps .and. q1(i) .gt. q1_min-eps .and. \
    	         q2(i) .lt. q2_min+eps .and. q2(i) .gt. q2_max-eps) then
			; unmixing region
			l_u = l_u+distance(q1(i), q2(i), Q1, Q2, R_q1, R_q2)*dA(i)
    	else
    		; overshooting region
    		l_o = l_o+distance(q1(i), q2(i), Q1, Q2, R_q1, R_q2)*dA(i)
    	end if
    	end if
	end do
	l_r = l_r/A0
	l_u = l_u/A0
	l_o = l_o/A0

    print(l_r)
    print(l_u)
    print(l_o)

    wks = gsn_open_wks("pdf", file_prefix+".correlation")

    res_xy = True
    res_xy@gsnDraw = False
    res_xy@gsnFrame = False
    res_xy@vpWidthF = 0.6
    res_xy@vpHeightF = 0.4
    res_xy@tmXMajorGrid = True
    res_xy@tmYMajorGrid = True
    res_xy@tmXMajorGridThicknessF = 1.0
    res_xy@tmYMajorGridThicknessF = 1.0
    res_xy@tmXMajorGridLineDashPattern = 2
    res_xy@tmYMajorGridLineDashPattern = 2
    res_xy@trXMinF = 0
    res_xy@trXMaxF = 1.1
    res_xy@trYMinF = 0
    res_xy@trYMaxF = 1
    res_xy@tmXBTickStartF = 0
    res_xy@tmXBTickEndF = 1.1
    res_xy@tmYLTickStartF = 0
    res_xy@tmYLTickEndF = 1

    res_xy@xyMarkLineMode = "Lines"
    res_xy@xyLineThicknessF = 3.0
    plot0 = gsn_csm_xy(wks, refer_q1, refer_q2, res_xy)

    res_line = True
    res_line@gsLineThicknessF = 3
    gsn_polyline(wks, plot0, (/0.1,1.0,1.0,0.1/), (/0.892,0.892,0.1,0.892/), res_line)

    res_xy@xyMarkLineMode = "Markers"
    res_xy@xyMarker = 1
    res_xy@xyMarkerColor = "red"
    plot1 = gsn_csm_xy(wks, q1, q2, res_xy)

    res_text = True
    res_text@txPerimOn = True
    res_text@txFont = 29
    res_text@txFontHeightF = 0.015
    res_text@txBackgroundFillColor = "white"
    diag_text = "l~B~r~N~ = "+sprintf("%4.2f", l_r*1.0d4)+"*10~S~-4~N~~C~"+ \
                "l~B~u~N~ = "+sprintf("%4.2f", l_u*1.0d7)+"*10~S~-7~N~~C~"+ \
                "l~B~o~N~ = "+sprintf("%4.2f", l_o)

    res_anno = True
    res_anno@amParallelPosF = -0.30
    res_anno@amOrthogonalPosF = 0.35

    legend = gsn_create_text(wks, diag_text, res_text)
    anno = gsn_add_annotation(plot0, legend, res_anno)

    overlay(plot0, plot1)
    draw(plot0)
    frame(wks)

end