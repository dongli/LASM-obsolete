load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

	root = "../build/Release"
	file_prefix = "lasm.deform.case4.240x120"

	f0 = addfile(root+"/"+file_prefix+".00000.nc", "r")
	f1 = addfile(root+"/"+file_prefix+".00300.nc", "r")

	num_lon = dimsizes(f0->lon)
	num_lat = dimsizes(f1->lat)

	t = fspan(0.1, 1, 19)
	A0 = new(dimsizes(t), double)
	A1 = new(dimsizes(t), double)
	lf = new(dimsizes(t), double)

	q0 = f0->q1/f0->q0
	q1 = f1->q1/f1->q0
	a = f0->volume

	A0 = 0
	A1 = 0
	do j = 0, num_lat-1
		do i = 0, num_lon-1
			do k = 0, dimsizes(t)-1
				if (q0(j,i) .ge. t(k)) then
					A0(k) = A0(k)+a(j,i)
				end if
				if (q1(j,i) .ge. t(k)) then
					A1(k) = A1(k)+a(j,i)
				end if
			end do
		end do
	end do

	lf = 0
	do k = 0, dimsizes(t)-1
		if (A0(k) .ne. 0) then
			lf(k) = A1(k)/A0(k)
		end if
	end do
	lf = lf*100

	wks = gsn_open_wks("pdf", file_prefix+".filament")

	res = True
	res@vpWidthF = 0.6
	res@vpHeightF = 0.4
	res@trXMinF = 0.1
	res@trXMaxF = 1.0
	res@trYMinF = 0
	res@trYMaxF = 140
	res@tmXBTickStartF = 0.1
	res@tmXBTickEndF = 1.0
	res@tmYLTickStartF = 0
	res@tmYLTickEndF = 140
	res@tmXBMode = "Manual"
	res@tmXBTickSpacingF = 0.1
	res@xyMarkLineMode = "MarkLines"
	res@tmXMajorGrid = True
	res@tmYMajorGrid = True
	res@tmXMajorGridThicknessF = 0.1
	res@tmYMajorGridThicknessF = 0.1
	res@tmXMajorGridLineDashPattern = 2
	res@tmYMajorGridLineDashPattern = 2
	res@xyLineThicknessF = 3
	res@xyMarkerSizeF = 10

	plot = gsn_csm_xy(wks, t, lf, res)

end
