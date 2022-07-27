;+
; :Description:
;    IDL interface for reading cpt file format file
; :Version:
;    v0.1 # init: May/23/2022 last modify Jul/27/2022
; :Params: fname
; :Keywords:
;    nptx:   variable receiving the value of count of Ptx
;    nparam: variable receiving the value of count of param inside Point
; :Returns: list of Ptx
; :Author: jtsung
;-
function cpt_readall, fname, nptx=ptxcount, nparam=paramcount
	compile_opt idl2, hidden
	
	;exit if file not exist
	if ~FILE_TEST(fname, /READ, /REGULAR) then begin
		PRINTF, -1, "Can NOT open ", fname
		RETURN, 0
	endif
	
	;open for read
	OPENR, lun, fname, /GET_LUN
	
	;magic number
	READU, lun, (mgc = BYTARR(13))
	realmgc = [   2,  76, 101, 114,  83,  65,  84,  64,  99, 112, 116,  10,   3]
	if ~ARRAY_EQUAL(mgc, realmgc) then begin
		PRINTF, -1, fname, " is NOT a cpt file!"
		FREE_LUN, lun
		RETURN, 0
	endif
	
	;version
	READU, lun, (ver = 0B)
	realmajor = 0B
	realminor = 1B
	realver = ISHFT(realmajor, 4) or realminor
	if ~(realver eq ver) then begin
		PRINTF, -1, "version NOT match!"
		FREE_LUN, lun
		RETURN, 0
	endif
	
	;meta
	READU, lun, (ptxcount = 0UL)
	READU, lun, (paramcount = 0B)
	
	;main data envelope
	ret = LIST(LENGTH=ptxcount)
	for iptx = 0UL, ptxcount-1 do begin
		
		;Pt
		READU, lun, (lon = 0.)
		READU, lun, (lat = 0.)
		READU, lun, (alt = 0US)
		READU, lun, (npoint = 0B)
		points = LIST(LENGTH=npoint)
		for ipoint = 0UL, npoint-1 do begin
			READU, lun, (time = FIX(0, TYPE=15))
			READU, lun, (param = MAKE_ARRAY(paramcount, /DOUBLE))
			points[ipoint] = {time: time, data: param}
		endfor
		pt = {lon: lon, lat: lat, alt: alt, point: points}
		
		;Px
		READU, lun, (time = FIX(0, TYPE=15))
		pixel = readpixel(lun)
		READU, lun, (nviciniy = 0B)
		vicinity = LIST(LENGTH=nviciniy)
		for ivicinity = 0UL, nviciniy-1 do vicinity[ivicinity] = readpixel(lun)
		px = {time: time, center: pixel, vicinity: vicinity}
		
		ret[iptx] = {pt: pt, px: px}
	endfor

	;ending
	READU, lun, (ending = MAKE_ARRAY(16, /BYTE))
	FREE_LUN, lun
	if ~ARRAY_EQUAL(ending, 0B) then begin
		PRINTF, -1, fname, " has NO ending, the results may be incorrect"
	endif
	
	RETURN, ret
end

function readpixel, lun
	compile_opt idl2, hidden

	READU, lun, (lon = 0.)
	READU, lun, (lat = 0.)
	READU, lun, (alt = 0US)
	READU, lun, (mask = 0B)
	READU, lun, (nchannel = 0B)
	READU, lun, (nlayer = 0B)
	
	channels = LIST(LENGTH=nchannel)
	for ichannel = 0UL, nchannel-1 do begin
		READU, lun, (centrewv = 0S)
		
		;obs
		READU, lun, (i = MAKE_ARRAY(nlayer, /DOUBLE))
		if centrewv lt 0 then begin
			READU, lun, (q = MAKE_ARRAY(nlayer, /DOUBLE))
			READU, lun, (u = MAKE_ARRAY(nlayer, /DOUBLE))
		endif
		
		;angles
		READU, lun, (sza = MAKE_ARRAY(nlayer, /DOUBLE))
		READU, lun, (vza = MAKE_ARRAY(nlayer, /DOUBLE))
		READU, lun, (saa = MAKE_ARRAY(nlayer, /DOUBLE))
		READU, lun, (vaa = MAKE_ARRAY(nlayer, /DOUBLE))
		
		;append channel
		if centrewv lt 0 then $
			channels[ichannel] = {i: i, q: q, u: u, sza: sza, vza: vza, saa: saa, vaa: vaa} $
		else $
			channels[ichannel] = {i: i, sza: sza, vza: vza, saa: saa, vaa: vaa}
	endfor
	
	READU, lun, (nextra = 0B)
	if nextra then begin
		READU, lun, (extra = MAKE_ARRAY(nextra, /DOUBLE))
		RETURN, {lon: lon, lat: lat, alt: alt, mask: mask, bands: channels, extra: extra}
	endif else begin
		RETURN, {lon: lon, lat: lat, alt: alt, mask: mask, bands: channels}
	endelse
end
