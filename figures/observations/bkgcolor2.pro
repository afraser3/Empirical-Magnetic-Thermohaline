pro bkgcolor2, xmin, xmax, ymin, ymax, nxboxes, nyboxes, x, y, data, dataerr, KICID=KICID

;Jamie Tayar OSU 6/25/2013
;color in the background as the average of the colors of the points in 
;  that section of the plot
; I think nxboxes, nyboxes should be floats??
; colortable must be defined before calling bkgcolor

;device, /close
             openw, unit1, 'fig5.dat', /get_lun, /append
for i=0, nxboxes-1 do begin
   for j=0, nyboxes-1 do begin

        bxmin=xmin+(xmax-xmin)/nxboxes*i
        bxmax=xmin+(xmax-xmin)/nxboxes*(i+1)
        bymin=ymin+(ymax-ymin)/nyboxes*j
        bymax=ymin+(ymax-ymin)/nyboxes*(j+1)
        avx=(bxmin+bxmax)/2.0
        avy=(bymin+bymax)/2.0

        thisbox=where(x gt bxmin and x le bxmax and y gt bymin and y le bymax)
        if thisbox[0] eq -1 then begin
;            openu, unit, 'fig5.dat', /get_lun
            printf, unit1, avx, avy, '      0', '     -0.1', '      0' 
;            free_lun, unit
        endif

        if thisbox[0] ne -1 then begin
             datacolor=median(data(thisbox))
             finddat1=where(data(thisbox) gt 0.0)
	     ndatgood=n_elements(finddat1)
             fraccolor=1.0*ndatgood/n_elements(data(thisbox))
	     if finddat1[0] eq -1 then fraccolor=-0.1


             printf, unit1, avx, avy, datacolor, fraccolor, n_elements(data(thisbox))
;	     wait, 100
;             free_lun, unit
;             standdev=stddev(data(thisbox))
             checka5=where(data(thisbox) gt 5.0)
             if checka5[0] ne -1 then above5=data(thisbox(where(data(thisbox) gt 5.0)))
             if checka5[0] ne -1 then above5err=dataerr(thisbox(where(data(thisbox) gt 5.0)))
;             order5=sort(above5)
;             order5err=sort(above5err)
;             if order5(0)-order5err gt 5.0 then begin
;             
;             bad=where(data(thisbox) gt 5.0)
;             endif else begin
             

;	     if n_elements(thisbox) gt 1 then standdev=medabsdev(data(thisbox))
;	     bad=where(data(thisbox) gt 5.0+2.0*standdev and data(thisbox) gt 5.0)
             bad=where(data(thisbox) gt 5.0)
;             if n_elements(where(data(thisbox) gt 5.0)) lt 3.0 then bad=where(data(thisbox) gt 5.0)

	     if bad[0] ne -1 then print, KICID(thisbox(bad)), data(thisbox(bad)), dataerr(thisbox(bad)), n_elements(thisbox)

             colorz=(datacolor-min(data))*250/(max(data)-min(data))
             polyfill, [bxmin, bxmin, bxmax, bxmax], $
                 [bymin, bymax, bymax, bymin], color=colorz 

;        if n_elements(thisbox) gt 1 and max(data(thisbox)) gt 5.0 then begin
;	set_plot, 'ps'
;	!p.font=0
;        bxminstr=strcompress(string(bxmin), /remove_all)
;        filname='~/Documents/Apogee/APGgraphs/'+'histx'+bxminstr+ $
;                       '_'+strcompress(string(bxmax), /remove_all)+'y' +$
;                       strcompress(string(bymin), /remove_all)+'_'+$
;                       strcompress(string(bymax), /remove_all)+'.ps'
;;;
;	device, file=filname, bits=8,$
;        	 /inches, /encapsulated, xsize=8, decomposed=0,$
;        	 ysize=5, /portrait,/color, set_font='Helvetica'
;       yhistmax=max([n_elements(data(thisbox))/3., n_elements(bad)] )
;         plot,  [min(data(thisbox)), max(data(thisbox))], [0, yhistmax], /nodata ;,   $
;            CLIP=[min(data(thisbox)), max(data(thisbox)),0, yhistmax], NOCLIP=0, $
;            XTITLE='data', $
;            YTITLE='Number', xstyle=1, $
;            ystyle=1, xthick=5, ythick=5, $
;            charsize=2.3, charthick=5
;           charsize=3, charthick=5
         ;, $
;        xrang=max(data(thisbox))-min(data(thisbox));
;	plothist, data(thisbox),xhist1, yhist, bin=xrang/20, /overplot
;        oplot, [5.0, 5.], [0, 10], color=255
;        if checka5[0] ne -1 then oploterror, above5, make_array(n_elements(above5), /float,/index), above5err, above5err*0.0, psym=1, color=150
;        xyouts,  data(thisbox(bad)), make_array(n_elements(above5), /float,/index), KICID(thisbox(bad))
;        device, /close
;        endif
        endif

   end
end


free_lun, unit1






end 
