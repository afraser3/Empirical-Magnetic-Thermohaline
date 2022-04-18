pro plotfig5bMixing


;bkg=sample

;xstring='A3PSeisMscale';'A3PSeisMscale' ;CorTeff16';'halfcorTeff16mbovytref' ;'ProtT12allbad';;'periodmaxH14v13';'SeisMscaleP12';'OmegaOmegacritP12' ;'SeisMprelim';
;ystring= 'CorFeH16';'CorFeH16' ; 'AldoResidual';'F8noFefractionaloffset';'SeisLoggscaleP12';'TeffAPGmMod18prelim'
;colorstring='C12C13b', 'CN16' ;

;xmin=0.8;-300;0.5;2.5;0.5;2.5;0.5;5000;0000000;3800;0.5 ;4400;3500 ;
;xmax=2.4;;5700;0;2.0; 99999999999999999;3.3;15;5200;30;3.0;0.5;5600;.0 ;2.5;1.3 ;4; 5000;5200;1.4 ;8000;
;ymin=-0.8;-1.0;1.0;-0.1;9.25;2.0;-0.2;-1;0;0.001;2.3;1.0;-1;-0.5;3;1.0;0000000;1.0 ;2 ;7 ;1.5 ;0. ;-3 ;2.2 ;
;ymax=0.6;3.1;4.0;1000000;4.0;0.4;12.61;5.0;0.2; 99999999999999999;.1;2.6;0.1;3.5;5.0;0.5;3.8;4.0  ;14 ;20 ;3.1;4.0 ;1.0 ;2.5 ;5.0;
;colmin= -9998;0.9  2.99;-0.1; -10000;
;colmax= 9999999999999999;1.1;2.5;  3.01; 0.1;-10;
;APGquick16savb, 'A3Cnotclumplogglt2gt1'
;APGquick16savb, 'A3Cnotclumplogggt3'   


c12c13='false'
;  mmin=0.8
;  mmax=2.4
;  loggmin=-0.8
;  loggmax=0.6
;  nelms=100
  mmin=0.8
  mmax=2.8
  loggmin=-1.5
  loggmax=0.5
  nelms=100


  delta=0.05
  delta2=0.09
  gs=10 ; grid size

!p.font=1
!p.thick=2
!x.thick=2
!y.thick=2
!z.thick=2
loadct, 0
  set_plot, 'ps'
  device, filename='Mfeh5c_mixingCRGBdr14v7.ps', /inches, /tt_font, set_font='Helvetica', xsize=8, ysize=8, font_size=10, /encapsulated 


;readcol, '~/Documents/Apogee/idl/fig5notclumplogg2-1C12C13.dat', ms, allloggs, kms1, junk, nstar1 ;logg3_b17
;readcol, '~/Documents/Apogee/idl/fig5notclumplogg2-1C12C13.dat', ms, allloggs, kms2, junk, nstar2

;readcol, '~/Documents/Apogee/idl/fig5notclumplogg3.dat', ms, allloggs, kms1, junk, nstar1
;readcol, '~/Documents/Apogee/idl/fig5notclumplogg2-1.dat', ms, allloggs, kms2, junk, nstar2

readcol, '~/Documents/Apogee/idl/fig5-CRGBlogggt28.dat', ms, allloggs, kms1, junk, nstar1
readcol, '~/Documents/Apogee/idl/fig5-CRGBlogglt23gt10.dat', ms, allloggs, kms2, junk, nstar2

;readcol, '~/Documents/Apogee/idl/fig5-A3Pdr14logg3.dat', ms, allloggs, kms1, junk, nstar1
;readcol, '~/Documents/Apogee/idl/fig5-A3Pdr14logg2-1.dat', ms, allloggs, kms2, junk, nstar2

;readcol, '~/Documents/Apogee/idl/fig5-EMPlogggt3.dat', ms, allloggs, kms1, junk, nstar1
;readcol, '~/Documents/Apogee/idl/fig5-EMP4logglt2gt1.dat', ms, allloggs, kms2, junk, nstar2

;readcol, '~/Documents/Apogee/idl/fig5-EMPdr14logggt3z1.dat', ms, allloggs, kms1, junk, nstar1
;readcol, '~/Documents/Apogee/idl/fig5-EMPdr14logglt2gt1z1.dat', ms, allloggs, kms2, junk, nstar2

;WARNING this does not check that the grids are the same!!!!

kms=kms1-kms2
if c12c13 eq 'true' then kms=kms2
bad=where(nstar1 lt 3 or nstar2 lt 3) ; 5 was standard
print, ms(56), allloggs(56), kms1(56), kms2(56), nstar1(56), nstar2(56)
kms[bad]=0.0
massmsun='Mass (M'+TexToIDL("_{\sun}")+' )'  ;TeXtoIDL("  sin   (km s^{-1})")
deltacnlabel=TexToIDL("\Delta [C/N]")
kms2=reform(kms, gs, gs)
ms2=reform(ms, gs, gs)
allloggs2=reform(allloggs, gs,gs)
;  plot,[ms[0], ms[1]], [kms[0],kms[1]], /nodata, /ynozero, xrange=[mmax, mmin],$
;    yrange=[loggmax, loggmin],$
;     ystyle=1,xstyle=1,$
;       xtitle=massmsun,$ 
;       ytitle='[Fe/H]', $
;       xthick=10, ythick=10, charsize=2.3, charthick=5  ;xthick=5, ythick=5, charsize=2.3, charthick=5  

xmin=mmin
xmax=2.0 ;mmax
ymin=-1.0 ;loggmin
ymax=loggmax

  plot,[ms[0], ms[1]], [kms[0],kms[1]], /nodata, /ynozero, xrange=[xmax, xmin],$
    yrange=[ymin,ymax],$
     ystyle=1,xstyle=1,$
       xtitle=massmsun,$ 
       ytitle='[Fe/H]', $
       xthick=10, ythick=10, charsize=2.7, charthick=5  ;xthick=5, ythick=5, charsize=2.3, charthick=5  
cgloadct, /brewer,18, /silent
if c12c13 eq 'true' then loadct, 0

;kms=-1.0*kms ;just to test if the unmixing is the problem delete later



goodcol=where( kms gt 0)
;badcol=where(kms eq 0)
minkms=min(kms(goodcol))
 for i=0, n_elements(kms)-1 do begin
       bxmin=ms[i]-delta
       bxmax=ms[i]+delta
       bymin=allloggs[i]-delta
       bymax=allloggs[i]+delta
       colorz=(kms[i]+max(kms))*255./(2.0*max(kms));255
       if kms[i] eq 0 then colorz=-9999

if c12c13 eq 'true' then  colorz=255-((kms[i]-minkms)*255./(max(kms)-minkms)) ;c12c13

;       if kms[i] lt 5.0 and kms[i] gt 0.0 then colorz=220
;       if kms[i] gt 5.0 and kms[i] lt 10.0 then colorz=160
;       if kms[i] gt 10.0 then colorz=99  
;       print, bymin, colorz
       if (bxmin gt xmin and bxmax lt xmax and bymin gt ymin and bymax lt ymax and colorz ne -9999) then polyfill, [bxmin, bxmin, bxmax, bxmax], $
                 [bymin, bymax, bymax, bymin], color=colorz 
 

 end

loadct, 0, /silent

 for i=0, n_elements(kms)-1 do begin
       bxmin=ms[i]-delta
       bxmax=ms[i]+delta
       bymin=allloggs[i]-delta
       bymax=allloggs[i]+delta
       bxmin2=ms[i]-delta2
       bxmax2=ms[i]+delta2
       bymin2=allloggs[i]-delta2
       bymax2=allloggs[i]+delta2
print, bymax2, bymin2
       ; outline polyfill region
       if i ne 0 and i ne (n_elements(kms)-1) then begin
           if( kms[i] ne 0 and kms[i-1] eq 0) then oplot, [bxmin2, bxmax2], [bymin2,bymin2], thick=2
           if( kms[i] ne 0 and kms[i+1] eq 0) then oplot, [bxmin2, bxmax2], [bymax2,bymax2], thick=2
   
       endif
       if i ge gs and i le n_elements(kms)-gs-1 then begin
           if( kms[i] ne 0 and kms[i-gs] eq 0) then oplot, [bxmin2, bxmin2], [bymin,bymax2], thick=2
           if( kms[i] ne 0 and kms[i+gs] eq 0) then oplot, [bxmax2, bxmax2], [bymin,bymax2], thick=2
       endif
       if i lt gs and i le n_elements(kms)-gs-1 then begin
           if( kms[i] ne 0 and kms[i+gs] eq 0) then oplot, [bxmax2, bxmax2], [bymin,bymax2], thick=2   
       endif
   end
xyouts, 1.7, -0.73, deltacnlabel, charsize=2.5

cgloadct, /brewer,18, /silent
;   axis, yaxis=1, /ynozero, xrange=[mmax, mmin],$
;    yrange=[loggmax, loggmin], ytickformat="(A1)",$
;     ystyle=1,xstyle=1,xthick=10, ythick=10, charsize=0, charthick=5 
;   axis, xaxis=0,  xrange=[mmax, mmin],xtickformat="(A1)",$
;    yrange=[loggmax, loggmin],$
;     ystyle=1,xstyle=1,xthick=10, ythick=10, charsize=2.3, charthick=5 
;   axis, yaxis=0, /ynozero, xrange=[mmax, mmin],$
;    yrange=[loggmax, loggmin], ytickformat="(A1)",$
;     ystyle=1,xstyle=1,xthick=10, ythick=10, charsize=0, charthick=5 
;   axis, xaxis=1,  xrange=[mmax, mmin],xtickformat="(A1)",$
;    yrange=[loggmax, loggmin],$
;     ystyle=1,xstyle=1,xthick=10, ythick=10, charsize=2.3, charthick=5 
;  contour, kms2, ms2, allloggs2, /overplot, $
;        /fill, c_color=[220,160,99,40], $
;         levels=[1.0, 5.0, 10.0]
;  contour,kms2, ms2, allloggs2 , /overplot , levels=[1.0,5.0,10.0];, c_labels=[1,1,1], c_charsize=0.1
rams=[0.86, 1.32, 1.66,1.17, 1.43, 1.25, 1.36, 1.56, 1.37, 0.88, 1.18, 1.53, -9999, 0.87, -9999, 1.03, -9999, -9999, -9999, 0.82, 2.07, 1.09, 0.87, 1.11, 0.91, 1.16]
raseislogg=[2.402, 2.443, 2.573, 2.362, 2.697, 2.425, 2.389, 2.486, 2.511, 2.337, 2.442, 2.493, -9999, 2.367, 3.2, 2.36, -9999, -9999, -9999, 2.374, 1.608, 2.379, 2.384, 2.375, 2.356, 2.399]
;oplot, rams, raseislogg, psym=sym(6),  symsize=2.0

;oplot, rams[14:25], raseislogg[14:25], psym=sym(5), color=20, symsize=1.5
;oplot, rams[14:25], raseislogg[14:25], psym=sym(10),  symsize=1.5
;oplot, rams[0:13], raseislogg[0:13], psym=sym(4), color=120, symsize=2.0


;oplot, rams[0:13], raseislogg[0:13], psym=sym(9),  symsize=2.0
;oploterror, [3.2, 3.2], [3.2, 3.2], [0.14, 0.14], [0.012, 0.012]


if c12c13 ne 'true' then COLORBAR,c_value=[255, 127, 0], /border, charsize=2.3,vertical=1,  nlevels=3, divisions=2, bottom=0, /invert,  $
         position=[0.35, 0.20, 0.40, 0.40],  tickname=[strmid(string(max(kms)), 4,5), $
                          strmid(string(.5*(max(kms)-max(kms))), 4, 6), $
			strmid(string(-1*max(kms)), 4, 5)] ;str(,'<5 '+TeXtoIDL("km s^{-1}"),'5-10 '+TeXtoIDL("km s^{-1}"),'>10 '+TeXtoIDL("km s^{-1}")]

if c12c13 eq 'true' then COLORBAR,c_value=[255, 127, 0], /border, vertical=1,  nlevels=3, divisions=2, bottom=0, $
         position=[0.45, 0.70, 0.50, 0.90],  tickname=[string(max(kms)), $
                          string(.5*(max(kms)-minkms)+minkms), $
			string(minkms)] ;str(,'<5 '+TeXtoIDL("km s^{-1}"),'5-10 '+TeXtoIDL("km s^{-1}"),'>10 '+TeXtoIDL("km s^{-1}")]
device, set_font='Helvetica Italic', /tt_font, font_size=10

xyouts, 3.83, 1.89, 'g', orientation=90, charsize=2.3, charthick=5
  device, /close
end

