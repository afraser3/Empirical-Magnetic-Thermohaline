pro APGquick16savb, startype

;Plot metallicity vs Mass for APOKASC data
;requires readcol2.pro, bkgcolor.pro, match2.pro, linfit.pro

restore, 'APGquick16varsb.sav'

; can be deleted next time saveAPGquick16 is run
ruidehmergercandidates=[2583386, 2972876, 3097360, 3216736, 3648674, 4350501, 4667909, 4820412, 5166068, 5180345, 5341131, 5385518, 5556743, 5639438, 5696938, 5820672, 5857618, 6024654, 6118479, 6143256, 6182668, 6200178, 6273090, 6437547, 6620586, 6633766, 6690139, 6864132, 7121674, 7376259, 7457184, 7499531, 7509923, 7778197, 8023523, 8055108, 8096716, 8108732, 8277879, 8301083, 8364786, 8517859, 8558329, 8708536, 8827808, 9227589, 9326896, 9329200, 9474201, 9512519, 9642805, 9784586, 9903598, 9907511, 9945389, 10132814, 10420335, 10645209, 10712314, 11186389, 11302371, 11304067, 11465942, 11555266, 11805876, 12254159, 12301741]

pgbestev=A3CEvstates
addin=where(A3CEvstates lt 0.5 and evstatepg gt 0.5)
pgbestev[addin]=evstatepg[addin]




ToffsetFecorrect=CorTeff16-TeffGS 
fix1=where(CorFeH16 gt 0)
ToffsetFecorrect[fix1]=ToffsetFecorrect[fix1]-140.
fix2=where(CorFeH16 gt -0.8 and CorFeH16 lt 0)
ToffsetFecorrect[fix2]=ToffsetFecorrect[fix2]-140-(140./0.8)*CorFeH16[fix2]



halfcorTeff14mbovytref=CorTeff14-(-382.5*CorFeH14+4607.)-(Logg14-2.5)/.0018

bx1=5150 ;5150.
bx2=4900 ;4850.
by1=2.65 ;2.6
by2=2.8 ;2.8
bx3=5100 ;5050.
by3=3.3 ;3.3
m1=(by2-by1)/(bx2-bx1)
m2=(by3-by2)/(bx3-bx2)
bx4=(m1*bx3+by1-by3-m2*bx1)/(m1-m2)
by4=m2*(bx4-bx1)+by1

A3PSeisMscalegood=A3PSeisMscale
A3PSeisLoggscalegood=A3PSeisLoggscale
bad=where(abs(A3PSeisLoggscale-CorLogg16) gt 0.5 or A3PSeisMscale gt 5 or A3PSeisMscale lt 0.4)
A3PSeisMscalegood[bad]=-9999
A3PSeisLoggscalegood[bad]=-9999


pgbestmass=masspg
addin=where(A3PSeisMscalegood gt 0 and pgbestmass lt 0)
pgbestmass[addin]=A3PSeisMscalegood[addin]

pgbestrot=protpg
addin=where(ProtT12 gt 0 and pgbestrot lt 0)
pgbestrot[addin]=ProtT12[addin]
addin=where(periodmaxH14v13 gt 0 and vsinimine13 gt 5)
pgbestrot[addin]=periodmaxH14v13[addin]

A3PSeisMscalegood2=A3PSeisMscale
A3PSeisLoggscalegood2=A3PSeisLoggscale
bad=where(abs(A3PSeisLoggscale-CorLogg16) gt 0.5 or A3PSeisMscale gt 3.5 or A3PSeisMscale lt 0.5)
A3PSeisMscalegood2[bad]=-9999
A3PSeisLoggscalegood2[bad]=-9999
gd=where(SeisMdw gt 0)
A3PSeisMscalegood2[gd]=SeisMdw[gd]
A3PSeisLoggscalegood2[gd]=SeisLoggdw[gd]

zerom=EMPSeisLoggscale
zero=where(EMPSeisLoggscale ne -9999)
zerom[zero]=0.0

xlist=EMPSeisMscale;dr14b.J_MAG_2M-dr14b.K_MAG_2M;NFe14;Age/1000;ProtT12allbad;Teff13;periodmaxH14v13;halfcorTeff13mbovytref;OCTSeisMscale;;SeisMS;RA12;KICID;CorTeff16-(-382.5*CorFeH16+4607.)-(Logg16-2.5)/.0018;
xlisterr=CorFeH14;vsinimine13err;SeisM14scale;EMPSeisLoggscale; ProtT12allbad;F8predict;Teff13;periodmaxH14v13err;;OCTSeisMscale;Teff13mbovytref;CN13;RA12;KICID ;SeisLoggerr;CorTeff12err;OmegaOmegacritP12 ; SeisMerr
ylist=CorFeH14;CN17b;dr14b.NUV_MAG_GALEX-dr14b.J_MAG_2M;CorAlphaFe16;CorTeff16-TeffGS;EMPSeisLoggscale;;CoralphaFe16;Logg14;F8logg-F8seislogg+0.078*CorFeH14-0.013;(F8-F8predictnoFe)/(F8predictnoFe); KICint;F8Logg;KICint;Logg13;ProtMc;delLoggOCT13;AlFe;;OCTSeisLoggscale;Teff13mbovytref;delLoggOCT13;;KICID;SeisLoggscaleP12;TeffAPGmMod;
ylisterr=CorTeff14;ProtT12;vscatter;EMPSeisLoggscale;ProtT12allbad;F8;KICint;F8Logg;OCTSeisLoggscale;halfcorTeff13mbovytref;AlFe;halfcorTeff13mbovytref;jthalfcorbovyoffsetbr;SeisLoggscaleP12;CorLogg12;SeisLoggscaleP12;
colorlist=CN14;;CorTeff16; CN14 ;/1000;CoralphaFe14;;Logg14;CorFeH14;KICint;ALOG10(periodmaxH14v13);KICint;noseisreason;CorLogg12;KICID;velocityT12P12;SeisMscaleP12;OmegaOmegacritP12
colorlisterr=CN14;EMPSeisMscale;ProtT12allbad;CorFeH14;KICint;ALOG10(periodmaxH14v13);KICint;noseisreason;CorLogg13;CN13;CorLogg12err;KICID;velocityT12P12;SeisMscaleP12;

;CorTeff14-(-382.5*CorFeH14+4607)-(CorLogg14-2.5)/.0018 ;

xstring='EMPSeisMscale';'A3PSeisMscale' ;CorTeff16';'halfcorTeff16mbovytref' ;'ProtT12allbad';;'periodmaxH14v13';'SeisMscaleP12';'OmegaOmegacritP12' ;'SeisMprelim';
ystring= 'CorFeH14z' ; 'AldoResidual';'F8noFefractionaloffset';'SeisLoggscaleP12';'TeffAPGmMod18prelim'
colorstring='CN14' ;;'EMPSeisMscale';'vsinimine13cor';'noseisreason';'velocityT12P12';'OmegaOmegacritP12';'SeisEvstates'; 

;print, 'wait'
;wait, 100

colortab=100 ;0;34 ; 100=colorblind rainbow

reversecolors='off' ; on, off, not set yet to work with goodcolors
goodcolors='on';'off'
bkg='sample'  ; options: on, off, hist, frac, mean, sample, make sure density is off
writename='off' ; options: on, off
reversex= 'off'
reversey= 'off'
ploterr= 'off'
writeevstate= 'off'
density='good'  ; options: on, good (colored boxes over red xs), goodr (reverse colors) goodrstretch(reverse colors+stretch color bar) ;make sure bkg is off
age='off'
xhist='off' ; on, off, log
rotkms='off'
fitline='off'
runmedian='off'
pltextradata='off' ;off, on =Mike Lum's cluster cn data, 'MH"=marc hons tess data
plttracks='off' ;off, on, rot
plttracksall='off'
plttracksal='off'
plttracksbf='off'
plttracksallm='on'
label='on' 
xplottype='logg' ;options=['Teff', 'logg', 'R', 'L', 'CN'] ; not sure if L and C/N are 100% working. also C/N solar zeropoint should probably be adapted for the chosen models
yplottype='CN' ;
cplottype='mass' ; options=['Teff', 'logg', 'R', 'mass', 'FeH', 'alpha']

readoccam='off'


xmin=0.8;-300;0.5;2.5;0.5;2.5;0.5;5000;0000000;3800;0.5 ;4400;3500 ;
xmax=2.8;;5700;0;2.0; 99999999999999999;3.3;15;5200;30;3.0;0.5;5600;.0 ;2.5;1.3 ;4; 5000;5200;1.4 ;8000;
ymin=-1.5;-1.0;1.0;-0.1;9.25;2.0;-0.2;-1;0;0.001;2.3;1.0;-1;-0.5;3;1.0;0000000;1.0 ;2 ;7 ;1.5 ;0. ;-3 ;2.2 ;
ymax=0.5;0.5;3.1;4.0;1000000;4.0;0.4;12.61;5.0;0.2; 99999999999999999;.1;2.6;0.1;3.5;5.0;0.5;3.8;4.0  ;14 ;20 ;3.1;4.0 ;1.0 ;2.5 ;5.0;
colmin= -9998;0.5;0.9  2.99;-0.1; -10000;
colmax= 9999999999999999;1.1;2.5;  3.01; 0.1;-10;

;;;;;;;;END CHANGE THINGS

;writecol, 'T12outall.txt', KICID, TMASSID, P12Teff, P12Tefferr, CorFeH12, CorFeH12err, SeisMscaleP12, SeisLoggscaleP12, ProtT12, ProtT12err

outdirectory= '~/Documents/Apogee/APGgraphs/'
outfile=startype+xstring+ystring+colorstring+'_16.ps'
if density ne 'off' then outfile=startype+xstring+ystring+colorstring+'_16den.ps'
outtext=startype
if reversex eq 'on' then begin

      x1=xmax
      x2=xmin
endif else begin
      x1=xmin
      x2=xmax
endelse      
if reversey eq 'on' then begin
      y1=ymax
      y2=ymin
endif else begin
      y1=ymin
      y2=ymax
endelse    

; SET UP PLOT
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
loadct, 0
   yrang=(y2-y1)
   xrang=(x2-x1)
   set_plot, 'ps'
   !p.font=0
;   device, file=outdirectory+outfile, bits=8,$
;         /inches, /encapsulated, xsize=8, decomposed=0,$
;         ysize=5, /portrait,/color, set_font='Helvetica'
   device, file=outdirectory+outfile, $
         /inches, /encapsulated,/tt_font, xsize=8, decomposed=0,$
         ysize=8, font_size=10,/color, set_font='Helvetica'
  
   
    plot,  [x1, x2], [y1,y2], /nodata,   $
         CLIP=[x1,x2,y1,y2], NOCLIP=0, $
         XTITLE=xstring, $
         YTITLE=ystring, xstyle=1, $
    
       ystyle=1,$
       xthick=10, ythick=10, charsize=2.3, charthick=5  ,$;xthick=5, ythick=5, charsize=2.3, charthick=5 

;         ystyle=1, xthick=5, ythick=5, $
;         charsize=2.3, charthick=5, $
;         charsize=3, charthick=5, $
;         xticks=2, yticks=3, $
;         yminor=10,$ ;7 
         xrange=[x1, x2], yrange=[y1, y2] 
         ;, $

    if label ne 'off' then xyouts,  x2-(.2*xrang), y2-(.1*yrang), outtext, charsize=1, $
         charthick=2

;oplot, [xmin,xmax], [ymin, ymax]
;oplot, [0,0], [ymin, ymax]
;oplot, [xmin,xmax], [xmin, xmax]
;oplot, [0, 1.0], [2.76, 10.36+2.76]

;oplot, [xmin, xmax], [0, 0], linestyle=1
;oplot, [xmin, xmax], [-35, -35], linestyle=1
;oplot, [xmin, xmax], [35, 35], linestyle=1
;oplot, [xmin, xmax], [230, -175]
;oplot, [xmin, xmax], [-100, -100], thick=4
;oplot, [xmin, xmax], [-66, -66], thick=4
;for flicker
;oplot, [xmin, xmax], [-0.2*xmin+0.78, -0.2*xmax+0.78]

;for todd
;oplot,[xmin, xmax], [2.59, 2.59]
;oplot,[xmin, xmax], [2.35, 2.35]
;oplot, [4480, 4480],[ymin, ymax] 

;for sam
;oplot, [xmin, xmax], [2.88, 2.88]
;oplot, [4634, 4634],  [ymin, ymax]

;for todd2
;oplot,[xmin, xmax], [2.59, 2.59]
;oplot,[xmin, xmax], [2.35, 2.35]
;oplot, [4480, 4480],[ymin, ymax] 
;oplot, [4457, 4457], [9.25, 12.61]
;oplot, [4782, 4782], [9.25, 12.61]
;oplot, [xmin,xmax], [-0.336, -0.336]
;oplot, [xmin,xmax], [3.22, 3.22]
;oplot, [10.93, 10.93],[ymin, ymax]
;oplot, [-0.336, -0.336], [ymin,ymax]
;oplot, [xmin, xmax], [4457.833-(382.5*0.251+4607.)-(3.22-2.5)/.0018 , 4457.833-(382.5*0.251+4607.)-(3.22-2.5)/.0018 ]
;oplot, [xmin, xmax], [4782.447-(382.5*0.251+4607.)-(3.22-2.5)/.0018 , 4782.447-(382.5*0.251+4607.)-(3.22-2.5)/.0018 ]
if colortab lt 50 then loadct, colortab
if colortab eq 100 then cgloadct, /brewer,31, /silent
;for corrected gravities
;oplot, [xmin, -0.5], [140,100], color=250
;oplot, [-0.5, xmax], [100, -200], color=250

;for cor Teff, Fe/H, uncor Logg
;oplot, [xmin, -0.35], [-40,-60], color=100, thick=4
;oplot, [-0.35, xmax], [-60, -280], color=100, thick=4

;for bump
;oplot, [4650,4470], [2.6, 2.75]
;oplot, [4550,4420],  [2.48,2.58] 

;for cor Teff, Fe/H, uncor Logg
a=[xmin, xmax]

;b=[-0.35, xmax]
;only rc to rgb
;oplot, a, 405.28*a+261.74, color=250, thick=4

;3 groups at same time

;oplot, a, 405.28*a+261.74, color=250, thick=4
;oplot, a, -183.08*a-139.13, color=250, thick=4
;oplot, a, 62.66*a+15.22, color=250, thick=4
;break at -0.35
;a=[xmin,-0.35]
;oplot, a, -62.00*a-86.27, color=0, thick=4
;oplot, b, -248.47*b-151.97, color=0, thick=4
;break at -0.40
;oplot, a, -53.82*a-80.97, color=0, thick=4
;oplot, b, -233.50*b-149.40, color=0, thick=4

;break at -0.30
;oplot, a, -77.58*a-95.44, color=0, thick=4
;oplot, b, -280.28*b-155.33, color=0, thick=4

;Todds star
;oplot, [xmin, xmax], [2.5913, 2.5913], thick=4
;oplot, [4480, 4480], [ymin, ymax], thick=4

;oplot, [xmin, xmax], [-254.74, -254.74], thick=4
;oplot, [0.0037, 0.0037], [ymin, ymax], thick=4

loadct,0
;     oplot, [xmin, xmax], [5.0,5.0]
;     oplot, [xmin, xmax], [xmin,xmax]
;      oplot, [-1.5, 1.0], [-300, 400]
;oplot, [1.829, 1.835, 1.8435, 1.8549, 1.8549, 1.9035, 1.9126, 1.9137, 1.9144, 1.9162, 1.9396, 1.9715, 1.9995, 2.0283, 2.0697, 2.1058, 2.1371, 2.183, 2.227, 2.2661, 2.307, 2.3350, 2.3729, 2.3934, 2.3923, 2.3843], [-2.192, -2.068, -1.893, -1.768, -1.589, -1.492, -1.346, -1.288, -1.190, -1.013, -0.910, -0.808, -.710, -0.572, -0.476, -0.397, -0.270, -0.17, -0.089, 0.012, 0.096, 0.211, 0.305, 0.401, 0.469], linestyle=1
;PLOT THINGS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if colortab lt 50 then loadct, colortab
if colortab eq 100 then cgloadct, /brewer,31, /silent
; Taking the S2 values preferentially, but if there's no S2 but there is an S1, use that.

;    good1=where(xlist gt xmin and xlist lt xmax and ylist gt ymin and ylist lt ymax and xlisterr lt 9998 and $ 
;                xlisterr gt -9998 and ylisterr lt 9998 and ylisterr gt -9998 and colorlist gt -9998 and colorlist lt 9998 )

    good1=where(xlist gt xmin and xlist le xmax and ylist gt ymin and ylist le ymax and $ 
                xlisterr gt -9998 and ylisterr gt -9998 and colorlist gt colmin and colorlist le colmax)
print, n_elements(good1)

    xlist=xlist(good1)

print, min(xlist)
    xlisterr=xlisterr(good1)
    ylist=ylist(good1)
    ylisterr=ylisterr(good1)
    colorlist=colorlist(good1)
    colorlisterr=colorlisterr(good1)
    SpecVsini=SpecVsini(good1)
    CorLogg=CorLogg(good1)
    CorLoggerr=CorLoggerr(good1)
    KICID=KICID(good1)
    Evstates=Evstates(good1)
    TMASSID=TMASSID(good1)
    CorFeH=CorFeH(good1)
    CorFeHerr=CorFeHerr(good1)
    CoralphaFe=CoralphaFe(good1)
    CoralphaFeerr=CoralphaFeerr(good1)
    APGLocID=APGLocID(good1)
    CorTeff=CorTeff(good1)
    CorTefferr=CorTefferr(good1)
    ASPCAPchi2=ASPCAPchi2(good1)
    SeisM=SeisM(good1)
    SeisMerr=SeisMerr(good1)
    SeisR=SeisR(good1)
    SeisRerr=SeisRerr(good1)
    SeisMscale=SeisMscale(good1)
    SeisLogg=SeisLogg(good1)
    SeisLoggerr=SeisLoggerr(good1)

    SeisMS=SeisMS(good1)
    SeisM10=SeisM10(good1)
    SeisM10err=SeisM10err(good1)
    SeisR10=SeisR10(good1)
    SeisR10err=SeisR10err(good1)
    SeisLogg10=SeisLogg10(good1)
    SeisLogg10err=SeisLogg10err(good1)

    allEvstates=allEvstates(good1)
    SeisEvstates=SeisEvstates(good1)
    Evstates=Evstates(good1)
    vsinimine=vsinimine(good1)
    vsinimineerr=vsinimineerr(good1)
    vsini3d=vsini3d(good1)
    Evstate=Evstate(good1)

    Bevstate=Bevstate(good1)
    CoralphaFe=CoralphaFe(good1)
    Agemod=Agemod(good1)
    ProtMc=ProtMc(good1)
    ProtMcerr=ProtMcerr(good1)
    CorVmicro=CorVmicro(good1)
    TeffAPGmMod=TeffAPGmMod(good1)
    CN=CN(good1)
    CM=CM(good1)
    CH=CH(good1)
    FeH2=FeH2(good1)
    NH=NH(good1)
    OH=OH(good1)
    CHerr=CHerr(good1)
    FeH2err=FeH2err(good1)
    NHerr=NHerr(good1)
    OHerr=OHerr(good1)
    Teff12=Teff12(good1)
    CorTeff12=CorTeff12(good1)
    CorTeff12err=CorTeff12err(good1)
    Logg12=Logg12(good1)
    Logg12err=Logg12(good1)
    CorLogg12=CorLogg12(good1)
    CorLogg12err=CorLogg12err(good1)
    CorFeH12=CorFeH12(good1)
    FeH12=FeH12(good1)
    CorFeH12err=CorFeH12err(good1)
    CoralphaFe12=CoralphaFe12(good1)
    CoralphaFe12err=CoralphaFe12err(good1)
    Vhelio=Vhelio(good1)
    Vscatt=Vscatt(good1)
    vsinimine12=vsinimine12(good1)
    vsinimine12err=vsinimine12err(good1)
    ProtT12=ProtT12(good1)
    ProtT12err=ProtT12err(good1)
    velocityT12=velocityT12(good1)
    velocityT14=velocityT14(good1)
    ProtC=ProtC(good1)
    ProtCerr=ProtCerr(good1)
    periodmax=periodmax(good1)
    periodmaxerr=periodmaxerr(good1)
    inclination=inclination(good1)
    inclinationerr=inclinationerr(good1)

    periodmaxP12=periodmaxP12(good1)
    periodmaxP12err=periodmaxP12err(good1)
    inclinationP12=inclinationP12(good1)
    inclinationP12err=inclinationP12err(good1)
    velocityT12P12=velocityT12P12(good1)
    periodmaxP12v10=periodmaxP12v10(good1)
    periodmaxP12v10err=periodmaxP12v10err(good1)
    inclinationP12v10=inclinationP12v10(good1)
    inclinationP12v10err=inclinationP12v10err(good1)
    Bevstate2=Bevstate2(good1)
    massB=massB(good1)

    numaxRG=numaxRG(good1)
    delnuRG=delnuRG(good1)
KICgMag=KICgMag(good1)
RA12=RA12(good1)
DEC12=DEC12(good1)
PropRAMot12=PropRAMot12(good1)
PropDecMot12=PropDecMot12(good1)
KICEBmV12=KICEBmV12(good1)
vsiniN=vsiniN(good1)
numaxT12=numaxT12(good1)
numaxT12err= numaxT12err(good1)
delnuT12= delnuT12(good1)
delnuT12err= delnuT12err(good1)
    SeisMscaleP12=SeisMscaleP12(good1)
    SeisMscaleP12err=SeisMscaleP12err(good1)
    SeisRscaleP12=SeisRscaleP12(good1)
    SeisRscaleP12err=SeisRscaleP12err(good1)
    SeisLoggscaleP12=SeisLoggscaleP12(good1)
    SeisLoggscaleP12err=SeisLoggscaleP12err(good1)

    SeisM14scale=SeisM14scale(good1)
    SeisM14err=SeisM14err(good1)
    SeisR14scale=SeisR14scale(good1)
    SeisR14err=SeisR14err(good1)
    SeisLogg14scale=SeisLogg14scale(good1)
    SeisLogg14err=SeisLogg14err(good1)

P12Teff=P12Teff(good1)
P12Tefferr=P12Tefferr(good1)
CEvstates=CEvstates(good1)
;OCTSeisMscale=OCTSeisMscale(good1)
;OCTSeisLoggscale=OCTSeisLoggscale(good1)
;OCTSeisRscale=OCTSeisRscale(good1)

;Teff13=Teff13(good1)
;CorTeff13=CorTeff13(good1)
;CorTeff13err=CorTeff13err(good1)

;CorLogg13=CorLogg13(good1)
;Logg13=Logg13(good1)
;CorLogg13err=CorLogg13err(good1)

;FeH13=FeH13(good1)
;CorFeH13=CorFeH13(good1)
;CorFeH13err=CorFeH13err(good1)

;alphaFe13=alphaFe13(good1)
;CoralphaFe13=CoralphaFe13(good1)
;CoralphaFe13err=CoralphaFe13err(good1)

;CN13=CN13(good1)
;CN13err=CN13err(good1)




;OFe13=OFe13(good1)
APGvsini16=APGvsini16(good1)
SeisDataflag=SeisDataflag(good1)
l1sup=l1sup(good1)
;svmevstate=svmevstate(good1)
ConEvstates=ConEvstates(good1)
MassD=MassD(good1)

;SYDSeisMscale=SYDSeisMscale(good1)
;SYDSeisLoggscale=SYDSeisLoggscale(good1)
;SYDSeisRscale=SYDSeisRscale(good1)

;SYDSeisMscalecor=SYDSeisMscalecor(good1)
;SYDSeisLoggscalecor=SYDSeisLoggscalecor(good1)
;SYDSeisRscalecor=SYDSeisRscalecor(good1)

;OCTSeisMscalecor=OCTSeisMscalecor(good1)
;OCTSeisLoggscalecor=OCTSeisLoggscalecor(good1)
;OCTSeisRscalecor=OCTSeisRscalecor(good1)

;OSSeisMscaleerr=OSSeisMscaleerr(good1)
;OSSeisRscaleerr=OSSeisRscaleerr(good1)
;OSSeisLoggscaleerr=OSSeisLoggscaleerr(good1)

CorTeff14=CorTeff14(good1)
CorLogg14=CorLogg14(good1)
CorFeH14=CorFeH14(good1)
CoralphaFe14=CoralphaFe14(good1)
Teff14=Teff14(good1)
Logg14=Logg14(good1)
FeH14=FeH14(good1)
CorTeff14err=CorTeff14err(good1)
CorLogg14err=CorLogg14err(good1)
CorFeH14err=CorFeH14err(good1)
CoralphaFe14err=CoralphaFe14err(good1)
CN14=CN14(good1)
CN14err=CN14err(good1)

CorTeff16=CorTeff16(good1)
CorLogg16=CorLogg16(good1)
CorFeH16=CorFeH16(good1)
CoralphaFe16=CoralphaFe16(good1)
Teff16=Teff16(good1)
Logg16=Logg16(good1)
FeH16=FeH16(good1)
actualFeH16=actualFeH16(good1)
MgFe16=MgFe16(good1)
CorTeff16err=CorTeff16err(good1)
CorLogg16err=CorLogg16err(good1)
CorFeH16err=CorFeH16err(good1)
CoralphaFe16err=CoralphaFe16err(good1)
CN16=CN16(good1)
CN16err=CN16err(good1)
APGvsini16=APGvsini16(good1)

;OCTdelnu=OCTdelnu(good1)
;OCTnumax=OCTnumax(good1)

GHBteff=GHBteff(good1)
;OCTSeisMscaleGHB=OCTSeisMscaleGHB(good1)
;OCTSeisLoggscaleGHB=OCTSeisLoggscaleGHB(good1)
;OCTSeisRscaleGHB=OCTSeisRscaleGHB(good1)

nvisits=nvisits(good1)
vscatter=vscatter(good1)
KICint=KICint(good1)
F8Logg=F8Logg(good1)
N_KEP_QUART=dr14b.N_KEP_QUART
N_KEP_QUART=N_KEP_QUART(good1)
TARGFLAGS=dr14b.TARGFLAGS
TARGFLAGS=TARGFLAGS(good1)
ASPCAPFLAGS=dr14b.ASPCAPFLAGS
ASPCAPFLAGS=ASPCAPFLAGS(good1)
H14R=H14R(good1)
H14M=H14M(good1)
H14logg=H14logg(good1)

periodmaxH14v13=periodmaxH14v13(good1)
periodmaxH14v13err=periodmaxH14v13err(good1)
inclinationH14v13=inclinationH14v13(good1)
inclinationH14v13err=inclinationH14v13err(good1)
periodmaxP12v13=periodmaxP12v13(good1)
periodmaxP12v13err=periodmaxP12v13err(good1)
inclinationP12v13=inclinationP12v13(good1)
inclinationP12v13err=inclinationP12v13err(good1)
vsinimine13cor=vsinimine13cor(good1)
vsinimine13err=vsinimine13err(good1)
noseisreason=noseisreason(good1)
periodmax13=periodmax13(good1)
periodmaxSDSSv13=periodmaxSDSSv13(good1)
periodmaxSDSSv13err=periodmaxSDSSv13err(good1)
prott12allbad=prott12allbad(good1)

   SeisMdw=SeisMdw(good1)
   SeisMdwerr=SeisMdwerr(good1)
   SeisRdw=SeisRdw(good1)
   SeisRdwerr=SeisRdwerr(good1)
   SeisLoggdw=SeisLoggdw(good1)
   SeisLoggdwerr=SeisLoggdwerr(good1)

F8seislogg=F8seislogg(good1)
CannonLogg=CannonLogg(good1)

;F8seismass=F8seismass(good1)
;F8seismasserr=F8seismasserr(good1)
;F8numax=F8numax(good1)
;F8numaxerr=F8numaxerr(good1)
;F8teff=F8teff(good1) 
;F8tefferr=F8tefferr(good1)
;F8feh=F8feh(good1)
;F8feherr=F8feherr(good1)
;F8=F8(good1)
;F8err=F8err(good1)
;F8Evstates=F8Evstates(good1)

F8teff14=F8teff14(good1) 
F8teff14err=F8teff14err(good1)
F8feh14=F8feh14(good1)
F8feh14err=F8feh14err(good1)
F8logg14=F8logg14(good1)
F8logg14err=F8logg14err(good1)

EMPSeisMscale=EMPSeisMscale(good1)
EMPSeisMscaleerr=EMPSeisMscaleerr(good1)
EMPSeisRscale=EMPSeisRscale(good1)
EMPSeisRscaleerr=EMPSeisRscaleerr(good1)
EMPSeisLoggscale=EMPSeisLoggscale(good1)
EMPSeisLoggscaleerr=EMPSeisLoggscaleerr(good1)
EMPSeisAge=EMPSeisAge(good1)
EMPSeisAgeerr=EMPSeisAgeerr(good1)

Gaiadistance=Gaiadistance(good1)


A3PSeisLoggscale=A3PSeisLoggscale(good1)
A3PSeisLoggscaleerr=A3PSeisLoggscaleerr(good1)
A3PSeisMscale=A3PSeisMscale(good1)
A3PSeisMscaleerr=A3PSeisMscaleerr(good1)
A3PSeisRscale=A3PSeisRscale(good1)
A3PSeisRscaleerr=A3PSeisRscaleerr(good1)

A3CEvstates=A3CEvstates(good1)
PWbinaryP=PWbinaryP(good1)

loggpg=loggpg(good1)
masspg=masspg(good1)
Protpg=Protpg(good1)
vexpectedpg=vexpectedpg(good1)
Rpg=Rpg(good1)
pgbestev=pgbestev(good1)
pgbestrot=pgbestrot(good1)

RA=RA(good1)
DEC=DEC(good1)
Kepmag=Kepmag(good1)
A3Pnumax=A3Pnumax(good1)
AgeGS=AgeGS(good1)

corerotcg=corerotcg(good1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;   print, 'M10', KICID[0:2], SeisM10[0:2], SeisM10err[0:2]



;writecol,'tugdualprots.txt', KICID, xlist, ylist

print, 'mean metallicity14', avg(CorFeH14)
print, 'mean Teff14', avg(CorTeff14)
;change this back when done
    periodmax3d=2*3.14159*SeisR*Rsun/vsini3d/60./60./24. ;days

;    wait, 100
 
;    if bkg eq 'on' then bkgcolor, xmin, xmax, ymin, ymax, 15.0, 15.0, xlist, ylist, colorlist, KICID=KICID
    if bkg eq 'on' then bkgcolor2, xmin, xmax, ymin, ymax,  30.0, 30.0, xlist, ylist, colorlist, colorlisterr, KICID=KICID
    if bkg eq 'mean' or bkg eq 'sample' then begin
          openw, unit1, 'fig5.dat', /get_lun
          printf, unit1, xstring, '  ', ystring, '  ', colorstring
          free_lun, unit1
    endif
;wait, 100
    if bkg eq 'mean' then bkgcolor2mean, xmin, xmax, ymin, ymax,  10.0, 10.0, xlist, ylist, colorlist, colorlisterr, KICID=KICID
    if bkg eq 'hist' then bkgcolorhistpub, xmin, xmax, ymin, ymax, 30.0, 30.0, xlist, ylist, colorlist, colorlisterr, KICID=KICID
    if density eq 'on' then densitycolor, xmin, xmax, ymin, ymax, 20.0, 20.0, xlist, ylist, colorlist, KICID=KICID, label

    if startype eq 'modelalphatestm11fhm06' then good=where(OCTSeisMscale gt 1.0 and OCTSeisMscale lt 1.2 and CorFeH13 gt -0.7 and CorFeH13 lt -0.5 and CEvstates eq 3.0)
    if startype eq 'modelalphatestm10fhm04' then good=where(OCTSeisMscale gt 0.9 and OCTSeisMscale lt 1.1 and CorFeH13 gt -0.5 and CorFeH13 lt -0.3 and CEvstates eq 3.0)
    if startype eq 'modelalphaemptestm10fhm06' then good=where(EMPSeisMscale gt 0.8 and EMPSeisMscale lt 1.2 and CorFeH14 gt -0.7 and CorFeH14 lt -0.5 and CEvstates eq 3.0)
    if startype eq 'modelalphaemptestm11fhm04' then good=where(EMPSeisMscale gt 1.0 and EMPSeisMscale lt 1.2 and CorFeH14 gt -0.5 and CorFeH14 lt -0.3 and CEvstates eq 3.0)
    if startype eq 'RGBsolFe' then good=where(CEvstates eq 3.0 and CorFeH13 gt -0.1 and CorFeH13 lt 0.1) 
    if startype eq 'svm2CL' then good=where(svmevstate eq '2CL')
    if startype eq 'l1norm' then good=where(l1sup eq 0)
    if startype eq 'l1sup' then good=where(l1sup eq 1)
    if startype eq 'MD1_58' then good=where(MassD ge 1.5 and MassD le 1.8)
    if startype eq 'MD1_7' then good=where(MassD ge 1.6 and MassD le 1.8)
    if startype eq 'MD2_0' then good=where(MassD ge 1.9 and MassD le 2.1)
    if startype eq 'SeisDataflag0' then good=where(SeisDataflag eq 0)
    if startype eq 'SeisDataflag9' then good=where(SeisDataflag eq 9)
    if startype eq 'OCTSeisLoggscalegt25' then good=where(OCTSeisLoggscale gt 2.5)
    if startype eq 'OCTSeisMscalelt3' then good=where(OCTSeisMscale lt 3.0)
    if startype eq 'OCTLogg25_27' then good=where(OCTSeisLoggscale lt 2.7 and OCTSeisLoggscale lt 2.5)
    if startype eq 'OCTSeisMscale1_7' then good=where(OCTSeisMscale lt 1.8 and OCTSeisMscale gt 1.6)
    if startype eq 'OCTSeisMscale1_3' then good=where(OCTSeisMscale lt 1.4 and OCTSeisMscale gt 1.2)
    if startype eq 'CRGB' then good=where(CEvstates eq 3.0)
    if startype eq 'CRGBvsinigt2' then good=where(CEvstates eq 3.0 and vsinimine12 gt 2.0)
    if startype eq 'Cclump' then good=where(CEvstates eq 1.0)
    if startype eq 'CRGBAGB' then good=where(CEvstates eq 3.5)
    if startype eq 'C2clump' then good=where(CEvstates eq 2.0)
    if startype eq 'solar2_3' then good=where(CorFeH13 gt -0.1 and CorFeH13 lt 0.1 and EMPSeisMscale ge 2.0 and EMPSeisMscale le 3.0)
    if startype eq 'C12clump' then good=where(CEvstates eq 1.5)

    if startype eq 'Callclump' then good=where(CEvstates ge 1.0 and CEvstates le 2.0) 
    if startype eq 'closeCallclump' then good=where(CEvstates ge 1.0 and CEvstates le 2.0 and Gaiadistance lt 1.0) 
    if startype eq 'closeCRGB' then good=where(CEvstates eq 3.0 and Gaiadistance lt 1.0) 
    if startype eq 'farCallclump' then good=where(CEvstates ge 1.0 and CEvstates le 2.0 and Gaiadistance ge 1.0) 
    if startype eq 'farCRGB' then good=where(CEvstates eq 3.0 and Gaiadistance ge 1.0)
    if startype eq 'CUReject' then good=where(CEvstates ge 4.5 and CEvstates le 6.0)

    if startype eq 'lowzrgb' then good=where(CorFeH13 gt -1.0 and CorFeH13 lt -0.5 and CEvstates eq 3)
    if startype eq 'highzrgb' then good=where(CorFeH13 gt 0.0 and CorFeH13 lt 0.5 and CEvstates eq 3)
    if startype eq 'gapcross' then good=where(Bevstate2 eq 'R' and massB gt 2.0)
    if startype eq 'logggt4'  then good=where(CorLogg12 gt 4.0)
    if startype eq 'logggt35'  then good=where(CorLogg13 gt 3.5)
    if startype eq 'all' then good=where(xlist ne -9999)
    if startype eq 'vsini' then good=where(SpecVsini gt 1.0 and SpecVsini ne 9999)
    if startype eq 'fast' then good=where(SpecVsini gt 10.0 and SpecVsini ne 9999) 
    if startype eq 'vsini5' then good=where(SpecVsini gt 5.0 and SpecVsini ne 9999)   
    if startype eq 'logglt3' then good=where(CorLogg13 lt 3.0)
    if startype eq 'logglt3fehgt0' then good=where(CorLogg13 lt 3.0 and CorFeH13 gt 0.0)
    if startype eq 'logglt3fehm05_0' then good=where(CorLogg13 lt 3.0 and CorFeH13 gt -0.5 and CorFeH13 lt 0.0)
    if startype eq 'logggt2' then good=where(SeisLogg gt 2.0)
    if startype eq 'logglt2' then good=where(SeisLogg lt 2.0)
    if startype eq 'specloggclump' then good=where(CorLogg gt 2.4 and CorLogg lt 2.9)
    if startype eq 'loggclump' then good=where(SeisLogg gt 2.4 and SeisLogg lt 2.9)
    if startype eq 'loggclumpnotRGB' then good=where(SeisLogg gt 2.3 and SeisLogg lt 2.7 and SeisR gt 8. and SeisR lt 12. and Evstates ne 3)
    if startype eq 'loggclumpnoEv' then good=where(SeisLogg gt 2.3 and SeisLogg lt 2.7 and Evstates ne 2 and Evstates ne 1 and Evstates ne 3)
    if startype eq 'dwarfsubgiant' then good=where(Evstates eq 1)
    if startype eq 'clump' then good= where(SeisEvstates gt 1.5 and SeisEvstates lt 2.5)
    if startype eq 'RGB' then good= where(SeisEvstates gt 2.5 and SeisEvstates lt 3.5)
    if startype eq 'allclump' then good= where(Evstates eq 2 or Bevstate eq 'clump' or Bevstate eq 'clump?' or (allEvstates gt 1.5 and allEvstates lt 2.5))
    if startype eq 'allRGB' then good= where(Evstates eq 3 or Bevstate eq 'RGB' or Bevstate eq 'RGB?' or (allEvstates gt 2.5 and allEvstates lt 3.5) )
    if startype eq 'lowFeH' then good=where(CorFeH13 lt -1.5)
    if startype eq 'FeHm35' then good=where(CorFeH16 lt -0.3 and CorFeH16 gt -0.4)
    if startype eq 'midlowFeH' then good=where(CorFeH13 lt -0.5 and CorFeH13 gt -1.3)
    if startype eq 'highFeH' then good=where(CorFeH13 gt 0.0)
    if startype eq 'highalpha14' then good=where(CoralphaFe14 gt 0.15)
    if startype eq 'highalpha14Cclump' then good=where(CoralphaFe14 gt 0.15 and CEvstates eq 1.0)
   if startype eq 'solar14' then good=where(CorFeH14 gt -.1 and  CorFeH14 lt 0.1 )
   if startype eq 'solar14Callclump' then good=where(CorFeH14 gt -.1 and  CorFeH14 lt 0.1 and CEvstates ge 1.0 and CEvstates le 2.0)
   if startype eq 'solar14CRGB' then good=where(CorFeH14 gt -.1 and  CorFeH14 lt 0.1 and CEvstates eq 3.0 )
   if startype eq 'CN14Callclump' then good=where(CorFeH14 gt -.5 and  EMPSeisLoggscale lt 3.1 and EMPSeisLoggscale gt 2.3 and EMPSeisMscale lt 1.7 and CEvstates ge 1.0 and CEvstates le 2.0 )
   if startype eq 'CN14CRGB' then good=where(CorFeH14 gt -.5 and  EMPSeisLoggscale lt 3.1 and EMPSeisLoggscale gt 2.3 and EMPSeisMscale lt 1.7 and CEvstates eq 3.0 )
   if startype eq 'solar14CRGB' then good=where(CorFeH14 gt -.1 and  CorFeH14 lt 0.1 and CEvstates eq 3.0 )
   if startype eq 'solar14CRGBm09' then good=where(CorFeH14 gt -.1 and  CorFeH14 lt 0.1 and CEvstates eq 3.0 and EMPSeisMscale gt 0.8 and EMPSeisMscale lt 1.0)
   if startype eq 'solar14CRGBm13' then good=where(CorFeH14 gt -.1 and  CorFeH14 lt 0.1 and CEvstates eq 3.0 and EMPSeisMscale gt 1.2 and EMPSeisMscale lt 1.4)
   if startype eq 'subsolar14Callclump' then good=where(CorFeH14 gt -.3 and  CorFeH14 lt -0.1 and CEvstates ge 1.0 and CEvstates le 2.0)
   if startype eq 'subsolar14CRGB' then good=where(CorFeH14 gt -.3 and  CorFeH14 lt -0.1 and CEvstates eq 3.0 )
   if startype eq 'subsolar14CRGBm09' then good=where(CorFeH14 gt -.3 and  CorFeH14 lt -0.1 and CEvstates eq 3.0 and EMPSeisMscale gt 0.8 and EMPSeisMscale lt 1.0)
   if startype eq 'subsolar14CRGBm13' then good=where(CorFeH14 gt -.3 and  CorFeH14 lt -0.1 and CEvstates eq 3.0 and EMPSeisMscale gt 1.2 and EMPSeisMscale lt 1.4)
   if startype eq 'supersolar14Callclump' then good=where(CorFeH14 gt 0.1 and  CorFeH14 lt 0.3 and CEvstates ge 1.0 and CEvstates le 2.0)
   if startype eq 'supersolar14CRGB' then good=where(CorFeH14 gt 0.1 and  CorFeH14 lt 0.3 and CEvstates eq 3.0 )
   if startype eq 'highalpha14solar' then good=where(CoralphaFe14 gt 0.15 and  CorFeH14 gt -0.2 and CorFeH14 lt 0.0)
    if startype eq 'highalpha14low' then good=where(CoralphaFe14 gt 0.15 and  CorFeH14 gt -1.0 and CorFeH14 lt -0.8)

    if startype eq 'highalpha14CRGBsolar' then good=where(CoralphaFe14 gt 0.15 and CEvstates eq 3.0 and CorFeH14 gt -0.2 and CorFeH14 lt 0.0)
    if startype eq 'highalpha14CRGBlow' then good=where(CoralphaFe14 gt 0.15 and CEvstates eq 3.0 and CorFeH14 gt -1.0 and CorFeH14 lt -0.8)
    if startype eq 'highalpha14CRGB' then good=where(CoralphaFe14 gt 0.15 and CEvstates eq 3.0)
    if startype eq 'highalpha14CRGB' then good=where(CoralphaFe14 gt 0.15 and CEvstates eq 3.0)
    if startype eq 'highalpha' then good=where(CoralphaFe13 gt 0.1)
    if startype eq 'lowalpha' then good=where(CoralphaFe13 le 0.1)
    if startype eq 'N6791' then good=where(APGLocID eq 4262 and CorFeH gt 0.15 and CorFeH lt 0.45)
    if startype eq 'real6819' then match2, KICID, real6819, n6819found, good
    if startype eq 'real6791' then match2, KICID, real6791, n6819found, good
    if startype eq 'MVmassivecrossers' then match2, KICID, MVmassivecrossers, MVfound, good

    if startype eq 'Izzardalpharich' then good=where(CoralphaFe16 gt (-0.06*CorFeH16+0.1))
    if startype eq 'highalpha16' then good=where(CoralphaFe16 gt (0.1))
    if startype eq 'highalpha16RGB' then good=where(CoralphaFe16 gt (0.1) and A3CEvstates eq 1)
    if startype eq 'lowalpha16RGB' then good=where(CoralphaFe16 lt (0.13) and A3CEvstates eq 1)
    if startype eq 'youngalpharich16' then good=where(CoralphaFe16 gt (0.1) and AgeGS lt 5.0 and AgeGS gt 0.)
    if startype eq 'fehlt05_16' then good=where(CorFeH16 lt (-0.5))
    if startype eq 'fehlt03_16' then good=where(CorFeH16 lt (-0.3))
    if startype eq 'CRGBm09' then good=where(CEvstates eq 3.0 and EMPSeisMscale gt 0.8 and EMPSeisMscale lt 1.0)
    if startype eq 'CRGBm11' then good=where(CEvstates eq 3.0 and EMPSeisMscale gt 1.0 and EMPSeisMscale lt 1.2)
    if startype eq 'CRGBm13' then good=where(CEvstates eq 3.0 and EMPSeisMscale gt 0.2 and EMPSeisMscale lt 1.4)
    if startype eq 'CRGBm15' then good=where(CEvstates eq 3.0 and EMPSeisMscale gt 0.4 and EMPSeisMscale lt 1.6)

if readoccam eq 'on' then begin
occam=mrdfits('~/Downloads/occam_member-DR16.fits', 1, header)          
oapg=occam.APOGEE_ID
thres=0.0 
ocg=where(occam.CLUSTER eq 'NGC 6819' and occam.RV_PROB gt thres and occam.FEH_PROB gt thres and occam.PM_PROB gt thres)
d6819b=oapg(ocg)                                                                                                              
ocg=where(occam.CLUSTER eq 'NGC 6811' and occam.RV_PROB gt thres and occam.FEH_PROB gt thres and occam.PM_PROB gt thres)
d6811b=oapg(ocg)   
ocg=where(occam.CLUSTER eq 'NGC 6791' and occam.RV_PROB gt thres and occam.FEH_PROB gt thres and occam.PM_PROB gt thres)
d6791b=oapg(ocg)   
endif

    if startype eq 'd6791' then match2, TMASSID, d6791b, n6791found, good
;    if startype eq 'd6791a' then match2, KICID, d6791, n6791found, good
    if startype eq 'd6819' then match2, TMASSID, d6819b, n6819found, good
    if startype eq 'd6811' then match2,TMASSID, d6811b, n6811found, good

    if startype eq 'N6819' then good=where(APGLocID eq 4263 and CorFeH gt -0.05 and CorFeH lt 0.25)
    if startype eq 'logggt25' then good=where(SeisLogg gt 2.5)
    if startype eq 'solar' then good=where(SeisM gt 0.8 and SeisM lt 1.2 and CorFeH gt -0.3 and CorFeH lt 0.3)
    if startype eq 'solarFeH' then good=where(CorFeH13 gt -0.1 and CorFeH13 lt 0.1)
    if startype eq 'solarRGB' then good=where(SeisM gt 0.8 and SeisM lt 1.2 and CorFeH gt -0.3 and $
                              CorFeH lt 0.3 and SeisEvstates gt 2.5 and SeisEvstates lt 3.5)
    if startype eq 'solarFeHRGB' then good=where(CorFeH12 gt -0.3 and CorFeH12 lt 0.3 and SeisEvstates ne 2.0 and ( SeisLogg lt 2.3 or Seislogg gt 3.0 or SeisEvstates eq 3.0))
    if startype eq '13MsunRGB' then good=where(SeisM gt 1.1 and SeisM lt 1.5 and SeisEvstates ne 2.0 and ( SeisLogg lt 2.3 or Seislogg gt 3.0 or SeisEvstates eq 3.0))
    if startype eq 'MLRGB' then good=where(SeisEvstates ne 2.0 and ( SeisLogg lt 2.3 or Seislogg gt 3.0 or SeisEvstates eq 3.0))

    if startype eq 'likelyclump' then good=where( SeisLogg gt 2.3 and SeisLogg lt 3.0 and TeffAPGmMod gt (CorFeH*280.0+120))
    if startype eq 'likelyRGB' then good=where(SeisLogg lt 2.3 or SeisLogg gt 3.0 or (SeisLogg gt 2.3 and SeisLogg lt 3.0 and TeffAPGmMod lt (CorFeH*280.0+120)))
    if startype eq 'vsini13' then good=where(vsinimine13cor gt 0.0)
    if startype eq 'vsini13good' then good=where(vsinimine13cor gt 5.0 and vsinimine13 gt vsinimine13err)
    if startype eq 'turnoff' then good=where(Logg13 ge 3.5 and Logg13 le 4.0)
    if startype eq 'hotstarCNlowfast' then good=where(SeisM14scale gt 1.5 and CN14 lt -0.25 and ((vsinimine13cor gt 7 and vsinimine13cor gt vsinimine13err) or (velocityT14 gt 7 and ProtT12 gt 30)) and (EMPSeisMscale gt 1.5 or EMPSeisMscale lt 0))
    if startype eq 'hotstars' then good=where(SeisM14scale gt 1.5 and ((vsinimine13cor gt 5 and vsinimine13cor gt vsinimine13err) or (ProtT12 gt 30)) and (EMPSeisMscale gt 1.5 or EMPSeisMscale lt 0))
    if startype eq 'hotstarCNlow' then good=where(SeisM14scale gt 1.5 and CN14 lt -0.25)
    if startype eq 'fehm03' then good=where(CorFeH14 gt -0.4 and CorFeH14 lt -0.2)
    if startype eq 'fehm025' then good=where(CorFeH14 gt -0.3 and CorFeH14 lt -0.2)
    if startype eq 'hotstartoofast' then match2, KICint, hottoofast, hottoofastb, good 
    if startype eq 'noseismo' then match2, KICint, noseismo, noseismofound, good   
    if startype eq 'noseismobinary' then begin
	match2, KICint, noseismo, noseismofound, gooda
        goodb=where(noseisreason(gooda) gt 0.5 and noseisreason(gooda) lt 4.5)
        good=gooda(goodb)
    endif
    if startype eq 'noseismorotation' then begin
	match2, KICint, noseismo, noseismofound, gooda
        goodb=where(noseisreason(gooda) gt 4.5 and noseisreason(gooda) lt 7.5)
        good=gooda(goodb)
    endif
    if startype eq 'noseismocontam' then begin
	match2, KICint, noseismo, noseismofound, gooda
        goodb=where(noseisreason(gooda) gt 7.5 and noseisreason(gooda) lt 10.5)
        good=gooda(goodb)
    endif
    if startype eq 'noseismorest' then begin
	match2, KICint, noseismo, noseismofound, gooda
        goodb=where(noseisreason(gooda) gt 10.5 and noseisreason(gooda) lt 14.5)
        good=gooda(goodb)
    endif
    if startype eq 'seismo' then good=where(CEvstates gt 0.0)
;    print, 'match2', good, noseismofound(10023)
;    print, KICID(10023), ' ',noseismo(81)
;    match2, KICID, noseismo, check1, check1g
;    print, 'directmatch', check1g
;print, 'match', where(noseismo eq 8718273),  where(noseismo eq '8718273')
;print, 'wtfmatch', where (noseismo eq KICID(10023)), where(KICID eq noseismo(81))

    if startype eq 'lithiumselect' then good=where(CEvstates gt 0 and CorTeff14 gt 4800 and CorTeff14 lt 4900 and EMPSeisMscale gt 1.0 and EMPSeisMscale lt 2.2)

    if startype eq 'coresurface' then match2, KICID, coresurface, csfound, good
    if startype eq 'halo' then match2, KICID, halo, halofound, good
    if startype eq 'weird' then match2, KICID, weirds, weirdfound, good
    if startype eq 'weirdfast' then match2, KICID, weirds, weirdfound, good
    if startype eq 'GO1' then match2, KICID, GO1, GO1found, good
    if startype eq 'highmal' then match2, KICID, highmal, malfound, good
    if startype eq 'GO1fast' then match2, KICID, GO1, GO1found, good
    if startype eq 'rafa' then match2, KICID, rafaspots, rafafound, good
    if startype eq 'fake5gym6' then match2, KICID, fake5gym6, fakefound, good
    if startype eq 'ericaCN' then match2, KICID, ericaCN, ericafound, good
    if startype eq 'visualoutlierA3' then match2, KICID, visualoutlierA3, voA3found, good
    if startype eq 'ruidehmergercandidates' then match2, KICID, ruidehmergercandidates, mergefound, good

    if startype eq 'scaleMbad' then good=where(abs(SeisMscale-SeisM) gt 0.4)
    if startype eq 'Teffoutlier' then match2, KICID, Teffoutlier, Tefffound, good
    if startype eq 'chi2bad' then good=where(ASPCAPchi2 gt 10)
    if startype eq 'therest' then good=where(ASPCAPchi2 lt 10. and SpecVsini lt 7.0)
    if startype eq 'myvsinilt10' then good=where(vsinimine lt 10.)
    if startype eq 'vsiniboxoutlier' then match2, KICID, vsiniboxoutlier, vsinifound, good
    if startype eq 'vsiniboxoutlier2' then match2, KICID, vsiniboxoutlier2, vsinifound2, good
    if startype eq 'newboxoutlier' then match2, KICID, newboxoutlier, vsinifound3, good
    if startype eq 'vsiniweird' then match2, KICID, allvsiniweird, vsinifound3, good
    if startype eq 'vsiniweirdvetted' then match2, KICID, vsiniweirdvetted, vsinifound4, good
    if startype eq 'vsiniweirdsigma' then match2, KICID, vsiniweirdsigma, vsinifound5, good
    if startype eq 'T12' then good=where(ProtT12 gt 0)
    if startype eq 'T12APG' then good=where(ProtT12 gt 0 and (strcmp(TMASSID, '2M', 2) eq 1))
    if startype eq 'T12multipleAPG' then good=where(ProtT12 gt 0 and Vscatt gt 0)
    if startype eq 'T12possibleRV' then good=where(ProtT12 gt 0 and Vscatt gt 0.5)
    if startype eq 'T12suggestiveRV' then good=where(ProtT12 gt 0 and Vscatt gt 0.25 and Vscatt le 0.5)
    if startype eq 'J12' then good=where(vsinimine12 gt 5.0)
    if startype eq 'vT12lt30' then good=where(velocityT12 lt 30.0) 
    if startype eq 'vsiniminegt5' then good=where(vsinimine gt 5.0 and vsinimineerr lt vsinimine)
    if startype eq 'vsinimine12gt5' then good=where(vsinimine12 gt 5.0 and vsinimine12err lt vsinimine12)
    if startype eq 'J10andT12' then good=where(ProtT12 gt 0 and vsinimine gt 5.0 and vsinimineerr lt vsinimine)
    if startype eq 'J12andT12' then good=where(ProtT12 gt 0 and vsinimine12 gt 5.0 and vsinimine12err lt vsinimine12)
    if startype eq 'J12orT12' then good=where(ProtT12 gt 0 or vsinimine12 gt 5.0)
    if startype eq 'J12notT12' then good=where(ProtT12 lt 0 and vsinimine12 gt 5.0)
    if startype eq 'T12notJ12' then good=where(ProtT12 gt 0 and vsinimine12 lt 5.0)
    if startype eq 'T12gt5notJ12' then good=where(velocityT12 gt 5.0 and vsinimine12 lt 5.0 and ProtT12 gt 0)
    if startype eq 'Mgt3' then good=where(SeisM gt 3.0)
    if startype eq 'vT12gt30' then good=where(velocityT12 gt 30.0 and ProtT12 gt 0)
    if startype eq 'gyrosample' then match2, KICID, gyrosample, gyrofound, good
    if startype eq 'CorsaroRC' then match2, KICID, CorsaroRC, Corfound, good
    if startype eq 'Corsaro6811' then match2, KICID, Corsaro6811, Corfound, good
    if startype eq 'Corsaro6819' then match2, KICID, Corsaro6819, Corfound, good
    if startype eq 'Corsaro6791' then match2, KICID, Corsaro6791, Corfound, good
    if startype eq 'paulstars' then match2, KICID, paulstars, paulfound, good

    if startype eq 'lowZweird' then match2, KICID, lowZweird, zfound, good
    if startype eq 'Enricostars' then match2, KICID, Enricostars, Efound, good
    if startype eq 'Deh15' then match2, KICID, Deh15, dfound, good
    if startype eq 'T12final' then match2, KICID, T12final, tffound, good
    if startype eq 'T12can' then match2, KICID, T12can, tcfound, good
    if startype eq 'comparison1' then match2, KICID, comparison1, cpfound, good
    if startype eq 'comparisonoct' then match2, KICID, comparisonoct, cofound, good
    if startype eq 'coleman' then match2, KICID, coleman, colfound, good
    if startype eq 'tanda' then match2, KICID, tanda, tanfound, good

    if startype eq 'delnureject' then match2, KICID, delnureject, tcfound, good
    if startype eq 'numaxreject' then match2, KICID, numaxreject, tcfound, good

    if startype eq 'Jasonpulse' then match2, KICID, Jasonpulse, Jfound, good   
    if startype eq 'T12P12likelyclump' then good=where(SeisMscaleP12 lt 2.0 and SeisMscaleP12 gt 0.5 and SeisLoggscaleP12 gt 2.3 and SeisLoggscaleP12 lt 2.6 and ProtT12 gt 0)
    if startype eq 'P12likelyclump' then good=where(SeisMscaleP12 lt 2.0 and SeisMscaleP12 gt 0.5 and SeisLoggscaleP12 gt 2.3 and SeisLoggscaleP12 lt 2.6)
    if startype eq 'T12P12likely2RC' then good=where(SeisMscaleP12 gt 2.0 and SeisMscaleP12 lt 3.0 and SeisLoggscaleP12 gt 2.4 and SeisLoggscaleP12 lt 2.8 and ProtT12 gt 0)
    if startype eq 'P12likely2RC' then good=where(SeisMscaleP12 gt 2.0 and SeisMscaleP12 lt 3.0 and SeisLoggscaleP12 gt 2.4 and SeisLoggscaleP12 lt 2.8)
    if startype eq 'supressedl2' then match2, KICID, supressedl2, sl2found, good
    if startype eq 'vT12lt200' then good=where(velocityT12P12 lt 200.0) 
    if startype eq 'vT12lt50' then good=where(velocityT12P12 lt 50.0) 
    if startype eq 'vT12gt10' then good=where(velocityT12P12 gt 10.0) 
    if startype eq 'goodvT12gt10' then good=where(velocityT12P12 gt 10.0 and velocityT12P12 lt 100.0)
    if startype eq 'goodvT14gt10' then good=where(velocityT14 gt 10.0 and velocityT14 lt 100.0)
    if startype eq 'goodvT14gt10m15cn25' then good=where(velocityT14 gt 10.0 and velocityT14 lt 100.0 and SeisM14scale gt 1.5 and CN14 lt -0.25 and ProtT12 gt 30)
    if startype eq 'mfitcheck' then good=where(Evstate eq 'RGB' and CorFeH gt -0.69 and CorFeH lt -0.51 )
    if startype eq 'SeisMscaleP12gt2' then good=where(SeisMscaleP12 gt 2.0) 
    if startype eq 'Kelt' then good=where(CorLogg12 lt 1.5 and KICgMag lt 12) 
    if startype eq 'mas2rc' then good=where(CorLogg12 gt 2.5 and CorLogg12 lt 3.2 and CorTeff12 gt 4600 and CorTeff12 lt 5100 and APGLocID eq 4401)
    if startype eq '4401' then good=where(APGLocID eq 4401)
    if startype eq '4403' then good=where(APGLocID eq 4403)
    if startype eq '4406' then good=where(APGLocID eq 4406)
    if startype eq 'unbiasedbox' then good=where(CorTeff12 gt 4900 and CorTeff12 lt 5200 and CorLogg12 gt 2.8 and CorLogg12 lt 3.2) 
    if startype eq 'vsiniweirdsigmanotRGB' then begin
         match2, KICID, vsiniweirdsigma, vsinifound5, good2b
         good2=where(SeisEvstates lt 2.5 or SeisEvstates gt 3.5)
         match2, good2b, good2, b2found, goodfound
        checkgoodf=where (goodfound ne -1)
        goodfound=goodfound(checkgoodf)
         good=good2b(goodfound)
    endif
    if startype eq 'vsiniweirdsigmaclump' then begin
         match2, KICID, vsiniweirdsigma, vsinifound5, good2b
         good2=where(SeisEvstates gt 1.5 and SeisEvstates lt 2.5)
         match2, good2b, good2, b2found, goodfound
        checkgoodf=where (goodfound ne -1)
        goodfound=goodfound(checkgoodf)
         good=good2b(goodfound)
    endif
    if startype eq 'ylowz' then match2, KICID, yvonnelowz, ylowz, good
    if startype eq 'myvsinilt5' then good=where(vsinimine lt 5.)
    if startype eq 'myvsinigt5' then good=where(vsinimine gt 5.)
    if startype eq 'myvsini12lt5' then good=where(vsinimine12 lt 5.)
    if startype eq 'myvsini12gt5' then good=where(vsinimine12 gt 5.)
    if startype eq 'solzmyvsini12gt5' then good=where(vsinimine12 gt 5. and CorFeH12 gt -0.4 and CorFeH12 lt 0.4)
    if startype eq 'myvsini12gt10' then good=where(vsinimine12 gt 10.)
    if startype eq 'myvsinigt5notRGB' then begin
         good2b=where(vsinimine gt 5.)
         good2=where(SeisEvstates lt 2.5 or SeisEvstates gt 3.5)
         match2, good2b, good2, b2found, goodfound
        checkgoodf=where (goodfound ne -1)
        goodfound=goodfound(checkgoodf)
         good=good2b(goodfound)
    endif
    if startype eq 'myvsinigt5clump' then begin
         good2b=where(vsinimine gt 5.)
         good2=where(SeisEvstates gt 1.5 and SeisEvstates lt 2.5)
         match2, good2b, good2, b2found, goodfound
        checkgoodf=where (goodfound ne -1)
        goodfound=goodfound(checkgoodf)
         good=good2b(goodfound)
    endif
    if startype eq 'vscatter05' then good=where(vscatter gt 0.5)
    if startype eq 'rotation' then good=where(SeisDataflag eq 9)
    if startype eq 'APG' then good=where(strcmp(TMASSID, '2M', 2) eq 1)
    if startype eq 'fast2' then match2, KICID, fast2, fast2found, good
    if startype eq 'mgt15' then good=where(SeisM gt 1.5)
    if startype eq 'mgt20' then good=where(SeisM gt 2.0)
    if startype eq 'mgt23' then good=where(SeisM gt 2.3)  
    if startype eq 'm22_24' then good=where(SeisM gt 2.2 and SeisM lt 2.4)
    if startype eq 'm24_26' then good=where(SeisM gt 2.4 and SeisM lt 2.6)
    if startype eq 'm26_28' then good=where(SeisM gt 2.6 and SeisM lt 2.8)
    if startype eq 'mgt20notRGB' then good=where(SeisM gt 2.0 and SeisEvstates lt 2.6) 
    if startype eq 'mgt23notRGB' then good=where(SeisM gt 2.3 and SeisEvstates lt 2.6)
    if startype eq 'mlt30' then good=where(SeisM lt 3.0)  
    if startype eq 'mlt30gt20' then good=where(SeisM lt 3.0 and SeisM gt 2.0) 
    if startype eq 'mlt30gt23' then good=where(SeisM lt 3.0 and SeisM gt 2.3) 
    if startype eq 'mlt20' then good=where(SeisM lt 2.0) 
    if startype eq 'mgt22' then good=where(SeisM gt 2.2)    
    if startype eq 'mgt22lt30' then good=where(SeisM gt 2.2 and SeisM lt 3.0)
    if startype eq 'mlt15' then good=where(SeisM lt 1.5)
    if startype eq 'mlt13' then good=where(SeisM lt 1.3)
    if startype eq 'alphamweird' then good=where(SeisM gt 1.3 and CoralphaFe gt 0.1)
    if startype eq 'gt13gy' then good=where(Agemod gt 13.0)
    if startype eq 'ZrichOld' then good=where(CorFeH gt 0.3 and Agemod gt 13.0 and CoralphaFe lt 0.1)
    if startype eq 'newalphaM' then match2, KICID, newalphaM, alphamfound, good
    if startype eq 'ProtMc' then good=where(ProtMc gt 0.0)
    if startype eq 'inbox' then good=where(ylist lt m1*(xlist-bx3)+by3 and ylist gt m1*(xlist[*]-bx1)+by1 and ylist lt m2*(xlist-bx3)+by3 and ylist gt m2*(xlist[*]-bx1)+by1)
    if startype eq 'inboxmgt20' then good=where(ylist lt m1*(xlist-bx3)+by3 and ylist gt m1*(xlist[*]-bx1)+by1 and ylist lt m2*(xlist-bx3)+by3 and ylist gt m2*(xlist[*]-bx1)+by1 and SeisM gt 2.0)
bovyint=2.4  ;2.3
bovygcut=2.7  ;2.75
bovyslope=0.6 ;0.527 ;1
    if startype eq '13bovy' then good=where(CorLogg13 lt 0.0018*(CorTeff13-(-382.5*bovyslope*CorFeH13+4607))+bovyint and CorLogg13 gt bovygcut)
    if startype eq '13bovy4401' then good=where(CorLogg13 lt 0.0018*(CorTeff13-(-382.5*bovyslope*CorFeH13+4607))+bovyint and CorLogg13 gt bovygcut and APGLocID eq 4401)
    if startype eq '13bovy4408' then good=where(CorLogg13 lt 0.0018*(CorTeff13-(-382.5*bovyslope*CorFeH13+4607))+bovyint and CorLogg13 gt bovygcut and APGLocID eq 4408)

    if startype eq '13bovymgt20' then good=where(CorLogg13 lt 0.0018*(CorTeff13-(-382.5*bovyslope*CorFeH13+4607))+bovyint and CorLogg13 gt bovygcut and OCTSeisMscale gt 2.0)
    if startype eq '13bovymgt23' then good=where(CorLogg13 lt 0.0018*(CorTeff13-(-382.5*bovyslope*CorFeH13+4607))+bovyint and CorLogg13 gt bovygcut and OCTSeisMscale gt 2.3)
    if startype eq '13bovySmgt20' then good=where(CorLogg13 lt 0.0018*(CorTeff13-(-382.5*bovyslope*CorFeH13+4607))+bovyint and CorLogg13 gt bovygcut and SYDSeisMscale gt 2.0)
    if startype eq '13bovySmgt23' then good=where(CorLogg13 lt 0.0018*(CorTeff13-(-382.5*bovyslope*CorFeH13+4607))+bovyint and CorLogg13 gt bovygcut and SYDSeisMscale gt 2.3)

    if startype eq 'LCflickerbad' then good=where(CEvstates eq 3.0 and F8 lt (-0.2*EMPSeisLoggscale+0.78) and EMPSeisLoggscale gt 2.5)
    if startype eq 'OCTSeisMscalegt20' then good=where(OCTSeisMscale gt 2.0)
    if startype eq 'OCTSeisMscalegt23' then good=where(OCTSeisMscale gt 2.3)

    if startype eq 'bovy' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut)
    if startype eq 'bovymgt20' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and SeisM gt 2.0)
    if startype eq 'bovymSgt20' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and SeisMS gt 2.0)
    if startype eq 'bovymgt23' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and SeisM gt 2.3)

    if startype eq 'bovymSgt23' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and SeisMS gt 2.3)
    if startype eq 'bovy4401' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and APGLocID eq 4401)
    if startype eq '14bovy5003' then good=where(CorLogg14 lt 0.0018*(CorTeff14-(-382.5*bovyslope*CorFeH14+4607))+bovyint and CorLogg14 gt bovygcut and APGLocID eq 5003)
    if startype eq '14bovy5010' then good=where(CorLogg14 lt 0.0018*(CorTeff14-(-382.5*bovyslope*CorFeH14+4607))+bovyint and CorLogg14 gt bovygcut and APGLocID eq 5010)
    if startype eq 'bovy4403' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and APGLocID eq 4403)
    if startype eq 'bovy4406' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and APGLocID eq 4406)
    if startype eq 'bovy4408' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and APGLocID eq 4408)

    if startype eq 'bovy4401dr10' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and APGLocID eq 4401 and CorLogg gt -1)
    if startype eq 'bovy4406dr10' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and APGLocID eq 4406 and CorLogg gt -1)
    if startype eq 'bovy44014406mgt20' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and SeisM gt 2.0 and (APGLocID eq 4401 or APGLocID eq 4406))
    if startype eq 'bovy44014406mgt23' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and SeisM gt 2.3 and (APGLocID eq 4401 or APGLocID eq 4406))
    if startype eq 'bovy44014406' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and (APGLocID eq 4401 or APGLocID eq 4406))
    if startype eq '44014406mgt20' then good=where(SeisM gt 2.0 and (APGLocID eq 4401 or APGLocID eq 4406))
    if startype eq '4401mgt20' then good=where(SeisM gt 2.0 and APGLocID eq 4401 )
    if startype eq '44014406mSgt20' then good=where(SeisMS gt 2.0 and (APGLocID eq 4401 or APGLocID eq 4406))
    if startype eq '4401mSgt20' then good=where(SeisMS gt 2.0 and APGLocID eq 4401 )
    if startype eq 'bovy44014406mSgt20' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and SeisMS gt 2.0 and (APGLocID eq 4401 or APGLocID eq 4406))
    if startype eq 'bovy44014403mSgt20' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and SeisMS gt 2.0 and (APGLocID eq 4401 or APGLocID eq 4403))
    if startype eq 'bovy4401mSgt20' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and SeisMS gt 2.0 and (APGLocID eq 4401 ))
    if startype eq 'bovy4403mSgt20' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and SeisMS gt 2.0 and (APGLocID eq 4403 ))
    if startype eq 'bovy4406mSgt20' then good=where(CorLogg12 lt 0.0018*(CorTeff12-(-382.5*bovyslope*CorFeH12+4607))+bovyint and CorLogg12 gt bovygcut and SeisMS gt 2.0 and (APGLocID eq 4406 ))
    if startype eq 'F8good' then good=where(F8seislogg lt 4.6 and F8seislogg gt 2.5)
    if startype eq 'APGvsini14lt10' then good=where(APGvsini14 lt 10)

    if startype eq '44014406' then good=where( (APGLocID eq 4401 or APGLocID eq 4406))

    if startype eq 'A3Pgood' then good=where(abs(A3PSeisLoggscale-CorLogg16) lt 0.5)
    if startype eq 'A3Pclumplogg' then good=where(A3PSeisLoggscale lt 2.6 and A3PSeisLoggscale gt 2.3)
    if startype eq 'A3Pnotclumplogg' then good=where(A3PSeisLoggscale gt 2.6 or A3PSeisLoggscale lt 2.3)
    if startype eq 'A3PGiraffe' then good=where(A3PSeisMscale gt 1.8 and A3PSeisMscale lt 2.0 and CorFeH16 lt -0.2)
    if startype eq 'speclogglt35' then good=where(CorLogg16 lt 3.5)
;oplot, CorTeff12, 0.0018*(CorTeff12-(-382.5*0.5*(-0.4)+4607))+bovyint, psym=3

    checkgood=where (good ne -1)
    good=good(checkgood)
    if startype eq 'GO1fast' then begin

          good2=where(SpecVsini(good) gt 5.0 and SpecVsini(good) ne 9999)
          good=good(good2)
    endif
    if startype eq 'weirdfast' then begin
          good2=where(SpecVsini(good) gt 10.0 and SpecVsini(good) ne 9999)
          good=good(good2)
    endif
    if startype eq 'CRGBlogggt3' then good=where(CEvstates eq 3.0 and EMPSeisLoggscale gt 3.0)
    if startype eq 'CRGBlogggt28' then good=where(CEvstates eq 3.0 and EMPSeisLoggscale gt 2.8)
    if startype eq 'CRGBlogglt23gt1' then good=where(CEvstates eq 3.0 and EMPSeisLoggscale gt 1.0 and EMPSeisLoggscale lt 2.3)
    if startype eq 'EMPlogggt3' then good=where( EMPSeisLoggscale gt 3.0)
    if startype eq 'EMPlogglt2gt1' then good=where(EMPSeisLoggscale gt 1.0 and EMPSeisLoggscale lt 2.0)

    if startype eq 'A3CRGBlowalpha' then good=where(A3CEvstates eq 1 and CoralphaFe16 lt 0.07 )
    if startype eq 'A3CRGBlogggt3' then good=where(A3CEvstates eq 1 and A3PSeisLoggscale gt 3.0 )
    if startype eq 'A3CRGBlogglt2gt1' then good=where(A3CEvstates eq 1 and A3PSeisLoggscale gt 1.0 and A3PSeisLoggscale lt 2.0 )
    if startype eq 'A3Cnotclumplogggt3' then good=where(A3CEvstates ne 2 and A3PSeisLoggscale gt 3.0 )
    if startype eq 'A3Cnotclumplogglt2gt1' then good=where(A3CEvstates ne 2 and A3PSeisLoggscale gt 1.0 and A3PSeisLoggscale lt 2.0 )
    if startype eq 'solar16' then good=where(CorFeH16 gt -0.1 and CorFeH16 lt 0.1)
    if startype eq 'A3CnotclumpSolar' then good=where(A3CEvstates lt 1.5 and CorFeH16 gt -0.1 and CorFeH16 lt 0.1)
    if startype eq 'CRGBSolar' then good=where(CEvstates eq 3.0 and CorFeH16 gt -0.1 and CorFeH16 lt 0.1)
    if startype eq 'A3CRGBSolar' then good=where(A3CEvstates eq 1 and CorFeH16 gt -0.1 and CorFeH16 lt 0.1)
    if startype eq 'A3CRGBSolarm15' then good=where(A3CEvstates eq 1 and CorFeH16 gt -0.1 and CorFeH16 lt 0.1 and A3PSeisMscale gt 1.4 and A3PSeisMscale lt 1.6)
    if startype eq 'A3CRGBSolarm10' then good=where(A3CEvstates eq 1 and CorFeH16 gt -0.1 and CorFeH16 lt 0.1 and A3PSeisMscale gt 0.9 and A3PSeisMscale lt 1.1)
    if startype eq 'A3CRGBm10m030' then good=where(A3CEvstates eq 1 and CorFeH16 gt -0.4 and CorFeH16 lt -0.2 and A3PSeisMscale gt 0.9 and A3PSeisMscale lt 1.1)
    if startype eq 'A3CRGBm10p030' then good=where(A3CEvstates eq 1 and CorFeH16 gt 0.2 and CorFeH16 lt 0.4 and A3PSeisMscale gt 0.9 and A3PSeisMscale lt 1.1)
    if startype eq 'MHbinary' then good=where(CorFeH16 gt 0.05 and CorFeH16 lt 0.15 and A3PSeisMscale gt 1.1 and A3PSeisMscale lt 1.3)
    if startype eq 'A3Cnotclumpm04' then good=where(A3CEvstates lt 1.5 and CorFeH16 gt -0.5 and CorFeH16 lt -0.3)
    if startype eq 'A3Callclump' then good=where(A3CEvstates eq 2)
    if startype eq 'A3Cnotclump' then good=where(A3CEvstates ne 2)
    if startype eq 'A3CRGB' then good=where(A3CEvstates eq 1)
    if startype eq 'A3CU' then good=where(A3CEvstates eq -1)
    if startype eq 'PWbinary' then good=where(PWbinaryP gt 0)
    if startype eq 'A3highmRGB' then good=where(A3CEvstates eq 1 and abs(A3PSeisLoggscale-CorLogg16) lt 0.5 and A3PSeisMscale gt 2.0 and A3PSeisMscale lt 5.0)
    if startype eq 'A3Callclumpfehs' then good=where(A3CEvstates eq 2 and (CorFeH16 lt -0.3 or CorFeH16 gt 0.3))

    if startype eq 'pgbestrotgt0' then good=where(pgbestrot gt 0)
   if startype eq 'A3Pnumaxgt10' then good=where(A3Pnumax gt 10)
   if startype eq 'A3Pagelt1' then good=where(AgeGS lt 1. and AgeGS gt 0.)
   if startype eq 'A3Pagegt1' then good=where(AgeGS gt 1.)
   if startype eq 'A3Pclumpagegt15' then good=where(AgeGS gt 15. and A3CEvstates eq 2)
   if startype eq 'A3Phighmnotclump' then good=where(A3CEvstates ne 2 and A3PSeisMscale gt 2.0 and A3PSeisMscale lt 5.0)
   if startype eq 'corerotcg' then good=where(corerotcg gt 0)
 ;if startype eq '4401' then

;writecol, 'pgapokasc.txt', KICID(good), masspg(good), Rpg(good), Protpg(good), vexpectedpg(good), pgbestev(good), A3PSeisMscale(good)
;writecol, 'subaru.txt', KICID(good), RA(good), DEC(good), RA(good)*0.0+2000, KepMag(good)


;print, "KICID, TMASSID, SeisM, ApgLocID, SeisLogg, SpecLogg, CorTeff, CorFeH, chi2, vsini, periodmax"
;print, KICID(good), TMASSID(good), APGLocID(good), SeisM(good), SeisLogg(good), CorLogg(good), $
;        CorTeff(good), CorFeH(good), ASPCAPchi2(good), SpecVsini(good), periodmax(good)

;writecol,'textout.txt', "KICID, TMASSID, SeisM, ApgLocID, SeisLogg, SpecLogg, CorTeff" , $
;          "CorFeH, chi2, vsini, periodmax, periodmax3d, Evstate"
;writecol,'lowzout.txt',fmt='(A20,1x, A30,1x, F11.6,1x,F11.6, 1x, F10.6, 1x, F10.6, 1x, F10.6, 1x, F10.6, 1x, F10.6, 1x, F10.6, 1x, F10.6)'  , KICID(good), TMASSID(good), CorTeff(good), CorTefferr(good),CorLogg(good), CorLoggerr(good), CorFeH(good), CorFeHerr(good), CoralphaFe(good), CoralphaFeerr(good)
;writecol,'textout.txt', KICID(good), TMASSID(good), APGLocID(good), SeisM(good), SeisR(good), $
;        SeisLogg(good), CorLogg(good), CorTeff(good), CorFeH(good), ASPCAPchi2(good), $
;        vsinimine(good), periodmax(good), periodmax3d(good), SeisEvstates(good)
;writecol,'textout.txt', KICID(good), TMASSID(good), APGLocID(good), SeisM(good), SeisR(good), $
;        SeisLogg(good), CorLogg(good), CorTeff(good), CorFeH(good), SeisEvstates(good),$
;        numaxRG(good),delnuRG(good)
;        vsinimine(good), vsinimineerr(good), periodmax(good), periodmaxerr(good), periodmax3d(good), ProtC(good), ProtCerr(good), inclination(good), inclinationerr(good)
;writecol,'textout.txt', KICID(good), TMASSID(good), APGLocID(good), SeisM(good), SeisR(good), $
;        SeisLogg(good), Logg12(good), Teff12(good), FeH12(good), $
;        vsinimine(good),  FeH2(good), FeH2err(good), CH(good), CHerr(good), NH(good), NHerr(good), OH(good), OHerr(good)
;writecol,'textout.txt', KICID(good), TMASSID(good), APGLocID(good), SeisM(good), SeisMerr(good), SeisR(good), SeisRerr(good), $
;        SeisLogg(good), SeisLoggerr(good), CorLogg12(good), CorTeff12(good),CorTeff12err(good), CorFeH12(good),CorFeH12err(good),$
;        vsinimine12(good), vsinimine12err(good), periodmax(good), periodmaxerr(good), SeisEvstates(good)
;writecol,'textoutprelim.txt', KICID(good), TMASSID(good), SeisM(good), SeisR(good), $
;        SeisLogg(good), CorLogg13(good),  CorTeff13(good), CorFeH13(good),$
;        CoralphaFe13(good), vsinimine12(good), vsinimine12err(good), $
;        ProtT12(good), ProtT12err(good), SeisEvstates(good), CN(good);, CH(good);, vsiniN(good)
;writecol,'textoutprelim.txt', KICID(good), TMASSID(good), SeisM(good), SeisR(good), $
;        SeisLogg(good), CorLogg13(good),  CorTeff13(good), CorFeH13(good),$
;        CoralphaFe13(good), vsinimine12(good), vsinimine12err(good), $
;        ProtT12(good), ProtT12err(good), SeisEvstates(good), CN(good);, CH(good);, vsiniN(good)
;writecol, 'Dantextoutvscatter.txt', KICID(good), TMASSID(good), N_KEP_QUART(good), CorFeH13(good), CorTeff13(good), CorLogg13(good), F8Logg(good), CN13(good), nvisits(good), vscatter(good), APGvsini14(good), vsinimine(good), vsinimine12(good),ProtT12(good), TARGFLAGS(good), ASPCAPFLAGS(good)
;writecol,'textoutprelim14.txt', KICID(good), TMASSID(good), CorLogg14(good),  CorTeff14(good), CorFeH14(good),$
;        CoralphaFe14(good)
;writecol,'textoutprelim13.txt', KICID(good), TMASSID(good), CorLogg13(good),  CorTeff13(good), CorFeH13(good),$
;        CoralphaFe13(good)
;writecol, 'basucor.txt', KICID(good), OCTSeisMscale(good),CorFeH13(good),SYDSeisMscale(good),SYDSeisMscalecor(good)

;writecol, 'kinout.txt', KICID(good), TMASSID(good), RA12(good), DEC12(good), PropRAMot12(good),PropDecMot12(good), KICEBmV12(good)
;writecol, 'textoutrot.txt', KICID(good), TMASSID(good), vsiniN(good), vsinimine12(good), vsinimine12err(good), $
;        periodmaxP12(good), periodmaxP12err(good), ProtT12(good), ProtT12err(good), inclinationP12(good), inclinationP12err(good)
;writecol, 'textoutrot.txt', KICID(good), TMASSID(good), vsiniN(good), vsinimine(good), vsinimineerr(good), $
;        periodmaxP12v10(good), periodmaxP12v10err(good), ProtT12(good), ProtT12err(good), inclinationP12v10(good), inclinationP12v10err(good)
;writecol, 'textoutrotH14v13.txt', KICID(good), TMASSID(good), APGvsini16(good), vsinimine13cor(good), vsinimine13err(good), $
;        periodmaxH14v13(good), periodmaxH14v13err(good), ProtMc(good), ProtMcerr(good), inclinationH14v13(good), inclinationH14v13err(good)
;writecol, 'diffrotsample.txt', KICID(good), SeisM10(good),vsinimine(good)
;writecol, 'svmEvstate.txt', KICID(good), svmevstate(good),ConEvstates(good)
;writecol, 'T12out.txt', KICID(good), TMASSID(good), P12Teff(good), P12Tefferr(good), CorFeH12(good), CorFeH12err(good), SeisMscaleP12(good), SeisLoggscaleP12(good), ProtT12(good), ProtT12err(good)
;writecol, 'T12dr12.txt', KICID(good), TMASSID(good), delnuRG(good), numaxRG(good), ProtT12(good), ProtT12err(good), SeisM(good), SeisR(good), SeisLogg(good), CorLogg12(good), CorTeff12(good), CorFeH12(good), vsinimine12(good), vsinimine12err(good)
;writecol, 'KELTout.txt', KICID(good), TMASSID(good), RA12(good), DEC12(good), KICgMag(good), CorLogg13(good), numaxT12(good), numaxT12err(good), delnuT12(good), delnuT12err(good), ProtT12(good), ProtT12err(good)
;writecol, 'unbiasedDR12.txt', KICID(good), TMASSID(good), CorTeff12(good),CorFeH12(good),CorLogg12(good),numaxT12(good),delnuT12(good), vsinimine(good), vsinimineerr(good), vsinimine12(good), vsinimine12err(good), vsiniN(good)
;writecol, 'comparisonoct.txt', KICID(good), TMASSID(good), numaxT12(good),delnuT12(good), SeisMS(good), SeisM(good), SeisLogg(good),CorLogg12(good),CorTeff12(good),P12Teff(good)
;writecol, 'comparisonoct2.txt',KICID(good),CorFeH12(good), CoralphaFe12(good), CN(good), ProtT12(good), ProtT12err(good),vsinimine(good), vsinimineerr(good), vsinimine12(good), vsinimine12err(good), vsiniN(good)

;writecol, 'MLtesterr3bdwarf.txt', KICID(good), SeisMdw(good), Logg13(good), CorFeH13(good), CoralphaFe13(good), CorTeff13(good), vsinimine13cor(good), SeisMdwerr(good), CorLogg13err(good), CorFeH13err(good), CoralphaFe13err(good), CorTeff13err(good), TMASSID(good)

;writecol, 'vsini13toprot.txt', KICID(good),TMASSID(good), vsinimine13cor(good),vsinimine13err(good),SeisMdw(good), SeisRdw(good),periodmaxSDSSv13(good), H14M(good), H14R(good), periodmaxH14v13(good)
;writecol, 'period14.txt', KICID(good),TMASSID(good), vsinimine13cor(good),vsinimine13err(good),EMPSeisMscale(good), EMPSeisRscale(good),velocityT14(good), SeisM14scale(good), SeisR14scale(good), ProtT12(good)
;if startype eq 'paulstars' then writecol, 'paulstars.txt', KICID(good), OCTSeisMscale(good), OCTSeisLoggscale(good), CorFeH13(good), CoralphaFe13(good), CorTeff13(good), CN13(good), OSSeisMscaleerr(good), OSSeisLoggscaleerr(good), CorFeH13err(good), CoralphaFe13err(good), CorTeff13err(good), TMASSID(good),vsinimine13cor(good), ProtT12(good),ProtT12err(good),ProtT12allbad(good), OCTSeisRscale(good), CorLogg13(good)

;writecol, 'MLtesterr3bdennis.txt', KICID(good), OCTSeisMscale(good), OCTSeisLoggscale(good), CorFeH13(good), CoralphaFe13(good), CorTeff13(good), CN13(good), OSSeisMscaleerr(good), OSSeisLoggscaleerr(good), CorFeH13err(good), CoralphaFe13err(good), CorTeff13err(good), TMASSID(good),OCTdelnu(good),OCTnumax(good)
;writecol, 'MLtesterr3bSYD.txt', KICID(good), SYDSeisMscale(good), SYDSeisLoggscale(good), CorFeH13(good), CoralphaFe13(good), CorTeff13(good), CN13(good), OSSeisMscaleerr(good), OSSeisLoggscaleerr(good), CorFeH13err(good), CoralphaFe13err(good), CorTeff13err(good), TMASSID(good)
;writecol, 'MLtesterr3bspec.txt', KICID(good), OCTSeisMscale(good), CorLogg13(good), CorFeH13(good), CoralphaFe13(good), CorTeff13(good), CN13(good), OSSeisMscaleerr(good), OSSeisLoggscaleerr(good), CorFeH13err(good), CoralphaFe13err(good), CorTeff13err(good), TMASSID(good)
;writecol, 'MLtesterr3bDR14.txt', KICID(good), EMPSeisMscale(good), EMPSeisLoggscale(good), CorFeH14(good), CoralphaFe14(good), CorTeff14(good), CN14(good), EMPSeisMscaleerr(good), EMPSeisMscaleerr(good)*0.0-9999, CorFeH14err(good), CoralphaFe14err(good), CorTeff14err(good), TMASSID(good)

;writecol, 'MLtesterr3bGHB.txt', KICID(good), OCTSeisMscaleGHB(good), OCTSeisLoggscaleGHB(good), CorFeH13(good), CoralphaFe13(good), GHBteff(good), CN13(good), OSSeisMscaleerr(good), OSSeisLoggscaleerr(good), CorFeH13err(good), CoralphaFe13err(good), CorTeff13err(good), TMASSID(good)
;writecol, 'MLtesterr2Aldo.txt', KICID(good), aldoOCTSeisM(good), OCTSeisLoggscale(good), CorFeH13(good), CoralphaFe13(good), CorTeff13(good), OFe13(good), aldoOCTSeisMerr(good), OSSeisLoggscaleerr(good), CorFeH13err(good), CoralphaFe13err(good), CorTeff13err(good), TMASSID(good)
;writecol, 'MLtesterr3Sarbanimcor.txt', KICID(good), OCTSeisMscalecor(good), OCTSeisLoggscale(good), CorFeH13(good), CoralphaFe13(good), CorTeff13(good), CN13(good), OSSeisMscaleerr(good), OSSeisLoggscaleerr(good), CorFeH13err(good), CoralphaFe13err(good), CorTeff13err(good), TMASSID(good)
;writecol, 'MLtesterr3Sarbanimloggcor.txt', KICID(good), OCTSeisMscalecor(good), OCTSeisLoggscalecor(good), CorFeH13(good), CoralphaFe13(good), CorTeff13(good), CN13(good), OSSeisMscaleerr(good), OSSeisLoggscaleerr(good), CorFeH13err(good), CoralphaFe13err(good), CorTeff13err(good), TMASSID(good)
;writecol, 'MLtesterr2SYD.txt', KICID(good), SYDSeisMscale(good), OCTSeisLoggscale(good), CorFeH13(good), CoralphaFe13(good), CorTeff13(good), OFe13(good), OSSeisMscaleerr(good), OSSeisLoggscaleerr(good), CorFeH13err(good), CoralphaFe13err(good), CorTeff13err(good), TMASSID(good)
;writecol, 'MLtesterr2seis.txt', KICID(good), OCTdelnu(good),OCTnumax(good), OCTSeisLoggscale(good), CorFeH13(good), CoralphaFe13(good), CorTeff13(good), OFe13(good), OSSeisMscaleerr(good), OSSeisLoggscaleerr(good), CorFeH13err(good), CoralphaFe13err(good), CorTeff13err(good), TMASSID(good)
;writecol, 'GHBout.txt', KICID(good), CorFeH13(good), CoralphaFe13(good), CorTeff13(good),  GHBteff(good), CorLogg13(good), Logg13(good), CorFeH13err(good), CoralphaFe13err(good), CorTeff13err(good),   TMASSID(good)
;writecol, 'ForGridModelingCor.txt', KICID(good),TMASSID(good), Logg13(good), CorLogg13(good), CorLogg13err(good), Teff13(good), CorTeff13(good), CorTeff13err(good), FeH13(good), CorFeH13(good), CorFeH13err(good), alphaFe13(good), CoralphaFe13(good), CoralphaFe13err(good) 
;writecol, 'rotvalues.txt', KICID(good), OCTSeisMscale(good), OCTSeisLoggscale(good), vsinimine13cor(good) ;
;writecol, 'noseismo.txt', KICID(good), TMASSID(good),CorTeff13(good), CorLogg13(good), noseisreason(good)

;writecol, 'loggcomparison_flicker.txt', KICID(good), TMASSID(good), Logg13(good), CorLogg13(good), Logg14(good), CorLogg14(good), CannonLogg(good), F8Logg(good), Teff14(good), FeH14(good)
;writecol, 'textoutrot14.txt', KICID(good), TMASSID(good),APGLocID(good), ProtT12allbad(good),ProtT12(good), ProtT12err(good), nvisits(good),vscatter(good), vsinimine(good), vsinimineerr(good), vsinimine13cor(good), vsinimine13err(good),  APGvsini14(good)

;writecol, 'Flicker_table.txt', KICID(good), TMASSID(good), F8numax(good), F8numaxerr(good), F8seismass(good), F8seismasserr(good), F8teff(good), F8tefferr(good), F8feh(good), F8feherr(good), F8(good), F8err(good), F8Evstates(good), CoralphaFe13(good), CoralphaFe13err(good), CN13(good),CN13err(good)

;writecol2, 'Flicker_table14.txt', KICID(good), TMASSID(good), F8numax(good), F8numaxerr(good), F8seismass(good), F8seismasserr(good), F8teff14(good), F8teff14err(good), F8feh14(good), F8feh14err(good), F8(good), F8err(good), F8Evstates(good), CoralphaFe14(good), CoralphaFe14err(good), CN14(good),CN14err(good), F8logg14(good), F8logg14err(good), vscatter(good)

;writecol2, 'textout16.txt', KICID(good), TMASSID(good), CEvstates(good), ProtT12(good), CorTeff16(good), CorLogg16(good), ylist(good), xlist(good), colorlist(good), vsinimine13(good), APGvsini16(good)
;writecol, 'MLtesterr3bDR16A3PalphaM_july.txt', KICID(good), A3PSeisMscale(good), A3PSeisLoggscale(good), CorFeH16(good), CoralphaFe16(good), CorTeff16(good), CN16(good), A3PSeisMscaleerr(good), A3PSeisLoggscaleerr(good), CorFeH16err(good), CoralphaFe16err(good), CorTeff16err(good), TMASSID(good)

if startype eq 'highalpha16' then writecol, 'youngalpharich_dr16.txt', KICID(good), A3PSeisMscale(good), A3PSeisLoggscale(good), CorFeH16(good), CoralphaFe16(good), CorTeff16(good), CN16(good),TMASSID(good), AgeGS(good), pgbestrot(good), PWbinaryP(good), vscatter(good)

;writecol, 'A3PMLtesterr3bDR16alphaM.txt', KICID(good), A3PSeisMscale(good), A3PSeisLoggscale(good), CorFeH16(good), CoralphaFe16(good), CorTeff16(good), CN16(good), A3PSeisMscaleerr(good), A3PSeisLoggscaleerr(good), CorFeH16err(good), CoralphaFe16err(good), CorTeff16err(good), TMASSID(good)

;writecol, 'A2MLtesterr3bDR16alphaM.txt', KICID(good), EMPSeisMscale(good), EMPSeisLoggscale(good), CorFeH16(good), CoralphaFe16(good), CorTeff16(good), CN16(good), EMPSeisMscaleerr(good), EMPSeisLoggscaleerr(good), CorFeH16err(good), CoralphaFe16err(good), CorTeff16err(good), TMASSID(good)

;writecol, 'A2MLtesterr3bDR14alphaM.txt', KICID(good), EMPSeisMscale(good), EMPSeisLoggscale(good), CorFeH14(good), CoralphaFe14(good), CorTeff14(good), CN14(good), EMPSeisMscaleerr(good), EMPSeisLoggscaleerr(good), CorFeH14err(good), CoralphaFe14err(good), CorTeff14err(good), TMASSID(good)
;writecol, 'MLtesterr3bDR16MgFe.txt', KICID(good), EMPSeisMscale(good), EMPSeisLoggscale(good), FeH16(good), MgFe16(good), CorTeff16(good), CN16(good), EMPSeisMscaleerr(good), EMPSeisLoggscaleerr(good), CorFeH16err(good), CoralphaFe16err(good), CorTeff16err(good), TMASSID(good)
spacearr=make_array(n_elements(KICID), /string, value= ' ')

writecol, 'CNmixing_A2DR14'+startype+'.txt',EMPSeisLoggscale(good), CorFeH14(good), CN14(good), EMPSeisMscale(good),spacearr(good), KICID(good)

if startype eq 'coresurface' then writecol, 'coresurface.txt', KICID(good), TMASSID(good),CorTeff13(good), CorLogg13(good), CorFeH13(good), CoralphaFe13(good), OCTSeisMscale(good),  SeisMscaleP12(good),SeisRscaleP12(good),ProtT12(good), vsinimine13cor(good), periodmaxH14v13(good)
if startype eq 'turnoff' then writecol, 'turnoff.txt', KICID(good), TMASSID(good),Teff13(good), Logg13(good), FeH13(good), vsinimine13cor(good),vsinimine13err(good),SeisMdw(good), SeisRdw(good),periodmaxSDSSv13(good), H14M(good), H14R(good), periodmaxH14v13(good)
if startype eq 'coleman' then writecol, 'colemanout.txt', KICID(good), TMASSID(good),CorTeff13(good), CorLogg13(good), CorFeH13(good), nvisits(good), vscatter(good)
    fraclist=intarr(n_elements(xlist))
    fraclist(good)=1
    if bkg eq 'frac' then fractioncolor, xmin, xmax, ymin, ymax, 30.0, 30.0, xlist, ylist, colorlist,fraclist, KICID=KICID

loadct, 0
; print, [bx1, bx2, bx3, bx4, bx1],[by1, by2,by3, by4, by1], m1,m2
; oplot, [bx1, bx2, bx3, bx4, bx1],[by1, by2,by3, by4, by1], thick=2
; oplot, [5200, 4900, 4900, 5200, 5200], [2.8, 2.8, 3.2, 3.2, 2.8], thick=2
if colortab lt 50 then loadct, colortab
if colortab eq 100 then cgloadct, /brewer,31, /silent

print, 'goods', n_elements(good)
print, 'mean metallicity14(good)', avg(CorFeH14(good))
print, 'mean Teff14(good)', avg(CorTeff14(good))

if label ne 'off' then xyouts,  x2-(.2*xrang), y2-(.2*yrang), n_elements(good), charsize=1, $
         charthick=2
;   print, KICID(good[0:2]), CorFeH12(good[0:2])

 loadct, 0 
 if startype eq 'vsiniweirdsigma' then begin
      match2, KICID, fast2, fast2found, goodrapid
      checkgoodr=where (goodrapid ne -1)
       goodrapid=goodrapid(checkgoodr)
      oplot, xlist(goodrapid), ylist(goodrapid), psym=6

 endif

 if density eq 'off' or density eq 'on' then oplot, xlist(good), ylist(good), psym=1 
 if writename eq 'on' then  xyouts, xlist(good), ylist(good), KICID(good)
 if writeevstate eq 'on' then  begin
	good2=where(Sevstate(good) ne '-9999')
        if good2[0] ne -1 then xyouts, xlist(good(good2)), ylist(good(good2)), Sevstate(good(good2))
	good2=where(Bevstate(good) ne '-9999')
        if good2[0] ne -1 then xyouts, xlist(good(good2)), ylist(good(good2)), Bevstate(good(good2))
	good2=where(Evstates(good) gt -9998)
        if good2[0] ne -1 then xyouts, xlist(good(good2)), ylist(good(good2)), Evstates(good(good2))
 endif
if colortab lt 50 then loadct, colortab
if colortab eq 100 then cgloadct, /brewer,31, /silent
    colorz=(colorlist(good)-min(colorlist))*255/(max(colorlist)-min(colorlist)+0.00001)
    if goodcolors eq 'on' then     colorz=(colorlist(good)-min(colorlist(good)))* $
                255./(max(colorlist(good))-min(colorlist(good)))
    if reversecolors eq 'on' then colorz=255-colorz

    if ploterr eq 'on' then begin
	 loadct, 0
	oploterror, xlist(good),  ylist(good),xlisterr(good), ylisterr(good), psym=1, color=150

	    FITEXY, xlist(good), ylist(good), intercept, slope, X_SIG=xlisterr(good), Y_SIG=ylisterr(good), $
	             sigmaab, chisq, q
;	    oplot, [xmin, xmax], [intercept+xmin*slope, intercept+xmax*slope], linestyle=2, color=0
;	     xyouts,  x1+(.1*xrang), y2-(.1*yrang), "slope="+string(slope)+'+/-'+string(sigmaab[1]), charsize=1, $
;	         charthick=2   
;	     xyouts,  x1+(.1*xrang), y2-(.2*yrang), "intercept="+string(intercept)+'+/-'+string(sigmaab[0]), charsize=1, $
;	         charthick=2   
if colortab lt 50 then loadct, colortab
if colortab eq 100 then cgloadct, /brewer,31, /silent
    endif
    usersym, [-.8, -.8, .8, .8, -.8],[-.8,.8,.8,-.8,-.8], /fill

    if density eq 'off' then begin 
	    plots, xlist(good), ylist(good),psym=7, symsize=1, color= [colorz,colorz,colorz]; was sym(5), symsize2

	    if (goodcolors eq 'off' and label ne 'off' and reversecolors eq 'off') then al_legend, [colorstring+string(max(colorlist)), $
                          colorstring+string(.5*(max(colorlist)+min(colorlist))), $
			colorstring+string(min(colorlist))],$
		        psym=[7,7,7], $
		        colors=[255,127,0],/right_legend, /bottom_legend
	    if (goodcolors eq 'on' and label ne 'off') then  al_legend, [colorstring+string(max(colorlist(good))), $
                          colorstring+string(.5*(max(colorlist(good))+min(colorlist(good)))), $
			colorstring+string(min(colorlist(good)))],$
		        psym=[7,7,7], $
		        colors=[255,127,0],/right_legend, /bottom_legend

	    if (goodcolors eq 'off' and label ne 'off' and reversecolors eq 'on') then al_legend, [colorstring+string(min(colorlist)), $
                          colorstring+string(.5*(max(colorlist)+min(colorlist))), $
			colorstring+string(max(colorlist))],$
		        psym=[7,7,7], $
		        colors=[255,127,0],/right_legend, /bottom_legend

    endif 
    if density eq 'on' then begin
	   colorz=colorz*0.0+255
	   plots, xlist(good), ylist(good),psym=7, symsize=2, color= [colorz,colorz,colorz] 

    endif
   print, 'label', label
   if density eq 'good' then densitycolor, xmin, xmax, ymin, ymax, 20.0, 20.0, xlist(good), ylist(good), colorlist(good), KICID=KICID(good), label
   if density eq 'goodr' then densitycolorrev, xmin, xmax, ymin, ymax, 20.0, 20.0, xlist(good), ylist(good), colorlist(good), KICID=KICID(good)
   if density eq 'goodrstretch' then densitycolorrevstretch, xmin, xmax, ymin, ymax, 20.0, 20.0, xlist(good), ylist(good), colorlist(good), KICID=KICID(good)
    if bkg eq 'sample' then bkgcolor2, xmin, xmax, ymin, ymax,  10.0, 10.0, xlist(good), ylist(good), colorlist(good), colorlisterr(good), KICID=KICID
   if fitline eq 'on' then begin
        linfit, xlist(good), ylist(good), 1.0+0.0*xlist(good), aout, bout
 	oplot, [xmin, xmax], [aout+xmin*bout, aout+xmax*bout]
 	xyouts, xmin+0.1*(xmax-xmin), ymax-0.1*(ymax-ymin), 'y='+strcompress(string(aout), /remove_all)+'+x'+strcompress(string(bout), /remove_all)
   endif
   if runmedian eq 'on' then begin
        loadct, 0
        xmed=xlist(good)
        ymed=ylist(good)
        aorder=sort(xmed)
        xorder=xmed(aorder)
        yorder=ymed(aorder)
        nmpoints=20.
        perpoint=fix(n_elements(xorder)/nmpoints)-1
        for i=0,nmpoints-1  do begin
           xmedplot=median(xorder[i*perpoint:i*perpoint+perpoint-1])
           ymedplot=median(yorder[i*perpoint:i*perpoint+perpoint-1])
           oplot, [xmedplot], [ymedplot], psym=5, symsize=2
        end
;	oplot, [0, xmax], [140, 140], thick=5
;	oplot, [xmin, -0.8], [0, 0], thick=5
;	oplot, [-0.8, 0], [0, 140], thick=5
if colortab lt 50 then loadct, colortab
if colortab eq 100 then cgloadct, /brewer,31, /silent
   endif
   if pltextradata eq 'on' then begin
        loadct, 34
        xextra=[.7, .8, 1.1, 1.6, 2.1, 2.5, 3.5, 8.0] ;cluster ages myr
        xextraerr=xextra*0.0
        yextra=[-.73, -.55, -0.6, -.49, -.36, -.51, -0.37, -.25] ; cluster cns 
        yextraerr=[.3/sqrt(2), .2/sqrt(2),  .2/sqrt(12), .17/sqrt(5), .2/sqrt(5), .25/sqrt(10), .4/sqrt(11), .2/sqrt(6)]; cluster cn err
	oploterror, xextra(good),  yextra(good),xextraerr(good), yextraerr(good), psym=1, color=255

if colortab lt 50 then loadct, colortab
if colortab eq 100 then cgloadct, /brewer,31, /silent
; fit Salaris 2015 relation
        cn=[-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0]
        met=0
        ageSalaris=(34.99+156.3*cn+20.52*met+298.6*(cn^2.)+78.16*cn*met+7.82*(met^2.)+ $
                 305.9*(cn^3.)+85.02*(cn^2.)*met+ $
                  22.93*cn*(met^2.)+0.8987*(met^3.)+$
                 141.1*(cn^4.)+21.96*(cn^3.)*met+$
                 16.14*(cn^2.)*(met^2.)+1.447*cn*(met^3.))
        oplot, ageSalaris, cn

    endif
    if pltextradata eq 'MH' then begin 
        readcol, '~/Documents/MixingLength/Garrett/new_outfiles/mhages_dr17al_specloggteff5.out', kic, mass, logg, feh, alfe, teff, cn, merr, loggerr, feherr, alfeerr, tefferr, tmass, intteff, intml, intage, intc12, intc13, intn14, intxsurf, format='A,F,F,F,F,F,F,F,F,F,F,F,A,F,F,F,F,F,F,F'
       loadct, 34, /silent
       if startype eq 'all' then good=where(xlist eq xlist)
       if startype eq 'A3CRGB' then good=where(teff-intteff lt 213)

       oplot, intage(good), feh(good), psym=sym(1), color=50 
   endif

    if age eq 'on' then begin
 
           openw,unit, 'ageout.txt', /get_lun
           printf, unit, $
		'KICID, 2mass, Mseis, Mmod, FeHapg, FeHmod, alphaapg, alphamod, loggseis, loggseiserr, loggspec, loggmod, Teffapg, Teffmodel, Teffmodmerr, Teffmodperr, evstateapg, Agemod, Agemerr, Ageperr, Agemaxmod'
           free_lun, unit
	   mavail=[0.6,0.7, 0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6]
           fehavail=[0.4,0.2,0.0,-0.2,-0.4,-0.6,-0.8,-1.0,-1.2,-1.4,-1.6,-1.8,-2.0]
	   alphaavail=[0.0,0.4]
           mbelow=['060','070','080', '090', '100', '110','120', '130','140', '150','160', '170', '180', '190', $
                   '200', '210','220', '230','240', '250','260']
           zhname=['p040', 'p020', 'm000', 'm020', 'm040', 'm060', 'm080', 'm100', 'm120', 'm140', 'm160', 'm180', 'm200']
           alname=['00','04']
           print, 'age1', n_elements(good)
           for i=0, n_elements(good)-1 do begin
               useless=min(abs(mavail-EMPSeisMscale(good(i))),ii)
               useless=min(abs(fehavail-CorFeH14(good(i))),k)
	       useless=min(abs(alphaavail-CoralphaFe14(good(i))),l)
;               if l eq 0 then trackname='/home/spitzer/tayar/EVOLUTION/output/YRECgrid/nodiff/nodiff_out/m'+mbelow(ii)+'zh'+$
;                    zhname(k)+'y264a13'+'_nodiff.track'
;               if l ne 0 then 
;               trackname='/home/spitzer/tayar/EVOLUTION/output/YRECgrid/nodiff/nodiff_out2/m'+mbelow(ii)+'zh'+zhname(k)+$
;                        'y264a23al'+alname[l]+'_nodiff.track'
               trackname='/home/jtayar/EVOLUTION/output/YRECgrid/nodiff/nodiff_out3/m'+mbelow(ii)+'zh'+zhname(k)+$
                        'y264a23al'+alname[l]+'_nodifv.track'
               result=file_test(trackname)
;               print, result
               if result eq 1 and OCTSeisMscale(good(i)) lt 2.8 then begin
                   print, trackname
          readcol2, trackname, skipline=14, $
                
                model, shells, Age, logLLsun, logRRsun, $
               logg, logTeff, Mccore, $
               Mcenv, R, T, Rho, P, $
               cappaenv, $
               ClogT, ClogRHO, $
               ClogP, CBETA, CETA,CX, $
               CY, CZ, LppI, LppII, LppIII, $
               LCNO, Ltripalpha, LHeC, Lgrav, $
               Lneutrinos, clSNU, GaSNU,$
               pp, pep, hep, Be7, B8, N13, O15, F17, diag1, diag2,$
               cHe3, cC12, cC13, cN14, cN15, cO16, cO17, cO18, $
               sHe3, sC12, sC13, sN14, sN15, sO16, sO17,sO18, sH2, $
               sLi6, sLi7, sBe9, $
               sX, sY, sZ, sZovrX, $
               Jtot, KErottot, totalI, CZI, Omsurface, Omcenter, $ ;FIXED 10/8
               Prot, Vrot, TauCZ, Hshellmfracbase, Hshellmfracmidpoint, $
               Hshellmfractop,  rfracbase, rfracmidpoint, rfractop, $
               Pphot, mass, $
           FORMAT="F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,D,D,D,D,D,F,F,F,F,F,F,F,F,F,F,F,F,F"
 ;this read in is wrong, se coreAngMomPlotter  ; FIXED 15/3/27
                   goodmod=where(CX lt 0.5)
                   if goodmod[0] ne -1 then begin
                   logggood=logg(goodmod)
                   logTeffgood=logTeff(goodmod)
                   Agegood=Age(goodmod)
;                   useless=min(abs(logggood-SeisLogg(good(i))),modindex)
;                   useless=min(abs(logggood-SeisLogg(good(i))-SeisLoggerr(good(i))),modindexm)        
;                   useless=min(abs(logggood-SeisLogg(good(i))+SeisLoggerr(good(i))),modindexp)  
                   useless=min(abs(logggood-OCTSeisLoggscale(good(i))),modindex)
                   useless=min(abs(logggood-OCTSeisLoggscale(good(i))-0.1),modindexm)        
                   useless=min(abs(logggood-OCTSeisLoggscale(good(i))+0.1),modindexp)             
           
                   loggmatch=logggood(modindex)
                   teffmatch=10^(logTeffgood(modindex))
                   teffmatchmerr=10^(logTeffgood(modindexm))-teffmatch
                   teffmatchperr=10^(logTeffgood(modindexp))-teffmatch
                   Agematch=Agegood(modindex)
                   Agemerr=Agegood(modindexm)-Agematch
                   Ageperr=Agegood(modindexp)-Agematch
                   Agemax=Agegood(n_elements(Agegood)-1)
                   openw,unit, 'ageout.txt', /get_lun, /append
                   print, teffmatch, SeisM(good(i))
                   printf, unit, KICID(good(i)),TMASSID(good(i)), OCTSeisMscale(good(i)), Mavail(ii), CorFeH13(good(i)), fehavail(k), CoralphaFe13(good(i)), alphaavail(l), $
                           OCTSeisLoggscale(good(i)), 0.1,CorLogg13(good(i)), loggmatch, CorTeff13(good(i)),teffmatch, teffmatchmerr, teffmatchperr, $
                           CEvstates(good(i)), agematch, Agemerr, Ageperr, agemax 
                   free_lun, unit
                   endif
               endif
           end
    endif
    
;if colortype eq 'alpha' then begin
;    al_legend, ['alpha= '+string(max(colorlist)), 'alpha= '+string(.5*(max(colorlist)+min(colorlist))),'alpha= '+string(min(colorlist))],$
;        psym=[8,8,8], $
;        colors=[255,127,0],/right_legend, /bottom_legend
;endif
;if colortype eq 'C' then begin
;     al_legend, ['C/Fe= '+string(max(colorlist)), 'C/Fe= '+string(.5*(max(colorlist)+min(colorlist))),'C/Fe= '+string(min(colorlist))],$
;        psym=[8,8,8], $
;        colors=[255,127,0],/right_legend, /bottom_legend
;endif
;if colortype eq 'N' then begin
;     al_legend, ['N/Fe= '+string(max(colorlist)), 'N/Fe= '+string(.5*(max(colorlist)+min(colorlist))),'N/Fe= '+string(min(colorlist))],$
;        psym=[8,8,8], $
;        colors=[255,127,0],/right_legend, /bottom_legend
;endif
    if rotkms eq 'on' then begin

    datadirectory= '~/EVOLUTION/output/ZRmod/'
    infiles=[ $               
               'm200w10e5Td6e4SbRGp2inFsbTjdT',  'm200w25e6Td1e3SbRGp2inFsbTjdT', 'm200w75e6Td1e3SbRGp2inFsbTjdT', $
               'm220w10e5Td6e4SbRGp2inFsbTjdT',  'm220w25e6Td1e3SbRGp2inFsbTjdT', 'm220w75e6Td1e3SbRGp2inFsbTjdT', $
               'm240w10e5Td6e4SbRGp2inFsbTjdT',  'm240w25e6Td1e3SbRGp2inFsbTjdT', 'm240w75e6Td1e3SbRGp2inFsbTjdT', $
               'm260w10e5Td6e4SbRGp2inFsbTjdT',  'm260w25e6Td1e3SbRGp2inFsbTjdT', 'm260w75e6Td1e3SbRGp2inFsbTjdT', $
               'm280w10e5Td6e4SbRGp2inFsbTjdT',  'm280w25e6Td1e3SbRGp2inFsbTjdT', 'm280w75e6Td1e3SbRGp2inFsbTjdT', $
               'm300w10e5Td6e4SbRGp2inFsbTjdT',  'm300w25e6Td1e3SbRGp2inFsbTjdT', 'm300w75e6Td1e3SbRGp2inFsbTjdT' $
            ]
    hlines=14 
    outfile= 'mZRcompare'
    outtext=''
    textperrun=[  'm20_w6' ,'M20_w1', 'M20_w3','m22_w6' ,'M22_w1', 'M22_w3', $
                  'm24_w6' ,'M24_w1', 'M24_w3', 'm26_w6' ,'M26_w1', 'M26_w3', $
                  'm28_w6' ,'M28_w1', 'M28_w3', 'm30_w6' ,'M30_w1', 'M30_w3']
    type='ZRmod'
    title='ZRmod'
    outdirectory='~/Documents/Apogee/APGgraphs/'
      GraphRotKmsModel, type, outdirectory, datadirectory, infiles, $
         outfile, outtext, textperrun, xmin, xmax, ymin, ymax, hlines, title
    cd, '~/Documents/Apogee/idl'
    endif


    if plttracks eq 'on' then begin
;          yalph=['y264a18al00']
          yalph=['y273a17al00']

;       if plttracksall eq 'on' then yalph=[yalph, 'y249a13al00','y279a13al00','y249a23al00','y279a23al00']
       if plttracksall eq 'on' then yalph=[yalph, 'y239a12al00','y290a12al00','y239a22al00','y290a22al00']
       if plttracksal eq 'on' then yalph=[yalph, 'y273a17al02', 'y273a17al04']  ;yalph=[yalph, 'y264a18al04']
       if plttracksbf eq 'on' then yalph=[yalph, 'y264a23al00']


;	   mavail=[0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6] ;diff
	   mavail=[0.6, 0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6] ;difv
           fehavail=[0.4,0.2,0.0,-0.2,-0.4,-0.6,-0.8,-1.0,-1.2,-1.4,-1.6,-1.8,-2.0]
	   alphaavail=[0.0,0.2,0.4]
;           mbelow=['080', '090', '100', '110','120', '130','140', '150','160', '170', '180', '190', $ ;diff
;                   '200', '210','220', '230','240', '250','260']
           mbelow=['060', '070','080', '090', '100', '110','120', '130','140', '150','160', '170', '180', '190', $ ;difv
                   '200', '210','220', '230','240', '250','260']
           zhname=['p040', 'p020', 'm000', 'm020', 'm040', 'm060', 'm080', 'm100', 'm120', 'm140', 'm160', 'm180', 'm200']
;           alname=['00','02','04']

          tmmin=min(A3PSeisMscale(good))
          tmmax=max(A3PSeisMscale(good))
          tzmin=min(CorFeH16(good))
          tzmax=max(CorFeH16(good))
;          tzmin=min(CorFeH14(good))
;          tzmax=max(CorFeH14(good))
;          tzmin=min(CorFeH12(good))
;          tzmax=max(CorFeH12(good))
          tms=[tmmin, tmmin, tmmax, tmmax]
          tzs=[tzmin, tzmax, tzmin, tzmax]
          if plttracksallm eq 'on' then begin
               tms=[mavail[0], mavail[0]]
               tzs=[tzmin, tzmax]
               for iii=1, n_elements(mavail)-1 do begin
                   if mavail[iii] gt tmmin and mavail[iii] lt tmmax then begin
                      tms=[tms, mavail[iii], mavail[iii]]
                      tzs=[tzs, tzmin, tzmax]
                   endif
               end
          endif


          plttypes=['Teff', 'logg', 'R', 'mass', 'FeH', 'alpha', 'L', 'CN']
          for j=0, n_elements(yalph)-1 do begin
           for i=0, n_elements(tms)-1 do begin
               useless=min(abs(mavail-tms(i)),ii)
               useless=min(abs(fehavail-tzs(i)),k)
;               trackname='/home/spitzer/tayar/EVOLUTION/output/YRECgrid/nodiff/nodiff_out2/m'+mbelow(ii)+'zh'+zhname(k)+$
;                        yalph(j)+'_nodiff.track'
;               trackname='/home/spitzer/tayar/EVOLUTION/output/YRECgrid/nodiff/nodiff_out3/m'+mbelow(ii)+'zh'+zhname(k)+$
;                        yalph(j)+'_nodifv.track'
               trackname='/home/jtayar/EVOLUTION/zuul/0/tayar/YRECgrid/nodiff_out4z/m'+mbelow(ii)+'fh'+zhname(k)+$
                        yalph(j)+'_grnodf.track'
print, trackname
               readcol2, trackname, /silent, skipline=14, $                
                model, shells, age, logLLsun, logRRsun, $
               logg, logTeff, Mccore, $
               Mcenv, R, T, Rho, P, $
               cappaenv, $
               ClogT, ClogRHO, $
               ClogP, CBETA, CETA,CX, $
               CY, CZ, LppI, LppII, LppIII, $
               LCNO, Ltripalpha, LHeC, Lgrav, $
               Lneutrinos, clSNU, GaSNU,$
               pp, pep, hep, Be7, B8, N13, O15, F17, diag1, diag2,$
               cHe3, cC12, cC13, cN14, cN15, cO16, cO17, cO18, $
               sHe3, sC12, sC13, sN14, sN15, sO16, sO17,sO18, sH2, $
               sLi6, sLi7, sBe9, $
               sX, sY, sZ, sZovrX, $
               Jtot, KErottot, totalI, CZI, Omenv, Omcore, $ ;FIXED 10/8
               Prot, Vrot, TauCZ, Hshellmfracbase, Hshellmfracmidpoint, $
               Hshellmfractop,  rfracbase, rfracmidpoint, rfractop, $
               Pphot, mass, $
           FORMAT="F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,D,D,D,D,D,F,F,F,F,F,F,F,F,F,F,F,F,F"


                

C12sun=0.00138155536
C13sun=0.0000155202526
N14sun=0.000409315048
N15sun=0.0
O16sun=0.00952570706
O17sun=0.0
O18sun=0.0000190982464
Xsun=0.718950491
bCH=ALOG10((sC12/12.+sC13/13.)/sX)-ALOG10((C12sun/12.+C13sun/13.)/Xsun)
bNH=ALOG10((sN14/14.+sN15/15.)/sX)-ALOG10((N14sun/14.+N15sun/15.)/Xsun)
CN=bCH-bNH
;print, 'waiting'
;wait, 100

           xptype=where(xplottype eq plttypes)
           yptype=where(yplottype eq plttypes)
           cptype=where(cplottype eq plttypes)

           if xptype eq '0' then xvars=10^logTeff
           if xptype eq '1' then xvars=logg
           if xptype eq '2' then xvars=10^logRRsun
           if xptype eq '6' then xvars=10^logLLsun
           if xptype eq '7' then xvars=CN
           if yptype eq '0' then yvars=10^logTeff
           if yptype eq '1' then yvars=logg
           if yptype eq '2' then yvars=10^logRRsun
           if yptype eq '6' then yvars=10^logLLsun
           if yptype eq '7' then yvars=CN
           if cptype eq '0' then cvars=10^logTeff
           if cptype eq '1' then cvars=logg
           if cptype eq '2' then cvars=10^logRRsun
           if cptype eq '3' then cvars=mavail(ii)
           if cptype eq '4' then cvars=fehavail(k)
           if cptype eq '5' and yalph(j) eq 'y273a17al04' then cvars=0.4
           if cptype eq '5' and yalph(j) eq 'y273a17al02' then cvars=0.2
           if cptype eq '5' and yalph(j) eq 'y273a17al00' then cvars=0.0
           if cptype eq '6' then cvars=10^logLLsun
           if cptype eq '7' then cvars=CN


           goodb=where(xvars gt xmin and xvars lt xmax and yvars gt ymin and yvars lt ymax and CX lt 0.5)
           if n_elements(goodb) gt 2 then begin

              colorz=(cvars-min(colorlist))*255/(max(colorlist)-min(colorlist))
              if goodcolors eq 'on' then  colorz=(cvars-min(colorlist(good)))* $
                     255/(max(colorlist(good))-min(colorlist(good)))
    
              if n_elements(colorz) eq 1 then colorz=xvars*0.0+colorz
              colorz=colorz(goodb)
  
              plots, xvars(goodb), yvars(goodb), psym=3, color= [colorz,colorz,colorz]
              xyouts, xvars(fix(goodb(n_elements(goodb)-.05*n_elements(goodb)))), yvars(fix(goodb(n_elements(goodb)-.05*n_elements(goodb)))), mbelow(ii)+'_'+zhname(k)+'_'+yalph(j), $
                   color=[colorz, colorz, colorz], charsize=0.6
           endif
           end
          end
    endif
    if plttracks eq 'rot' then begin
           yalph=['MrhFTFw1OvN'] ;'MrhFTFw1', 

;2_8MrhFTFw1yrec8
           mavail=[2.4, 2.6, 2.8, 3.0]
           mbelow=['2_4', '2_6', '2_8', '3_0']



;	   mavail=[0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6] ;diff
;	   mavail=[0.6, 0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6] ;difv
;           fehavail=[0.4,0.2,0.0,-0.2,-0.4,-0.6,-0.8,-1.0,-1.2,-1.4,-1.6,-1.8,-2.0]
;	   alphaavail=[0.0,0.4]
;           mbelow=['080', '090', '100', '110','120', '130','140', '150','160', '170', '180', '190', $ ;diff
;                   '200', '210','220', '230','240', '250','260']
;           mbelow=['060', '070','080', '090', '100', '110','120', '130','140', '150','160', '170', '180', '190', $ ;difv
;                   '200', '210','220', '230','240', '250','260']
;           zhname=['p040', 'p020', 'm000', 'm020', 'm040', 'm060', 'm080', 'm100', 'm120', 'm140', 'm160', 'm180', 'm200']
;           alname=['00','04']

          tmmin=min(EMPSeisMscale(good))
          tmmax=max(EMPSeisMscale(good))
;          tzmin=min(CorFeH(good))
;          tzmax=max(CorFeH(good))

          tms=[tmmin, tmmin, tmmax, tmmax]
;          tzs=[tzmin, tzmax, tzmin, tzmax]
          if plttracksallm eq 'on' then begin
               tms=[mavail[0], mavail[0]]
;               tzs=[tzmin, tzmax]
               for iii=1, n_elements(mavail)-1 do begin
 
                   tms=[tms, mavail[iii], mavail[iii]]
;                   tzs=[tzs, tzmin, tzmax]
               end
          endif


          plttypes=['Teff', 'logg', 'R', 'mass', 'FeH', 'alpha', 'L', 'CN']
          for j=0, n_elements(yalph)-1 do begin
           for i=0, n_elements(tms)-1 do begin
               useless=min(abs(mavail-tms(i)),ii)
;               useless=min(abs(fehavail-tzs(i)),k)
;               trackname='/home/spitzer/tayar/EVOLUTION/output/YRECgrid/nodiff/nodiff_out2/m'+mbelow(ii)+'zh'+zhname(k)+$
;                        yalph(j)+'_nodiff.track'
;               trackname='/home/spitzer/tayar/EVOLUTION/output/YRECgrid/nodiff/nodiff_out3/m'+mbelow(ii)+'zh'+zhname(k)+$
;                        yalph(j)+'_nodifv.track'
               trackname='/home/jtayar/EVOLUTION/output/RotGrid_rh/'+mbelow(ii)+yalph[j]+'.track'


               readcol2, trackname, /silent, skipline=14, $                
       		         model, shells, Age, logLLsun, logRRsun, $
           		        logg, logTeff, Mccore, $
		                Mcenv, R, T, Rho, P, $
		                cappaenv, $
		                ClogT, ClogRHO, $
		                ClogP, CBETA, CETA,CX, $
		                CY, CZ, LppI, LppII, LppIII, $
      	 	 	       LCNO, Ltripalpha, LHeC, Lgrav, $
 		               Lneutrinos, $
		               pp, pep, hep, Be7, B8, N13, O15, F17, $
   		               cHe3, cC12, cC13, cN14, cN15, cO16, cO17, cO18, $
      		               sHe3, sC12, sC13, sN14, sN15, sO16, sO17,sO18, sH2, $
        		       sLi6, sLi7, sBe9, $
       			        sX, sY, sZ, sZovrX, $
        		       Jtot, KErottot, totalI, CZI, Omsurface, Omcenter, $ ;FIXED
        		       Prot, Vrot, TauCZ, Hshellmfracbase, Hshellmfracmidpoint, $
       		               Hshellmfractop,  rfracbase, rfracmidpoint, rfractop, $
            		       Pphot, mass, $
             FORMAT="F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F "

C12sun=0.00138155536
C13sun=0.0000155202526
N14sun=0.000409315048
N15sun=0.0
O16sun=0.00952570706
O17sun=0.0
O18sun=0.0000190982464
Xsun=0.718950491
bCH=ALOG10((sC12/12.+sC13/13.)/sX)-ALOG10((C12sun/12.+C13sun/13.)/Xsun)
bNH=ALOG10((sN14/14.+sN15/15.)/sX)-ALOG10((N14sun/14.+N15sun/15.)/Xsun)
CN=bCH-bNH

           xptype=where(xplottype eq plttypes)
           yptype=where(yplottype eq plttypes)
           cptype=where(cplottype eq plttypes)

           if xptype eq '0' then xvars=10^logTeff
           if xptype eq '1' then xvars=logg
           if xptype eq '2' then xvars=10^logRRsun
           if xptype eq '6' then xvars=10^logLLsun
           if xptype eq '7' then xvars=CN
           if yptype eq '0' then yvars=10^logTeff
           if yptype eq '1' then yvars=logg
           if yptype eq '2' then yvars=10^logRRsun
           if yptype eq '6' then yvars=10^logLLsun
           if yptype eq '7' then yvars=CN
           if cptype eq '0' then cvars=10^logTeff
           if cptype eq '1' then cvars=logg
           if cptype eq '2' then cvars=10^logRRsun
           if cptype eq '3' then cvars=mavail(ii)
;           if cptype eq '4' then cvars=fehavail(k)
;           if cptype eq '5' and yalph(j) eq 'y264a18al04' then cvars=max(colorlist)
;           if cptype eq '5' and yalph(j) ne 'y264a18al04' then cvars=00

           goodb=where(xvars gt xmin and xvars lt xmax and yvars gt ymin and yvars lt ymax and CX lt 0.5)
           if n_elements(goodb) gt 2 then begin

              colorz=(cvars-min(colorlist))*255/(max(colorlist)-min(colorlist))
              if goodcolors eq 'on' then  colorz=(cvars-min(colorlist(good)))* $
                     255/(max(colorlist(good))-min(colorlist(good)))
    
              if n_elements(colorz) eq 1 then colorz=xvars*0.0+colorz
              colorz=colorz(goodb)
  
              plots, xvars(goodb), yvars(goodb), psym=3, color= [colorz,colorz,colorz]
              xyouts, xvars(fix(goodb(n_elements(goodb)-.05*n_elements(goodb)))), yvars(fix(goodb(n_elements(goodb)-.05*n_elements(goodb)))), mbelow(ii), $;+'_'+zhname(k)+'_'+yalph(j), $
                   color=[colorz, colorz, colorz], charsize=0.6
           endif
           end
          end
    endif


device, /close
;print, 'medloggerr', median(SeisLoggerr(good)), 'm', median(SeisMerr(good)), 'feh', median(CorFeHerr(good))
;print, 'meanloggerr', mean(SeisLoggerr(good)), 'm', mean(SeisMerr(good)), 'feh', mean(CorFeHerr(good))
loadct, 0

    if xhist eq 'on' or xhist eq 'log' then begin
	set_plot, 'ps'
	!p.font=0
	device, file=outdirectory+'hist'+outfile, bits=8,$
        	 /inches, /encapsulated, xsize=8, decomposed=0,$
        	 ysize=5, /portrait,/color, set_font='Helvetica'
print, n_elements(xlist(good))

         if xhist ne 'log' then plot,  [x1, x2], [0, n_elements(xlist(good))/3.], /nodata,   $
            CLIP=[x1,x2,0, n_elements(xlist(good))/3.], NOCLIP=0, $
            XTITLE=xstring, $
            YTITLE='Number', xstyle=1, $
            ystyle=1, xthick=5, ythick=5, $
            charsize=2.3, charthick=5, $
;           charsize=3, charthick=5, $
;         xticks=2, yticks=3, $
;         yminor=10,$ ;7 
         xrange=[x1, x2] 
         ;, $

        if xhist eq 'log' then plot,  [x1, x2], [0.9, n_elements(xlist(good))], /nodata,   $
            CLIP=[x1,x2,0, n_elements(xlist(good))], NOCLIP=0, $
            XTITLE=xstring, $
            YTITLE='Number', xstyle=1, $
            ystyle=1, xthick=5, ythick=5, /ylog,$
            charsize=2.3, charthick=5, $
;           charsize=3, charthick=5, $
;         xticks=2, yticks=3, $
;         yminor=10,$ ;7 
         xrange=[x1, x2] 

	plothist, xlist(good),xhist1, yhist, bin=xrang/20., /overplot
        meanhist=strcompress(string(mean(xlist(good))),/remove_all)
        stddhist=strcompress(string(stddev(xlist(good))), /remove_all)
        xyouts,  0.9*abs(x2-x1)+x1, 0.9*max(yhist), outtext, charsize=1, $
         charthick=2
        xyouts, 0.1*abs(x2-x1)+x1, 1.0*max(yhist), 'mean='+meanhist, charsize=1, $
         charthick=2
        xyouts, 0.5*abs(x2-x1)+x1, 1.0*max(yhist), 'sigma='+stddhist, charsize=1, $
         charthick=2
        device, /close
    endif


cd, '~/Documents/Apogee/idl'
set_plot, 'x'
end

