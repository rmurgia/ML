
;; Read in the output of SpecExtract.

pro read_spectra

base = '../'
filename = base+'spec2048_n5000_z2.000.dat'


ztime  = 0.0d
omega0 = 0.0d
omegaL = 0.0d
omegab = 0.0d
h100   = 0.0d
box100 = 0.0d
Xh     = 0.0d
nbins  = 0L
numlos = 0L

openr,1,filename

;; Header
readu,1,ztime,omega0,omegaL,omegab,h100,box100,Xh
readu,1,nbins,numlos


;; Sight-line locations (axis, x,y,z)
iaxis = lonarr(numlos)
xaxis = dblarr(numlos)
yaxis = dblarr(numlos)
zaxis = dblarr(numlos)
readu,1,iaxis,xaxis,yaxis,zaxis

;; Sight-line scale (comoving kpc/h, km/s)
posaxis = dblarr(nbins)
velaxis = dblarr(nbins)
readu,1,posaxis,velaxis

;; Normalised gas density, log(rho/<rho>)
density = dblarr(numlos*nbins)  
readu,1,density

;; HI data (HI/H, T in K, vpec in km/s, optical depth)
H1frac   = dblarr(numlos*nbins) 
temp_H1  = dblarr(numlos*nbins)
vpec_H1  = dblarr(numlos*nbins)
tau_H1  = dblarr(numlos*nbins)
readu,1,H1frac,temp_H1,vpec_H1,tau_H1

ind = where(finite(tau_H1) eq 0)
if(ind(0) ne -1) then begin
   print,'Bad elements in tau_H1!',n_elements(ind)
endif

close,1


print
print,filename
print,'z:      ',ztime
print,'nbins:  ',nbins
print,'numlos: ',numlos
print,'OmegaM: ',omega0
print,'OmegaL: ',omegaL
print,'Omegab: ',omegab
print,'Hubble: ',h100
print,'BoxSize:',box100
print,'Xh:     ',Xh
print,'<F>:    ',mean(exp(-tau_H1[*]))
print


;; Plot transmission in first sight-line
PLOTLOS = 12
window,0,xsize=600,ysize=400
Device,Retain=2
!p.font=-1
plot,velaxis,exp(-tau_H1[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]),xstyle=1,ystyle=1,charsize=1.75,yrange=[-0.1,1.1],ytitle='F',xtitle='vH [km/s]'

;; Plot temperature-density plane, 10000 random points
ind = floor(n_elements(density)*randomu(10, 10000,/uniform))
window,1,xsize=400,ysize=400
Device,Retain=2
!p.font=-1
plot,alog10(density(ind)),alog10(temp_H1(ind)),xstyle=1,ystyle=1,charsize=1.75,xtitle='log(1+delta)',ytitle='log(T/K)',psym=3,xrange=[-2,3],yrange=[3,6]


end
