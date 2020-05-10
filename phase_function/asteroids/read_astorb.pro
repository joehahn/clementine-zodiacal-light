pro read_astorb,number,name,H,diameter,a,e,i,O,w,M, $
  DATA=data,FILE=file

;Read Ted Bowell's `astorb.dat' file containing asteroid orbit
;elements, diameters, etc. Bowell's file is obtained at the URL
;ftp://ftp.lowell.edu/pub/elgb/astorb.html .

;Outputs are asteroid number, name, absolute magnitude,
;diameter in km, semimajor axis in AU, eccentricity, inclination,
;longitude of ascending node, argument of perihelion, and mean
;anomaly. All angles are in degrees. Option output is DATA,
;which is a string array containing the contents of the input
;file. The optional input parameter is a string naming
;Bowell's input file, which defaults to FILE='astorb.dat'

;set input filename
if (keyword_set(file) eq 0) then file='astorb.dat'
print,'reading '+file

;Get number of entries in the datafile.
spawn_str='wc -l '+file
spawn,spawn_str,wc_result
N=long(wc_result)
N=N(0)

;read data
data=strarr(N)
get_lun,unit
openr,unit,file
readf,unit,data
close,1

;asteroid number
number=long(strtrim(strmid(data,0,5),2))

;asteroid name
name=strtrim(strmid(data,6,18),2)

;absolute magnitude
H=float(strtrim(strmid(data,41,5),2))

;diameter in km
diameter=float(strtrim(strmid(data,58,5),2))

;mean anomaly in degrees
M=float(strtrim(strmid(data,114,10),2))

;argument of perihelion in degrees
w=float(strtrim(strmid(data,125,10),2))

;longitude of ascending node in degrees
O=float(strtrim(strmid(data,136,10),2))

;inclination in degrees
i=float(strtrim(strmid(data,146,10),2))

;eccentricity in degrees
e=float(strtrim(strmid(data,157,10),2))

;semimajor axis in AU
a=float(strtrim(strmid(data,168,12),2))

end
