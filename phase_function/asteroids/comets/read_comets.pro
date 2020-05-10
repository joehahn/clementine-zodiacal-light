pro read_comets,name,q,e,i,O,w,DATA=data

;Read the Minor Planet Center's comet datafile containing comet
;orbit elements. This file is obtained at the URL
;http://cfa-www.harvard.edu/iau/Ephemerides/Comets/Soft00Cmt.txt

;Outputs are comet number, name, perihelion in AU,
;eccentricity, inclination, longitude of ascending node, and
;argument of perihelion. All angles are in degrees. Option
;output is DATA, which is a string array containing the contents
;of the input file. The optional input parameter is a string
;naming the input file, which defaults to FILE='comets.dat'

;First, get orbits for non-returning, long-period comets
;set input filename
file='Cat2.dat'
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
;comet name
name=strtrim(strmid(data,0,15),2)
;perihelion in AU
q=float(strtrim(strmid(data,39,8),2))
;eccentricity
e=float(strtrim(strmid(data,49,8),2))
;argument of perihelion in degrees
w=float(strtrim(strmid(data,65,8),2))
;longitude of ascending node in degrees
O=float(strtrim(strmid(data,75,8),2))
;inclination in degrees
i=float(strtrim(strmid(data,85,8),2))

;Next, get orbits for returning, periodic comets
;set input filename
file='Cat3.dat'
print,'reading '+file
;Get number of entries in the datafile.
spawn_str='wc -l '+file
spawn,spawn_str,wc_result
N2=long(wc_result)
N2=N2(0)
;read data
data=strarr(N2)
get_lun,unit
openr,unit,file
readf,unit,data
close,1
;comet name
name2=strtrim(strmid(data,0,15),2)
;perihelion time
time2=strtrim(strmid(data,20,10),2)
;perihelion in AU
q2=float(strtrim(strmid(data,39,8),2))
;eccentricity
e2=float(strtrim(strmid(data,49,8),2))
;argument of perihelion in degrees
w2=float(strtrim(strmid(data,65,8),2))
;longitude of ascending node in degrees
O2=float(strtrim(strmid(data,75,8),2))
;inclination in degrees
i2=float(strtrim(strmid(data,85,8),2))

;Next, select the most recent apparitions of the 
;periodic comets
;another name
another_name2=strtrim(strmid(data,142,28),2)
j=where(another_name2 ne '')-1
jj=where(j gt 0)
if (jj(0) ne -1) then j=j(jj)
name3=name2(j)
q3=q2(j)
e3=e2(j)
w3=w2(j)
O3=O2(j)
i3=i2(j)

;combine orbit elements
name=[name,name3]
q=[q,q3]
e=[e,e3]
w=[w,w3]
O=[O,O3]
i=[i,i3]

end
