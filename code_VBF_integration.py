#!/bin/env python

import collections
import matplotlib
import pyVBF
import math
import numpy as np
import pylab as pl
import sympy
from sympy import limit, Symbol, oo,Ellipse
from collections import deque
from matplotlib import patches
from astropy.units import Quantity 	#referenced from https://github.com/cta-observat					ory/ctapipe/blob/master/ctapipe/reco/hillas.pyL7					1 1
from astropy import units as u


__all__ = [
    'MomentParameters',
    'HighOrderMomentParameters',
    'hillas_parameters', ]

#################################################
# function to compute integration and return it
################################################
def computeAreaUnderCurve(t_array):
    # no need to calculate dx because we can make width of each rect
    #  1
    # using reimann sums  x.i * dx
    # using variable sumn
    # referenced from http://personal.bgsu.edu/~zirbel/5920/calculator/numerical_integration.pdf 
	sumn = 0.0
	l = 0
	while (l < (len(t_array) - 1 )) :
		sumn += ( ((t_array[l] + t_array[l+1]) / 2) * 1)
		l = l + 1 
    		# sumation should yield area under the curve 
        	##### now return sumn which now holds the area under the curve
	return sumn
####################################################
####################################################
# create array to hold 11,283*40 values for histogram
def integrationWindow(e_array):
	alen = len(e_array)
	new_array = [0] * 40 * alen
    # loops to create array of all elements and return it 
	for i in range(len(e_array)):
		for j in range(len(e_array[i])):
			new_array[i * 40 + j] = e_array[i][j]
		 
	return new_array
####################################################
####################################################
########## get lowest value in array
def lowestArray(array):
    #find lowest value in the array
	temp = array[0]
	index1 = 0
	for p in range(len(array)):
		if (array[p] < temp ):
			temp = array[p]
			index1 = p
	## return the lowest value and the index ti was found at 
	return temp,index1
######################################################
def highestArray(array):
    #find highest value in array
	temp = array[0]
	for p in range(len(array)):
		if (array[p] > temp ):
			temp = array[p]
			index = p
	# return the highest value and the index it was found at 
	return temp,index
######################################################
#####################################################
# Funtion to create histogram from all events in file
# CARE200010.vbf
def allEventsHist(fn):
	#	
	# Plot trace
	# Read a vbf simulation file and plot a trace for a single PM in a single        event
	# Open the file, print how many events it contains
	f = pyVBF.VBankFileReader(fn)
	print "Events in file in function :", f.numPackets()
	nn = f.numPackets()
    	## initialize loop variable
    	n = 0
    	## intialize hist to acculmuate event values
    	fig2 = pl.figure()
    
    	hist,binC = np.histogram([] , bins=100, range=[0,300],density=False)
  	
	while (n < 2 ): 
    		niceShower = n  # event ID
    		vp = f.readPacket(niceShower) # read the packet containing the event
        	
		### check to see if packet has an event 
       		if (vp.hasArrayEvent() ):
    			ae = vp.getArrayEvent() # get the array event
    			event = ae.getEvent(0)  # get the event for the first telescope (					the only one in the SCT files)
    			qi = integrationWindow(event)
			binB = np.linspace(0 ,500, 100)
    			hist1,bin1 = np.histogram(qi,bins=100,range=[0,300],density=False)
			hist = hist + hist1
		
    			n = n+1
		else :
			n = n+1
	# take hist array and find highest value which is the pedestal value
	hv,ind = highestArray( hist )
    	ax = fig2.add_subplot(111)
    	ax.bar(bin1[:-1] , hist, width=(bin1[1] - bin1[0]), color='blue')
	ax.set_xlim(0,300)
    	pl.yscale('log')
    	pl.show()
	return bin1[ind]
#####################################################
#####################################################






# # Plot trace
# 
# Read a vbf simulation file and plot a trace for a single PM in a single event

# Open the file, print how many events it contains
f = pyVBF.VBankFileReader("CARE200010.vbf")
#print "Events in file:", f.numPackets()

# An event with a nice shower
niceShower = 8819 # event ID
nns = int(raw_input("Enter Event NUmber based on number of events in package\n"))
notniceShower = nns
vp = f.readPacket(notniceShower) # read the packet containing the event
ae = vp.getArrayEvent() # get the array event
if vp.getArrayEvent():
	print "Has EVENT\n"

event = ae.getEvent(0) # get the event for the first telescope (the only one in the SCT files)


gps = event.getGPSTime()
nicePM = 8130 # a PM with signal

t = event[notniceShower]
s = [x for x in range(len(t))]


############################################################
### Troubleshooting Area
## get cut value
valy = "Y"
inp = str(raw_input( "Print histogram, enter Y or N\n"))
if inp == valy:
	name = "CARE200010.vbf"
	cut = allEventsHist(name)
	print cut 
	print "*****************************"
############################################################
### start analysis
# get pedestal, average first 10 values of t and store in to variable
############################################################
ped = 0
for i in range(0,10):
	ped += t[i]
ped = ped / 10

#the area under the pedestal is 
areaped = ped * len(t)
#############################################################
# now that we have pedestal and the area under the curve subtract the pedestal from the area under the curve to aquire the area we want

auc = computeAreaUnderCurve(t)
wa = auc - areaped
print "The desired area under the curve is " , wa

## wa or wanted area now holds the desired area 
## untested until i configure bashr

theinput = str(raw_input("Do you want to print the charge of pixel\n"))
if theinput == "Y":
	fig1 = pl.figure()

	ax = fig1.add_subplot(111)
	ax.step(s,t)     #check
	ax.set_xlabel("Time [samples]")
	ax.set_ylabel("Signal [ADC]")
	ax.set_ylim(0,50)

	pl.step(s, t)
	pl.xlabel("Time [samples]")
	pl.ylabel("Signal [ADC]")

	#pl.ylim(0,50)
	#pl.plot([0,40], [ped,ped], 'k-',lw=1)
	#pl.text(0,16,"Pedestal",color='red')
	pl.show()

MomentParameters = collections.namedtuple(
    "MomentParameters",
    "size,cen_x,cen_y,length,width,r,phi,psi,miss")

#####################m###################################
######## Clean pixels
#### create a cut value for charges in order to produce clean image
########################################################

# #### Basic functions to read camera geo and draw

def getCameraGeo(filename = "camgeo.dat"):
    camgeo = np.genfromtxt(filename)
    xpix = camgeo[:,1]
    ypix = camgeo[:,2]
    idpix = camgeo[:,0]

    return (xpix, ypix)

#def drawCamera(ax, charges, cmap=pl.get_cmap("plasma")):
#	(xpix, ypix) = getCameraGeo()
#	ax.scatter(xpix, ypix, s = 30, c=charges, lw=0, cmap=cmap) # marker='s'

# modified draw function 
def drawCamera(ax, charges, cmap=pl.get_cmap("plasma")):
	(xpix, ypix) = getCameraGeo()
	xpix1 = [0] * len(xpix)
	ypix1 = [0] * len(ypix)
        charges1 = [0] * len(charges)
	# boolean that holds whether pixel has a charge
	there = [0] * len(charges)
	# create variable to hold a set distance for a pixel distance cut
	distance = (xpix[1] - xpix[0])**2
	distance1 = (ypix[1] - ypix[0])**2
	#d = math.sqrt(distance + distance1)
	d = 10
	print math.sqrt(distance + distance1)
	w = 0
	ccval = int(raw_input("Enter Charge Value\n"))
        while (w < len(charges)):
		# apply charge cut to eliminate background noise
		# charge cut of about 1200 leaves behind the brightest pixels and the shower
        	if (charges[w]  >= ccval):
                                #print charges[w]
                                #print "charge ^"
  			xpix1[w] = xpix[w]
                	ypix1[w] = ypix[w]
                	charges1[w] = charges[w]
			there[w] = 1
                	w = w+1
                else:
			there[w] = 0
                	w = w+1
	if sum(there) >= 5:
	# create cut based on distance and if a pixel is by itself
	### if a pixel is by itself delete it 
        ##  if distance < x && distance > x
        #   and if distance < y && distance > y
        #   set charge1 val to 0 and xpix1 && ypix1 to zero
		grid = np.array(zip(xpix,ypix,charges1) , dtype={'names' :['x','y','c'], 'formats' :['f4' , 'f4' , 'f4']})
	
		aa = [0] * len(xpix1)
		x = 0
		while (x < len(charges1) ):

			cut = (grid['x'] < xpix1[x] + d) & (grid['x'] > xpix1[x] - d) & (grid['y'] < ypix1[x] + d) & (grid['y'] > ypix1[x] - d)
			aa[x] = np.sum(grid[cut]['c'] > 0)
			if (aa[x] <= 1 ):
				charges1[x] = 0
				xpix1[x] = 0
                        	ypix1[x] = 0
				there[x] = 0

			x = x + 1
		# test 
		print"---------------------------------"
		print "Sum of booleans is " ,sum(there)
		print "--------------------------------"
		low, lv = lowestArray(charges1)
		print "lowest values in charges", low,lv
		print "--------------------------------"
	# hp[0] = size , hp[1] = x center ,  hp[2] = y ceneter, 
	# hp[3] = length , hp[4] = width, hp[5] = r , hp[6] = phi
	# hp[7] =  psi , hp[8] = miss 
	# Structure to hold hillas parameters ^

	# direction function to calculate hillas parameter based on info of ellipse
		hp = direction(charges1,xpix1,ypix1)
		print hp
		print "##########################################################"

	# use sypy ellipse to create an ellipse
		ee = Ellipse(xpix1,ypix1)

	# test the ellipse to see if it works properly 
		print "ellipse center is ",ee.center

	# construct line to represent direction and hillas parameters 
		xc = hp.cen_x		# center of shower 
		yc = hp.cen_y
		l = hp.length 
		wid = hp.width 				# length of shower
	# phi is angle between center of image and centroid
		angphi = (hp.phi * 180 ) / math.pi			# angle of the shower in radians
		angpsi = ( hp.psi * 180 ) / math.pi 			# angle in degrees
	
	# attempting to compute points on line 
		xf = xc + (l * math.cos(angphi))	
		yf = yc	+ (l * math.sin(angphi))
		xf1 = xc - (l * math.cos(angphi)) 
        	yf1 = yc - (l * math.sin(angphi))

	# test points
		txf = xc + (l * math.cos(angpsi))
        	tyf = yc + (l * math.sin(angpsi))
		tx1f = xc - (l * math.cos(angpsi))
		ty1f = yc - (l * math.sin(angpsi))

		print "********************"
		print xf,yf
		print xf1,yf1
		print "*********************"
		print "test points "
		print txf,tyf
		print tx1f,ty1f
		print '*******************'
		xi = (2 * xc) - xf
		yi = (2 * yc) - yf
	
		# test point
		txi = (2 * xc) - txf
        	tyi = (2 * yc) - tyf


		print "********************"
        	print xi,yi
		print "*******************"
		print "test intial points "
		print txi,tyi
		# change this part !!!!!!! look for predefined function to fit ellipse to data
		# psipy - determine over pixels what to include 
		bright = max(charges1)
		b,ind = highestArray(charges1)
		bx = xpix[ind]
		by = ypix[ind]

		print "brightest val is " , bright
		print "brightest pixel is ", ind,b
		print "the angle phi is ",hp[6],angphi
		print " angle of psi is ", angpsi
		# changing between xpix , ypix and xpix1 , ypix1 shows shower in whole or piece 
		ax.scatter(xpix, ypix, s = 30, c=charges1, lw=0, cmap=cmap) # marker='
		# attempting to draw line to signify direction
		#ax.plot([xi,xf] , [yi,yf],color='k',linewidth=2)
		#ax.plot([txi,txf] , [tyi,tyf],color='k',linewidth=2)
	
		ax.text(xc,yc, r'\C')
		ax.text(bx,by, r'B')
		#ellipse = patches.Ellipse((xc,yc),wid,l,angpsi,lw=2,fill=False)
		e1 = patches.Ellipse((xc,yc),wid*2,l*2,-angphi,lw=2,fill=False)
		# test stuff
	
		angtheta = 180 - 90 - (-angphi)
		angzeta = 180 - 90 - (angpsi)
	
		# calculate side inbetween phi and 90 degrees
		sideA = (yc / math.sin(angphi) ) * math.sin(angphi)
		pp = sideA
		# calculate midpoint of image axis and draw line to centroid from that point
		mdx = (xc + pp ) / 2
		mdy = (0 + yc) / 2
		# calculate missing angle in triangle with phi in it 
		phi1 = 180 - 90 - (-angphi)
		# calucate angle inbetween phi1 and psi , since phi1 and psi are basically equall
		zeta = 90 - phi1
		# calculate theta or angle at center 
		eta = 180 - 90 - zeta
		# calculate alpha
		alpha = 180 - 90 - angpsi
		print "phi1, zeta,eta, alpha "
		print phi1,zeta,eta,alpha
		print "****************************"
		print angtheta,angzeta,sideA
		#ax.plot([-50,250],[0,0],color='r',ls='-',linewidth=2)
		ax.plot([0,xc],[0,yc],color='r',linewidth=2)
		#ax.plot([xc,xc],[0,yc],color='r',linewidth=2)
		#ax.plot([pp,xc],[0,yc],color='r',linewidth=2)
		#ax.plot([pp,xc],[0,0],color='r',linewidth=2)
		#ax.plot([0,mdx],[0,mdy],color='r',linewidth=2)
		e2 = patches.Ellipse((xc,yc),wid,l,eta,lw=2,fill=False)
		e3 = patches.Ellipse((xc,yc),wid,l,alpha,lw=2,fill=False)
		ax.add_patch(e2)
		#ax.add_patch(e3)

		#ax.add_patch(ellipse)
		# return the charges1 array to be used in getting hilas parameters 
		# in direction function
		return charges1
	return sum(there)
################################################################################################################################################################

###################################################################
###############################################################
#######################################
def direction(charges1,xpix1,ypix1):

	(xpix, ypix) = getCameraGeo()
	#unit = u.mm				# get units 

	pix_x = Quantity(np.asanyarray(xpix1,dtype=np.float64)).value
	pix_y = Quantity(np.asanyarray(ypix1,dtype=np.float64)).value
	image  = np.asanyarray(charges1, dtype=np.float64)
	assert pix_x.shape == image.shape
	assert pix_y.shape == image.shape

	size = image.sum()
	mdata = np.row_stack([pix_x,pix_y,pix_x * pix_x,pix_y * pix_y,pix_x * pix_y]) * image
	ms = mdata.sum(axis=1) / size

	vx2 = ms[2] - ms[0] **2
	vy2 = ms[3] - ms[1] **2
	vxy = ms[4] - ms[0] * ms[1]

	dd = vy2 - vx2
	zz = np.sqrt(dd **2 + 4.0 * vxy **2)

	# miss
	uu = 1.0 + dd/zz
	vv = 2.0 - uu
	miss = np.sqrt(((uu * ms[0] **2 + vv * ms[1] **2) / 2.0) - (ms[0] * ms[1] * 2.0 * vxy / zz))

	# shower shape parameters
	width = np.sqrt(vx2 + vy2 - zz)
	length = np.sqrt(vx2 + vy2 + zz)
	azwidth = np.sqrt(ms[2] + ms[3] - zz)

	# rotation angle of ellipse relative to centroid
	tanpsi_numer = (dd + zz) * ms[1] + 2.0 * vxy * ms[0]
	tanpsi_denom = (2 * vxy * ms[1]) - (dd - zz) * ms[0]
	psi = ((np.pi / 2.0) + np.arctan2(tanpsi_numer,tanpsi_denom))
	# polar coordinates of centroid
	rr = np.hypot(ms[0],ms[1])
	phi = np.arctan2(ms[1],ms[0])

	return MomentParameters(size=size,cen_x=ms[0],cen_y=ms[1],length=length,width=width,r=rr,phi=phi,psi=psi,miss=miss)

	
###################################################################################################################################################


sums = [np.sum(event[i]) for i in range(len(event))]
fig = pl.figure(figsize=(12,12))
ax2 = fig.add_subplot(111)    
shower = drawCamera(ax2,sums)
pl.show()




