#!/bin/env python

import collections
import matplotlib
import pyVBF
import math
import numpy as np
import pylab as pl
import sympy
from sympy import limit, Symbol, oo
from collections import deque
from matplotlib import patches
from astropy.units import Quantity  # referenced from https://github.com/cta-observatory/ctapipe/blob/master/ctapipe/reco/hillas.pyL7					1 1
from astropy import units as u
from astropy.coordinates import Angle


__all__ = [
    'MomentParameters',
    'HighOrderMomentParameters',
    'hillas_parameters', ]

MomentParameters = collections.namedtuple(
    "MomentParameters",
    "size,cen_x,cen_y,length,width,r,phi,psi,miss")

class CTA:
    def __init__(self, pfile, camfile,tarray):
        self.pfile = pfile
        self.camfile = camfile
        self.tarray = tarray
	# initalize variables for the for loop 
	a = 0
	b = 0
        f = pyVBF.VBankFileReader(self.pfile)
	# store total number of packets just in case you need it 
	self.totalEvents = int(f.numPackets())
	self.events = [0] * self.totalEvents
	print self.totalEvents
	# create a dictionary to store events that exist 
	self.dictt = {}
	for a in range(0,self.totalEvents):
       		vp = f.readPacket(a)
                if vp.hasArrayEvent():
			self.dictt[b] = a
                	#x = {b:a}
			b  += 1
	self.dicttlen = len(self.dictt)
    def computeAreaUnderCurve(self):
        # no need to calculate dx because we can make width of each rect
        #  1
        # using reimann sums  x.i * dx
        # using variable sumn
        # referenced from http://personal.bgsu.edu/~zirbel/5920/calculator/numerical_integration.pdf
        sumn = 0.0
        l = 0
        while l < (len(self.tarray) - 1):
            sumn += (((self.tarray[l] + self.tarray[l + 1]) / 2) * 1)
            l += 1
            # summation should yield area under the curve
            # now return sumn which now holds the area under the curve
        return sumn
    # create array to hold 11,283*40 values for histogram
    def findLowest(self,pix,charge):
	# using the pix and charge discern what the smallest value
	# in the array is
	temp = pix[0] 
	index = 0	
	for a in range(len(charge)):
		if charge[a] > 0 and pix[a] != 0 :
			if pix[a] < temp:
				temp = pix[a]
				#print "value of low",temp
				index = a
	return temp,index	
    def findHighest(self,pix,charge):
        # using the pix and charge discern what the smallest value
        # in the array is
        temp = -10
        index = 0
        for a in range(len(charge)):
                if charge[a] > 0 and pix[a] != 0:
                        if pix[a] > temp:
                                temp = pix[a]
				print "value of highest",temp
                                index = a
        return temp,index
    def integrationWindow(self,earray):
        alen = len(earray)
        new_array = [0] * 40 * alen
        # loops to create array of all elements and return it
        for i in range(len(earray)):
            for j in range(len(earray[i])):
                new_array[i * 40 + j] = earray[i][j]
        return new_array
    def lowestArray(self,array):
        # find lowest value in the array
        temp = array[0]
        index1 = 0
        for p in range(len(array)):
            if (array[p] < temp ):
                temp = array[p]
                index1 = p
        # return the lowest value and the index ti was found at
        return temp,index1
    def highestArray(self,array):
        # find highest value in array
        temp = array[0]
        index = 0
        for p in range(len(array)):
            if array[p] > temp :
                temp = array[p]
                index = p
        # return the highest value and the index it was found at
        return temp,index
    def allEventsHist(self):
        # Plot trace
        # Read a vbf simulation file and plot a trace for a single PM in a single event
        # Open the file, print how many events it contains
        f = pyVBF.VBankFileReader(self.file1)
        print "Events in file in function :", f.numPackets()
        nn = f.numPackets()
        ## initialize loop variable
        n = 0
        ## intialize hist to acculmuate event values
        fig2 = pl.figure()

        hist, binC = np.histogram([], bins=100, range=[0, 300], density=False)

        while (n < 2):
            niceShower = n  # event ID
            vp = f.readPacket(niceShower)  # read the packet containing the event

            ### check to see if packet has an event
            if vp.hasArrayEvent():
                ae = vp.getArrayEvent()  # get the array event
                event = ae.getEvent(0)  # get the event for the first telescope (					the only one in the SCT files)
                qi = self.integrationWindow(event)
                binB = np.linspace(0, 500, 100)
                hist1, bin1 = np.histogram(qi, bins=100, range=[0, 300], density=False)
                hist += hist1
                n += 1
            else:
                n += 1
        # take hist array and find highest value which is the pedestal value
        hv, ind = self.highestArray(hist)
        ax = fig2.add_subplot(111)
        ax.bar(bin1[:-1], hist, width=(bin1[1] - bin1[0]), color='blue')
        ax.set_xlim(0, 300)
        pl.yscale('log')
        pl.show()
        return bin1[ind]
    def getCameraGeo(self ):
        #  filename="camgeo.dat"
        self.camfile = np.genfromtxt(filename)
        xpix = self.camfile[:, 1]
        ypix = self.camfile[:, 2]
        idpix = self.camfile[:, 0]
        return xpix, ypix
    def hillasparameters(self,charges,x,y):
	unit = Quantity(x).unit
	x = Quantity(np.asanyarray(x,dtype=np.float64)).value
	y = Quantity(np.asanyarray(y,dtype=np.float64)).value
	charges = np.asanyarray(charges,dtype=np.float64)
	assert x.shape == charges.shape
        assert y.shape == charges.shape
	# sum of all charges which gives size
        size = charges.sum()
	# find average x and y values of coordinates
	meanx = np.sum(charges * x) / size
	meany = np.sum(charges * y) / size
	# find major axis value 
	# y = a*x + b
	sxx = np.sum(charges * (x - meanx) **2) / size
	syy = np.sum(charges * (y - meany) **2) / size
	sxy = np.sum(charges * (x - meanx) * (y - meany)) / size
	sxxx = np.sum(charges * (x - meanx) **3) / size
	syyy = np.sum(charges * (y - meany) **3) / size
	sxyy = np.sum(charges * (x - meanx) * (y - meany) **2) / size
	sxxy = np.sum(charges * (x - meanx) **2 * (y - meany)) / size
	sx4 = np.sum(charges * (x - meanx) **4 ) / size
	sy4 = np.sum(charges * (y - meany) **4) / size

	d0 = syy - sxx
	d1 = 2 * sxy
	
	# val = d^2 + 4 * sxy^2 
	d2 = d0 + np.sqrt(d0 * d0 + d1 * d1)
	a = d2 / d1
	# angle between ellipse major axis and x axis of camera
	delta = np.arctan(a)
	b = meany - a * meanx
	
	# used for higher order parameters
	cos_delta = 1 / np.sqrt(1 + a **2)
	sin_delta = a * cos_delta
	
	# compute hillas 
	w2 = (syy + a * a * sxx - 2 * a * sxy) / (1 + a **2)
	l2 = (sxx + a * a * syy + 2 * a * sxy) / ( 1 + a **2)
	
	if w2 < 0 :
		width = 0
	else:
		width = np.sqrt(w2)
	if l2 < 0 :
		length = 0
	else:
		length = np.sqrt(l2)
	# calculate miss 
	miss = np.abs(b / np.sqrt( 1 + a **2))
	r = np.sqrt(meanx * meanx + meany * meany)
	phi = np.arctan2(meany,meanx)

	# compute higher order moments ? 
		
	return MomentParameters(size=size, cen_x=meanx, cen_y=meany, length=length, width=width, r=r, phi=phi,psi=delta,miss=miss)

    def direction(self, charges1, xpix1, ypix1):
        (xpix, ypix) = self.getCameraGeo()
        # unit = u.mm				# get units

        pix_x = Quantity(np.asanyarray(xpix1, dtype=np.float64)).value
        pix_y = Quantity(np.asanyarray(ypix1, dtype=np.float64)).value
        image = np.asanyarray(charges1, dtype=np.float64)
        assert pix_x.shape == image.shape
        assert pix_y.shape == image.shape

        size = image.sum()
        mdata = np.row_stack([pix_x, pix_y, pix_x * pix_x, pix_y * pix_y, pix_x * pix_y]) * image
        ms = mdata.sum(axis=1) / size

        vx2 = ms[2] - ms[0] ** 2
        vy2 = ms[3] - ms[1] ** 2
        vxy = ms[4] - ms[0] * ms[1]

        dd = vy2 - vx2
        zz = np.sqrt(dd ** 2 + 4.0 * vxy ** 2)

        # miss
        uu = 1.0 + dd / zz
        vv = 2.0 - uu
        miss = np.sqrt(((uu * ms[0] ** 2 + vv * ms[1] ** 2) / 2.0) - (ms[0] * ms[1] * 2.0 * vxy / zz))

        # shower shape parameters
        width = np.sqrt(vx2 + vy2 - zz)
        length = np.sqrt(vx2 + vy2 + zz)
        azwidth = np.sqrt(ms[2] + ms[3] - zz)

        # rotation angle of ellipse relative to centroid
        tanpsi_numer = (dd + zz) * ms[1] + 2.0 * vxy * ms[0]
        tanpsi_denom = (2 * vxy * ms[1]) - (dd - zz) * ms[0]
        psi = ((np.pi / 2.0) + np.arctan2(tanpsi_numer, tanpsi_denom))
        # polar coordinates of centroid
        rr = np.hypot(ms[0], ms[1])
        phi = np.arctan2(ms[1], ms[0])

        return MomentParameters(size=size, cen_x=ms[0], cen_y=ms[1], length=length, width=width, r=rr, phi=phi, psi=psi,miss=miss)
    def drawCamera(self, x, charges, cmap=pl.get_cmap("plasma")):
        (xpix, ypix) = self.getCameraGeo()
        xpix1 = [0] * len(xpix)
        ypix1 = [0] * len(ypix)
        charges1 = [0] * len(charges)
        # create variable to hold a set distance for a pixel distance cut
        distance = (xpix[1] - xpix[0]) ** 2
        distance1 = (ypix[1] - ypix[0]) ** 2
        # d = math.sqrt(distance + distance1)
        d = 10
        print math.sqrt(distance + distance1)
        w = 0
	# boolean that holds whether pixel has a charge
	self.boolCharges = [0] * len(charges)
	chargeCut = int(raw_input("Please Enter a Charge Cut FAM"))
        while (w < len(charges)):
            # apply charge cut to eliminate background noise
            # charge cut of about 1200 leaves behind the brightest pixels and the shower
            if (charges[w] >= chargeCut):
                xpix1[w] = xpix[w]
                ypix1[w] = ypix[w]
                charges1[w] = charges[w]
                self.boolCharges[w] = 1
                w += 1
            else:
                self.boolCharges[w] = 0
                w += 1
        if sum(self.boolCharges) >= 5:
            # create cut based on distance and if a pixel is by itself
            ### if a pixel is by itself delete it
            ##  if distance < x && distance > x
            #   and if distance < y && distance > y
            #   set charge1 val to 0 and xpix1 && ypix1 to zero
            grid = np.array(zip(xpix, ypix, charges1), dtype={'names': ['x', 'y', 'c'], 'formats': ['f4', 'f4', 'f4']})

            aa = [0] * len(xpix1)
            x = 0
            while x < len(charges1):
                cut = (grid['x'] < xpix1[x] + d) & (grid['x'] > xpix1[x] - d) & (grid['y'] < ypix1[x] + d) & (
                    grid['y'] > ypix1[x] - d)
                aa[x] = np.sum(grid[cut]['c'] > 0)
                if (aa[x] <= 1):
                    charges1[x] = 0
                    xpix1[x] = 0
                    ypix1[x] = 0
                    self.boolCharges[x] = 0
                x += 1
            # test
            print"---------------------------------"
            print "Sum of booleans is ", sum(self.boolCharges)
            print "--------------------------------"
            low, lv = self.lowestArray(charges1)
            print "lowest values in charges", low, lv
            print "--------------------------------"
            # hp[0] = size , hp[1] = x center ,  hp[2] = y ceneter,
            # hp[3] = length , hp[4] = width, hp[5] = r , hp[6] = phi
            # hp[7] =  psi , hp[8] = miss
            hp = self.hillasparameters(charges1, xpix1, ypix1)
            print hp
            print "##########################################################"
            # construct line to represent direction and hillas parameters
            xc = hp.cen_x  # center of shower
            yc = hp.cen_y
            l = hp.length
            wid = hp.width  # length of shower

            # phi is angle between center of image and centroid
            # phi is angle between center of image and centroid
            angphi = (hp.phi * 180) / math.pi  # angle of the shower in radians
            angpsi = (hp.psi * 180) / math.pi  # angle in degrees

	    # find y values for limits of ellipse
	    ylow,i1 = self.lowestArray(ypix1)
	    yhigh,i2 = self.findHighest(ypix1,charges1)
	    xylow = xpix[i1]
	    xyhigh = xpix[i2]

	    # find x values for limits of ellipse
	    xhigh,i3 = self.highestArray(xpix1)
	    xlow,i4 = self.findLowest(xpix1,charges1)
	    yxhigh = ypix1[i3]
	    yxlow = ypix1[i4]
	    print "Values of y high and low\n",ylow,i1,yhigh,i2
	    print "Values of x high and low \n",xhigh,i3,xlow,i4
	    print "y high and min \n",max(ypix1),min(ypix1)

            # attempting to compute points on line
            xf = xc + (l * math.cos(angphi))
            yf = yc + (l * math.sin(angphi))
            xf1 = xc - (l * math.cos(angphi))
            yf1 = yc - (l * math.sin(angphi))
            # attempting to compute points on line
            xf = xc + (l * math.cos(angphi))
            yf = yc + (l * math.sin(angphi))
            xf1 = xc - (l * math.cos(angphi))
            yf1 = yc - (l * math.sin(angphi))

            # test points
            txf = xc + (l * math.cos(angphi))
            tyf = yc + (l * math.sin(angphi))
            tx1f = xc - (l * math.cos(angphi))
            ty1f = yc - (l * math.sin(angphi))
	    print "angle phi and psi\n" , angphi,angpsi
            print "********************"
            print xf, yf
            print xf1, yf1
            print "*********************"
            print "test points "
            print txf, tyf
            print tx1f, ty1f
            print '*******************'
            xi = (2 * xc) - xf
            yi = (2 * yc) - yf

            # test point
            txi = (2 * xc) - txf
            tyi = (2 * yc) - tyf

            print "********************"
            print xi, yi
            print "*******************"
            print "test intial points "
            print txi, tyi
            # change this part !!!!!!! look for predefined function to fit ellipse to data
            # psipy - determine over pixels what to include
            bright = max(charges1)
            b, ind = self.highestArray(charges1)
            bx = xpix[ind]
            by = ypix[ind]

            print "brightest val is ", bright
            print "brightest pixel is ", ind, b
            print "the angle phi is ", hp[6] 
            print " angle of psi is ", angpsi
		
	    fig = pl.figure(figsize=(12, 12))
	    ax = fig.add_subplot(111)
            # changing between xpix , ypix and xpix1 , ypix1 shows shower in whole or piece
            ax.scatter(xpix, ypix, s=30, c=charges1, lw=0, cmap=cmap)  # marker='
            # attempting to draw line to signify direction
            # ax.plot([xi,xf] , [yi,yf],color='k',linewidth=2)
            # ax.plot([txi,txf] , [tyi,tyf],color='k',linewidth=2)

            ax.text(xc, yc, r'\C')
            ax.text(bx, by, r'B')
            # ellipse = patches.Ellipse((xc,yc),wid,l,angpsi,lw=2,fill=False)
            e1 = patches.Ellipse((xc, yc), wid, l, angphi+180, lw=2, fill=False)
            # test stuff

            angtheta = 180 - 90 - (angphi)
            angzeta = 180 - 90 - (angphi)

            # calculate side inbetween phi and 90 degrees
            sideA = (yc / math.sin(angphi)) * math.sin(angphi)
            pp = sideA
            # calculate midpoint of image axis and draw line to centroid from that point
            mdx = (xc + pp) / 2
            mdy = (0 + yc) / 2
            # calculate missing angle in triangle with phi in it
            phi1 = 180 - 90 - (angphi)
            # calucate angle inbetween phi1 and psi , since phi1 and psi are basically equall
            zeta = 90 - phi1
            # calculate theta or angle at center
            eta = 180 - 90 - zeta
            # calculate alpha
            alpha = 180 - 90 - angphi
            print "phi1, zeta,eta, alpha "
            print phi1, zeta, eta, alpha
            print "****************************"
            print angtheta, angzeta, sideA
            # ax.plot([-50,250],[0,0],color='r',ls='-',linewidth=2)
            ax.plot([0, xc], [0, yc], color='r', linewidth=2)
            # ax.plot([xc,xc],[0,yc],color='r',linewidth=2)
            # ax.plot([pp,xc],[0,yc],color='r',linewidth=2)
            # ax.plot([pp,xc],[0,0],color='r',linewidth=2)
            # ax.plot([0,mdx],[0,mdy],color='r',linewidth=2)
            e2 = patches.Ellipse((xc, yc), l*10, wid*5, eta, lw=2, fill=False)
            e3 = patches.Ellipse((xc, yc), wid, l, alpha, lw=2, fill=False)
            ax.add_patch(e2)
            ax.add_patch(e3)

            # ax.add_patch(ellipse)
            # return the charges1 array to be used in getting hilas parameters
            # in direction function
            return charges1
        return 0

    def searchEvent(self,value):
	first = 0
	last = len(self.dictt) - 1
	found = False
	while first <= last and not found:
		mp = int((first + last)/2)
		if self.dictt[mp] == value:
			found = True
		else :
			if value < self.dictt[mp]:
				last = mp - 1
			else :
				first = mp + 1
	return found
			
		
    # Main Starts here
    # Plot trace
    # Read a vbf simulation file and plot a trace for a single PM in a single event
    # Open the file, print how many events it contains

f = pyVBF.VBankFileReader("CARE200010.vbf")

print "Events in file:", f.numPackets()
    # An event with a nice shower & not nice shower
#niceShower = 8819  # event ID
#notniceShower = 1500
#vp = f.readPacket(notniceShower)  # read the packet containing the event
#ae = vp.getArrayEvent()  # get the array event

#if vp.getArrayEvent(): 
#	print "Has EVENT"
#event = ae.getEvent(0)  # get the event for the first telescope (the only one in the SCT files)

#gps = event.getGPSTime()
#nicePM = 8130  # a PM with signal

#t = event[notniceShower]
#s = [x for x in range(len(t))]
    # instatiuate class
file1 = "CARE200010.vbf"
filename = "camgeo.dat"
t = []
c = CTA(file1, filename,t)
# check to see if user wants to run a certain run 
ev = str(raw_input("Do you wish to run a certain run ?\n"))
if (ev.capitalize() == "Y"):
	ew = int(raw_input("Please enter event number\n"))
	p = c.searchEvent(ew)
	if p == True :
		vp = f.readPacket(ew)
		ae = vp.getArrayEvent()
		event = ae.getEvent(0)
		gps = event.getGPSTime
		sums = [np.sum(event[i]) for i in range(len(event))]
		fig = pl.figure(figsize=(12, 12))
		ax2 = fig.add_subplot(111)
		shower = c.drawCamera(ax2, sums)
		pl.show()
	else :
		print "File does not exist in packet" 
else:
# else cycle through all the events 
	a = 0
	while a < f.numPackets():
		p = c.searchEvent(a)
        	if p == True :
			vp = f.readPacket(a)
	                ae = vp.getArrayEvent()
	                event = ae.getEvent(0)
	                gps = event.getGPSTime
			sums = [np.sum(event[i]) for i in range(len(event))]
			fig = pl.figure(figsize=(12, 12))
			ax2 = fig.add_subplot(111)
			shower = c.drawCamera(ax2, sums)
			pl.show()
			a += 1
		else:
			a += 1 
