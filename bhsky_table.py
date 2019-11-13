""" Calculates (b, phi, Amp) table for passing photons from the sky
"""

import sys
import math
import numpy as np

def out(rtop, b, emax, emin, delnew, s, eobs):
    sum = 0.0
    if (rtop < Rad): rtop = Rad
    if (b == 0.0): b = bmin / 1000.0

    de = (emax - emin) / float(iter)
    e = emin - (de / 2.0)

    for j in range(0, iter-1):
        e = e + de

        if e < math.pi / 4.0:
            top = math.tan(e) * de
            bot = ((s ** 2) / ((b ** 2) * math.cos(e) ** 2)) - 1 + (Rs * math.cos(e) / s)
        else:
            top = de
            bot = ((s ** 2) / ((b ** 2) * (math.sin(e) ** 2))) - (math.cos(e) ** 2 / math.sin(e) ** 2) + (Rs / s) * (
                        (math.cos(e) ** 3) / (math.sin(e) ** 2))
        bot = math.sqrt(bot)
        term = top / bot

        if delnew > math.pi - delesc:
            sum = sum + term
            continue

        if e < eobs and delnew < math.pi / 2.0:
            sum = sum + 2.0 * term
            continue

        if e > eobs:
            sum = sum + term
    return sum


def comp(deltest, b, btest, delold, phideg, mphi):
	if Rad >= 1.50 * Rs:  # Star outside photon sphere
		
		if deltest > (math.pi - delesc):  # Infinity to observer
			b = btest
			delnew = deltest
			s = Rad
			eobs = math.acos(Rad / Dobs)
			emin = eobs
			emax = math.pi / 2.0
			sum = out(0, b, emax, emin, delnew, s, eobs)

		if delearth < deltest <= (math.pi - delesc):
			b = btest
			delnew = deltest
			rtop = 2.00 * b * math.cos(math.acos(-math.sqrt(27) * (Rs / 2) / b) / 3) / math.sqrt(3.00)
			s = rtop
			eobs = math.acos(s / Dobs)
			emin = 0.0
			emax = math.pi / 2.0
			sum = out(rtop, b, emax, emin, delnew, s, eobs)

		if deltest <= delearth:
			if phideg > 360.00:
				return sum, delnew, delold, b
			if mphi < eps:
				return sum, delnew, delold, b

			conv = 1.2
			delnew = (delold * delearth ** (conv - 1)) ** (1 / conv)
			b = Bbig * math.sin(delnew)
			rtop = 2 * b * math.cos(math.acos(-math.sqrt(27) * (Rs / 2) / b) / 3) / math.sqrt(3)
			s = rtop
			eobs = math.acos(s / Dobs)
			emin = 0
			emax = math.pi / 2
			converge = 'on'
			sum = out(rtop, b, emax, emin, delnew, s, eobs)

	elif Rad < 1.5 * Rs:
		delnew = 0.0
		sum = 0.0
		if Dobs > 1.5 * Rs:

			if deltest > (math.pi - delesc):
				if phideg > 360:
					return sum, delnew, delold, b
				b = btest
				delnew = deltest
				s = Rad
				eobs = math.acos(s / Dobs)
				emin = eobs
				emax = math.pi / 2.00
				sum = out(0, b, emax, emin, delnew, s, eobs)

			if delearth <= deltest < (math.pi - delesc):
				if phideg > 360:
					return sum, delnew, delold, b
				b = btest
				delnew = deltest
				rtop = 2.0 * b / math.sqrt(3.0) * math.cos(math.acos(-math.sqrt(27) * (Rs / 2) / b) / 3)
				s = rtop
				eobs = math.acos(s / Dobs)
				emin = 0
				emax = math.pi / 2
				sum = out(rtop, b, emax, emin, delnew, s, eobs)

			if deltest < delearth:
				if phideg > 360.0:
					return sum, delnew, delold, b
				if mphi < eps:
					return sum, delnew, delold, b
				conv = 1.010
				delnew = (delold * delearth ** (conv - 1.0)) ** (1.0 / conv)
				b = Bbig * math.sin(delnew)
				rtop = 2 * b / math.sqrt(3) * math.cos(math.acos(-math.sqrt(27) * (Rs / 2) / b) / 3)
				s = rtop
				eobs = math.acos(s / Dobs)
				emin = 0
				emax = math.pi / 2.0
				converge = 'on'
				sum = out(rtop, b, emax, emin, delnew, s, eobs)

		if Dobs <= 1.50 * Rs:
			if deltest > delesc and phideg < 300.0:
				b = btest
				delnew = deltest
				s = Rad
				eobs = math.acos(s / Dobs)
				emin = eobs
				emax = math.pi / 2
				sum = out(0, b, emax, emin, delnew, s, eobs)

			if deltest < delesc or phideg >= 300.0:
				if phideg > 360.0:
					return sum, delnew, delold, b
				if mphi < eps:
					return sum, delnew, delold, b
				conv = 1.010
				delnew = (delold * delesc ** (conv - 1)) ** (1 / conv)
				b = Bbig * math.sin(delnew)
				s = Rad
				eobs = math.acos(s / Dobs)
				emin = eobs
				emax = math.pi / 2.0
				converge = 'on'
				sum = out(0, b, emax, emin, delnew, s, eobs)

	return sum, delnew, delold, b


iter = 100000
eps = sys.float_info.epsilon
f1 = open("star.dat", "r")
f2 = open("star.out", "w+")

star_data = f1.readlines()

test = (str(star_data[1]))
test1 = (str(star_data[2]))
test2 = (str(star_data[3]))

Rs = (np.float(test[20:29]))
Rad = (np.float(test1[20:29]))
Dobs = (np.float(test2[20:29]))

f1.close()

Rinf = Rad / math.sqrt(1.00 - (Rs / Rad))
bmin = Rs * 1.500 * math.sqrt(3.00)
Bbig = Dobs / math.sqrt(1.00 - (Rs / Dobs))

sinarg = 1.500 * math.sqrt(3.00) * Rs * math.sqrt(1.00 - (Rs / Dobs)) / Dobs

if sinarg > 1.00:
    sinarg = 1.00

delesc = math.asin(sinarg)

if Rad > 1.500 * Rs:
    Earg = (Rad / Dobs) * math.sqrt((1.00 - (Rs / Dobs)) / (1.00 - (Rs / Rad)))
else:
    Earg = 1.500 * (Rs / Dobs) * math.sqrt((1.00 - (Rs / Dobs)) / (1.00 - (Rs / (1.500 * Rs))))

delearth = math.asin(Earg)

if Dobs < 1.500 * Rs:
    delesc = math.pi - delesc
    delearth = math.pi - delearth

converge = 'off'

mphi = math.pi / 180.0
delmax = math.pi
delmin = delesc
ddelbig = (delmin - delmax) / 179.47543232470
deltest = delmax - ddelbig

phideg = 0.0
delold = 0.0
b = 0.0
phinew = 0.0



for i in range(0, 1800):
	if phideg > 400:
		break
	if i <= 1: ddel = ddelbig  
	else: ddel = mdel

	if mphi >= 5.0 * math.pi / 180.0:
		ddel = ddel / 2.0
	deltest = deltest + ddel
	btest = Bbig * math.sin(deltest)
	if b == 0:
		b = bmin / 1000.0

	sum, delnew, delold, b = comp(deltest, b, btest, delold, phideg, mphi)

	phiold = phinew
	phinew = sum

	phideg = phinew * 180.0 / math.pi
	phiout = 180.0 - phideg
	mphi = phinew - phiold

	delout = delnew * 180.0 / math.pi
	mdel = (delnew - delold)
	delold = delnew

	if delnew * 180.0 / math.pi >= 179.9990:
		Amp = 1.0
		top = mdel * math.sin(delnew)
		bot = mphi * math.sin(math.pi - phinew)
		if (bot == 0):
			break
		if phideg == 0 and delout == 0: continue
		f2.write("{:.5f}".format(phideg) + "   " + "{:.5f}".format(delout) + "    " + "{:.5f}".format(Amp) + "\n")
		#print("{:.5f}".format(phideg) + "   " + "{:.5f}".format(delout) + "    " + "{:.5f}".format(Amp))
		if converge == 'on':
			sum, delnew, delold, b = comp(deltest, b, btest, delold, phideg, mphi)
		continue

	top = mdel * math.sin(delnew)
	bot = mphi * math.sin(math.pi - phinew)
	blueshift = 1.0
	if(bot==0):
		break
	if phideg == 0 and delout == 0: continue
	Amp = abs(top / bot) * blueshift

	if phiold == 0: Amp = 1.0
	if Amp >= 100: Amp = 100.0
	f2.write("{:.5f}".format(phideg) + "   " + "{:.5f}".format(delout) + "    " + "{:.5f}".format(Amp) + "\n")
	#print("{:.5f}".format(phideg) + "   " + "{:.5f}".format(delout) + "    " + "{:.5f}".format(Amp))
	if converge == 'on':
		sum, delnew, delold, b = comp(deltest, b, btest, delold, phideg, mphi)

f2.close()
