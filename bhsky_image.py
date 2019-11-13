import math
import numpy as np
import matplotlib.pyplot as plt

def interpA(theta, Aout, nth):
    for i in range(0, nth):
        if thint[i] < theta:
            continue
        term = (theta - thint[i-1])/(thint[i] - thint[i-1])
        Aout = Aint[i-1] + term*(Aint[i] - Aint[i-1])
        return Aout
    Aout = 0.0
    return Aout


def interpD(theta, delout, nth):
    for i in range(0, nth):
        if thint[i] < theta:
            continue
        term = (theta - thint[i - 1]) / (thint[i] - thint[i - 1])
        delout = delint[i - 1] + term * (delint[i] - delint[i - 1])
        return delout
    delout = -1.0
    return delout


def jump(th1, npt2, del2, Amp2, xin, yin, nth):
    th2 = 360.0 - th1
    del2 = interpD(th2, del2, nth)
    if del2 < earth:
        return x1, y1, A1, x2, y2, A2, npt2
    Amp2 = interpA(th2, Amp2, nth)
    apixel = math.sqrt(psize*Amp2)
    if apixel < amagcut and smag > smagcut:
        return x1, y1, A1, x2, y2, A2, npt2

    thin = del2*math.pi/180.0
    phin = math.atan2(yin, xin)
    phin = phin + cphiRot

    xin_orig = math.sin(thin)*math.cos(phin)
    yin_orig = math.sin(thin)*math.sin(phin)
    zin_orig = math.cos(thin)

    xin_orig = - xin_orig
    yin_orig = - yin_orig

    xin_view = xin_orig
    yin_view = yin_orig * math.cos(cthetaRot) - zin_orig * math.sin(cthetaRot)
    zin_view = yin_orig * math.sin(cthetaRot) + zin_orig * math.cos(cthetaRot)
    if zin_view < 0.0:
        return x1, y1, A1, x2, y2, A2, npt2

    r2in_view = math.sqrt(xin_view ** 2 + yin_view ** 2)
    thin_view = math.acos(zin_view)
    phin_view = math.atan2(yin_view, xin_view)

    thcomp = thin_view * 180.0 / math.pi
    if thin_view > 70.0:
        return x1, y1, A1, x2, y2, A2, npt2

    npt2 = npt2 + 1
    x2.insert(npt2, math.sqrt(1.0 - math.cos(thin_view)) * (xin_view / r2in_view) / camnorm)
    y2.insert(npt2, math.sqrt(1.0 - math.cos(thin_view)) * (yin_view / r2in_view) / camnorm)

    A2.insert(npt2, apixel)

    if abs(x2[npt2]) > 1.0 or abs(y2[npt2]) > 1.0:
        npt2 = npt2 - 1
        return x1, y1, A1, x2, y2, A2, npt2

    return x1, y1, A1, x2, y2, A2, npt2


f1 = open("star.dat", "r")
f2 = open("star.out", "r")
f3 = open("bhsky1.out", "w+")
f4 = open("bhsky2.out", "w+")

x1 = []
y1 = []
x2 = []
y2 = []

A1 = []
A2 = []

del2 = 0
Amp2 = 0

star_data = f1.readlines()
f1.close()

test = [[]]
for i in range(0, 9):
    test.insert(i, str(star_data[i+1]))

Rs = (np.float(test[0][20:29]))
Rad = (np.float(test[1][20:29]))
Dobs = (np.float(test[2][20:29]))
smagcut = (np.float(test[3][20:29]))
cthetaRot = (np.float(test[4][20:29]))
cphiRot = (np.float(test[5][20:29]))
sthetaRot = (np.float(test[6][20:29]))
sphiRot = (np.float(test[7][20:29]))

name = sphiRot
#print(" Rs        =", Rs)
#print(" Rad       =", Rad)
#print(" Dobs      =", Dobs)
#print(" smagcut   =", smagcut)
#print(" cthetaRot =", cthetaRot)
print(" cphiRot   =", cphiRot)
#print(" sthetaRot =", sthetaRot)
#print(" sphiRot   =", sphiRot)
name = cphiRot
#sphiRot = -sphiRot + 80.0

f = math.pi/180.00
zblue = math.sqrt(1.0 - (Rs/Dobs))
cthetaRot = cthetaRot*math.pi/180.00
cphiRot = cphiRot*math.pi/180.00
sthetaRot = sthetaRot*math.pi/180.00

#print('zblue =', zblue)

thint = np.loadtxt("star.out", skiprows=0, ndmin=2)[:, 0]
delint = np.loadtxt("star.out", skiprows=0, ndmin=2)[:, 1]
Aint = np.loadtxt("star.out", skiprows=0, ndmin=2)[:, 2]
f2.close()
nth = len(thint)
#print("Reading in star data..")

aread = np.loadtxt("new.num")[:, 0]
dread = np.loadtxt("new.num")[:, 1]
sread = np.loadtxt("new.num")[:, 2]

nptot = len(aread)

npt1 = -1
npt2 = -1

camnorm = math.sqrt(1.00 - math.cos(math.pi/4.00))
amagcut = 10.00**((-smagcut/2.5120)+2.00)
Rearthapp = Rad*math.sqrt((1.00-(Rs/Dobs))/(1.00 - (Rs/Rad)))
earth = math.atan2(Rearthapp, Dobs)*180.00/math.pi

#print("Earth angle subtended =", earth)

del1 = 0.0
Amp1 = 0.0

for i in range(0, nptot):
    alpha = aread[i] + sphiRot
    delta = dread[i]
    smag = sread[i]
    
    xin0 = -math.sin(delta*f+(math.pi/2.0))*math.cos(alpha*f+(math.pi/2.0))
    zin0 = math.sin(delta * f + (math.pi / 2.0)) * math.sin(alpha * f + (math.pi / 2.0))
    yin0 = -math.cos(delta*f+(math.pi/2.0))

    xin = xin0
    yin = yin0*math.cos(sthetaRot) - zin0*math.sin(sthetaRot)
    zin = yin0*math.sin(sthetaRot) + zin0*math.cos(sthetaRot)

    psize = 10.0**(-smag/2.5120 + 2.0)
    r3 = math.sqrt(xin**2 + yin**2 + zin**2)
    r2 = math.sqrt(xin**2 + yin**2)
    th1 = math.acos(zin/r3)*180.0/math.pi

    del1 = interpD(th1, del1, nth)
    if del1 < earth:
        x1, y1, A1, x2, y2, A2, npt2 = jump(th1, npt2, del2, Amp2, xin, yin, nth)
        #print('yo oy')
        continue

    Amp1 = interpA(th1, Amp1, nth)
    apixel = math.sqrt(psize*Amp1)
    if apixel < amagcut:
        x1, y1, A1, x2, y2, A2, npt2 = jump(th1, npt2, del2, Amp2, xin, yin, nth)
        #print('yo oy2')
        continue
        
    
    thin = del1*math.pi/180.0
    phin = math.atan2(yin, xin)
    phin = phin + cphiRot

    xin_orig = math.sin(thin)*math.cos(phin)
    yin_orig = math.sin(thin)*math.sin(phin)
    zin_orig = math.cos(thin)

    xin_view = xin_orig
    yin_view = yin_orig*math.cos(cthetaRot) - zin_orig*math.sin(cthetaRot)
    zin_view = yin_orig * math.sin(cthetaRot) + zin_orig * math.cos(cthetaRot)
    if zin_view < 0.0:
        x1, y1, A1, x2, y2, A2, npt2 = jump(th1, npt2, del2, Amp2, xin, yin, nth)
        #print('yo yo')
        continue

    r2in_view = math.sqrt(xin_view**2 + yin_view**2)
    thin_view = math.acos(zin_view)
    phin_view = math.atan2(yin_view, xin_view)

    thcomp = thin_view*180.0/math.pi
    if thcomp > 70.0:
        x1, y1, A1, x2, y2, A2, npt2 = jump(th1, npt2, del2, Amp2, xin, yin, nth)
        #print('yo')
        continue
    npt1 = npt1 + 1
    #print ('incoming2')
    x1.insert(npt1, math.sqrt(1.0 - math.cos(thin_view))*(xin_view/r2in_view)/camnorm)
    y1.insert(npt1, math.sqrt(1.0 - math.cos(thin_view)) * (yin_view / r2in_view) / camnorm)

    if abs(x1[npt1]) > 1.0 or abs(y1[npt1]) > 1.0:
        npt1 = npt1 - 1
        x1, y1, A1, x2, y2, A2, npt2 = jump(th1, npt2, del2, Amp2, xin, yin, nth)
        continue

    A1.insert(npt1, apixel)

    x1, y1, A1, x2, y2, A2, npt2 = jump(th1, npt2, del2, Amp2, xin, yin, nth)
    continue
print ('npt1 =', npt1)

X1 = []
Y1 = []
s1 = []

X2 = []
Y2 = []
s2 = []

for i in range(0, npt1 + 1):
    f3.write("{:.5f}".format(x1[i]) + "    " + "{:.5f}".format(y1[i]) + "    " + "{:.5f}".format(A1[i]) + "\n")
    X1.append(x1[i])
    Y1.append(y1[i])
    s1.append(float(A1[i]*A1[i]))

f3.close()

if npt2 == 0:
    npt2 = 2
    for i in range(0, 1):
        x2[i] = -2.0
        y2[i] = -2.0
        A2[i] = 0.0

#print("Writing out star data...")
print('npt2 =', npt2)

for k in range(0, npt2+1):
	f4.write("{:.5f}".format(x2[k]) + '    ' + "{:.5f}".format(y2[k]) + '    ' + "{:.5f}".format(A2[k]) + '\n')
	X2.append(x2[k])
	Y2.append(y2[k])
	s2.append(float(A2[k]*A2[k]))

f4.close()
r = Dobs/Rs

fig = plt.figure(figsize = [16, 12], dpi = 100)
ax = fig.add_subplot(111, autoscale_on = False, autoscalex_on = False, autoscaley_on = False)
ax.set_title('r/Rs = {:.2f}'.format(r))

newx = [-1, -1, 1, 1]
newy = [-1, 1, -1, 1]

ax.scatter(newx, newy, c= 'black', marker='o')
ax.scatter(X1, Y1, s=s1, c= 'white', marker='o')
ax.scatter(X2, Y2, s=s2, c= 'white', marker='o')
plt.autoscale(False)

ax.set_facecolor('xkcd:black')
ax.axis('scaled')
plt.savefig('insps {:.2f}.png'.format(name), dpi = 100)
#plt.show()
