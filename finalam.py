import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from math import sqrt
def master(x_ar,a1=0,c1=0,s1=1,a2=0,c2=0,s2=1,a3=0,c3=0,s3=1,a4=0,c4=0,s4=1,a5=0,c5=0,s5=1,a6=0,c6=0,s6=1,a7=0,c7=0,s7=1,a8=0,c8=0,s8=1,a9=0,c9=0,s9=1,a10=0,c10=0,s10=1,a11=0,c11=0,s11=1):
    return [_1gaussian(i,a1,c1,s1)+_1gaussian(i,a2,c2,s2)+_1gaussian(i,a3,c3,s3)+_1gaussian(i,a4,c4,s4)+_1gaussian(i,a5,c5,s5)+_1gaussian(i,a6,c6,s6)+_1gaussian(i,a7,c7,s7)+_1gaussian(i,a8,c8,s8)+_1gaussian(i,a9,c9,s9)+_1gaussian(i,a10,c10,s10)+_1gaussian(i,a11,c11,s11) for i in x_ar]
def _1gaussian(x, a1, c1, s1):
    return (a1*(1/((s1)*(np.sqrt(2*np.pi))))*(np.exp(-((x-c1)**2/((2*(s1))**2)))))
def fwhm(x):
    return 2*sqrt(2*np.log(2))*(x)
def straightline(x,m,c):
    return x*m+c
def energy(x):
    return calfac*(x)+inter
def plot1g(x, var):
    for i in range(int(len(var)/3)):
        print(*list(var[i*3+j] for j in range(3)))
        plt.plot(x, _1gaussian(x, *list(var[i*3+j] for j in range(3))),'--',lw=0.6)
a=np.loadtxt('am1.txt')
b=np.loadtxt('background.txt')
a[:,1]=a[:,1]-b[:,1]
for i in range(len(a[:,1])):
    if(i<90) : a[i,1]=0

x_ar,y_ar,maxch,maxc,I=a[:2000,0],a[:2000,1],[],[],[]
plt.plot(x_ar,y_ar,'g')

for i in range(len(x_ar)):
    if(y_ar[i]==max(y_ar[i-min(1,i):i+1])):
        if(y_ar[i]>10) and (y_ar[i]-y_ar[i-3])>4 and y_ar[i]-y_ar[i+3]>4:
            maxc.append(y_ar[i])
            maxch.append(x_ar[i])
print("max. count values-",maxc)
print("max. count at channel values-",maxch)

calfac=(59.54-13.94)/(928.4730782104126-217.5)
inter=59.54-calfac*928.4730782104126
ee=-1*inter/calfac
print("calibration factor-",calfac,"intercept-",inter)
E=[]
for i in range(len(maxch)):
    E.append(calfac*maxch[i])
print("energy values according to peaks-",E)

x1=[94,181,211,277,320,406,560,925]
x2=[109,195,228,286,341,418,576,934]
C=[[100,102,104],[183,186,188,190,193],[218,220,223,0],[278,280,0],[325,327,335,0,0],[407,409,413,416],[562,566,568,571,0],[926,928,0]]
A=[]
channel,E,ER,FWHM,channelerr,FWHMerr,ERerr,Eerr=[],[],[],[],[],[],[],[]
for i in range(len(x1)):
    print("peak ",i+1)
    x_ar=a[int(x1[i]):int(x2[i]),0]
    y_ar=a[int(x1[i]):int(x2[i]),1]
    p0=[]
    bound=[],[]
    for k in range(len(C[i])):
        p0.append(max(y_ar))
        #if k==1: p0.append(max(y_ar)*100)
        #else : p0.append(1)
        if C[i][k]!=0 :
#           p0.append(max(y_ar))
           p0.append(C[i][k])
        else :
#           p0.append(1)
           p0.append(C[i][0])
        p0.append(0.2)
        bound[0].append(0)
        bound[1].append(10000000)
        if C[i][k]!=0 :
            bound[0].append(C[i][k]-0.5)
            bound[1].append(C[i][k]+0.5)
        else :
            bound[0].append(x_ar[0])
            bound[1].append(x_ar[-1])
        bound[0].append(0)
        bound[1].append(10)
    print(p0,bound)
    
#    bound=[0]*len(C[i])+[t-3 for t in C[i]]+[0]*len(C[i]),[10000]*len(C[i])+[t+3 for t in C[i]]+[10]*len(C[i])
    popt_2gauss,pcov_2gauss=scipy.optimize.curve_fit(master,x_ar,y_ar,p0,bounds=bound)
    plot1g(x_ar,popt_2gauss)
    plt.plot(x_ar,master(x_ar,*popt_2gauss))
    b=int((len(popt_2gauss))/3)
    maxchannel=0
    maxchannelerr=0
    for j in range(b):
        E.append(energy(popt_2gauss[3*j+1]))
#        if(popt_2gauss[3*j]==max(popt_2gauss[3*j] for j in range(b))):
        channel.append(popt_2gauss[3*j+1])
        channelerr.append(pcov_2gauss[3*j+1,3*j+1]**0.5)
        Eerr.append(energy(ee+pcov_2gauss[3*j+1,3*j+1]**0.5))
#        logE.append(np.log(energy(popt_2gauss[b+j])))
        FWHM.append(fwhm(popt_2gauss[3*j+2]))
        FWHMerr.append(fwhm(pcov_2gauss[3*j+2,3*j+2]**0.5))
        ER.append((fwhm(popt_2gauss[3*j+2])/(popt_2gauss[3*j+1]))*100)
#        logER.append(np.log(ER[j]))
        ERerr.append(ER[-1]*(sqrt((FWHMerr[-1]/FWHM[-1])**2+(Eerr[-1]/E[-1])**2)))
#        logERerr.append(np.log(1+(ERerr[j]/ER[j])))
        print("FWHM",j+1,"-",FWHM[-1],"(+/-)",FWHMerr[-1],"\nEnergy",j+1,"-",E[-1],"(+/-)",Eerr[-1],"\nEnergy resolution",j+1,"-",ER[-1],"(+/-)",ERerr[-1],'%')
#    channel.append(maxchannel)
#    channelerr.append(maxchannelerr)
plt.show()
f.open("newfe.txt",'w')
f.write(E,Eerr,FWHM)

''' for k in range(2):
        p0.append(max(y_ar))
        
        p0.append(1)
        bound[0].append(0)
        bound[0].append(x_ar[0])
        bound[0].append(0)
        bound[1].append(10000000)
        bound[1].append(x_ar[-1])
        if max(y_ar>100):
            bound[1].append(4)
        else:
            bound[1].append(2)
'''

'''

def _4gaussian(x, a1,a2,a3,a4,c1,c2,c3,c4,s1,s2,s3,s4):
    gauss=_1gaussian(x, a1, c1, s1)+_1gaussian(x, a2, c2, s2)+_1gaussian(x, a3, c3, s3)+_1gaussian(x, a4, c4, s4)
    return gauss
def _3gaussian(x, a1,a2,a3,c1,c2,c3,s1,s2,s3):
    gauss=_1gaussian(x, a1, c1, s1)+_1gaussian(x, a2, c2, s2)+_1gaussian(x, a3, c3, s3)
    return gauss
def _2gaussian(x, a1,a2,c1,c2,s1,s2):
    gauss=_1gaussian(x, a1, c1, s1)+_1gaussian(x, a2, c2, s2)
    return gauss

def _5gaussian(x, a1,a2=0,a3=0,a4=0,a5=0,a6=0,a7=0,a8=0,a9=0,c1,c2=0,c3=0,c4=0,c5=0,c6=0,c7=0,c8=0,c9=0,s1,s2=1,s3=1,s4=1,s5=1,s6=1,s7=1,s8=1,s9=1):
    gauss=_1gaussian(x, a1, c1, s1)+_1gaussian(x, a2, c2, s2)+_1gaussian(x, a3, c3, s3)+_1gaussian(x, a4, c4, s4)+_1gaussian(x, a5, c5, s5)
    return gauss

#x3=[2,5,4,1,1,1,1,4,3,1]
#p=[
#logER,,logE,logERerr
channel,E,ER,FWHM,channelerr,FWHMerr,ERerr,Eerr=[],[],[],[],[],[],[],[]
for i in range(len(x1)):
    print("peak ",i+1)
    x_ar=a[int(x1[i]):int(x2[i]),0]
    y_ar=a[int(x1[i]):int(x2[i]),1]
    if (x3[i]==1):
        p0=[max(y_ar),(x1[i]+x2[i])/2.0,1]
        popt_2gauss,pcov_2gauss=scipy.optimize.curve_fit(_1gaussian,x_ar,y_ar,p0,bounds=([0]*3,[10**10]*3))
        plt.plot(x_ar,_1gaussian(x_ar, *popt_2gauss),'--')
    elif (x3[i]==2):
        p0=[max(y_ar),max(y_ar),(2*x1[i]+x2[i])/3.0,(x1[i]+2*x2[i])/3.0,1,1]
        popt_2gauss,pcov_2gauss=scipy.optimize.curve_fit(_2gaussian,x_ar,y_ar,p0,bounds=([0]*6,[10**10]*6))
        plt.plot(x_ar,_2gaussian(x_ar, *popt_2gauss),'--')
    elif (x3[i]==3):
        p0=[max(y_ar)/2.0,max(y_ar),max(y_ar)/2.0,(3*x1[i]+x2[i])/4.0,(2*x1[i]+2*x2[i])/4.0,(x1[i]+3*x2[i])/4.0,1,1,1]
        popt_2gauss,pcov_2gauss=scipy.optimize.curve_fit(_3gaussian,x_ar,y_ar,p0,bounds=([0]*9,[10**10]*9))
        plt.plot(x_ar,_3gaussian(x_ar, *popt_2gauss),'--')
    elif (x3[i]==4):
        p0=[max(y_ar)/2.0,max(y_ar),max(y_ar),max(y_ar)/2.0,(4*x1[i]+x2[i])/5.0,(3*x1[i]+2*x2[i])/5.0,(2*x1[i]+3*x2[i])/5.0,(x1[i]+4*x2[i])/5.0,1,1,1,1]
        popt_2gauss,pcov_2gauss=scipy.optimize.curve_fit(_4gaussian,x_ar,y_ar,p0,bounds=([0]*12,[10**10]*12))
        plt.plot(x_ar,_4gaussian(x_ar, *popt_2gauss),'--')
    elif (x3[i]==5):
        p0=[max(y_ar)/2.0,max(y_ar),max(y_ar),max(y_ar),max(y_ar)/2.0,(5*x1[i]+x2[i])/6.0,(4*x1[i]+2*x2[i])/6.0,(3*x1[i]+3*x2[i])/6.0,(2*x1[i]+4*x2[i])/6.0,(x1[i]+5*x2[i])/6.0,1,1,1,1,1]
        popt_2gauss,pcov_2gauss=scipy.optimize.curve_fit(_5gaussian,x_ar,y_ar,p0,bounds=([0]*15,[10**10]*15))
        plt.plot(x_ar,_5gaussian(x_ar, *popt_2gauss),'--')
    plot1g(x_ar, popt_2gauss)
#    print(popt_2gauss)
    b=int((len(popt_2gauss))/3)
    maxchannel=0
    maxchannelerr=0
    for j in range(b):
        E.append(energy(popt_2gauss[3*j+1]))
#        if(popt_2gauss[3*j]==max(popt_2gauss[3*j] for j in range(b))):
        channel.append(popt_2gauss[b+j])
        channelerr.append(pcov_2gauss[b+j,b+j]**0.5)
        Eerr.append(energy(pcov_2gauss[b+j,b+j]**0.5))
#        logE.append(np.log(energy(popt_2gauss[b+j])))
        FWHM.append(fwhm(popt_2gauss[2*b+j]))
        FWHMerr.append(fwhm(pcov_2gauss[2*b+j,2*b+j]**0.5))
        ER.append((fwhm(popt_2gauss[2*b+j])/(popt_2gauss[b+j]))*100)
#        logER.append(np.log(ER[j]))
        ERerr.append(ER[-1]*(sqrt((FWHMerr[-1]/FWHM[-1])**2+(Eerr[-1]/E[-1])**2)))
#        logERerr.append(np.log(1+(ERerr[j]/ER[j])))
        print("FWHM",j+1,"-",FWHM[-1],"(+/-)",FWHMerr[-1],"\nEnergy",j+1,"-",E[-1],"(+/-)",Eerr[-1],"\nEnergy resolution",j+1,"-",ER[-1],"(+/-)",ERerr[-1],'%')
#    channel.append(maxchannel)
#    channelerr.append(maxchannelerr)
plt.plot(a[:,0],a[:,1],'k',lw=1.2)
plt.show()
plt.plot(E,FWHM,'bo')
plt.errorbar(E,FWHM,yerr=FWHMerr,xerr=Eerr,fmt="none")
plt.show()
plt.plot(E,ER,'bo')
plt.errorbar((E),(ER),yerr=(ERerr),fmt="none")
plt.show()
p0=[-0.5,0]
vaars,errorstline=scipy.optimize.curve_fit(straightline,np.log(E),np.log(ER),p0)
logERerr=list(np.log(1+ERerr[i]/ER[i]) for i in range(len(ER)))
print('parameters of straight line fitting-',vaars)
plt.errorbar(np.log(E),np.log(ER),yerr=logERerr,fmt="none")
plt.plot(np.log(E),straightline(np.log(E),*vaars),lw=0.5)
#plt.plot(np.log(E),straightline(np.log(E), -0.838,*vaars[1:2]),'k',lw=0.5)
plt.plot(np.log(E),np.log(ER),'bo')
a=vaars[0]
b=vaars[1]
c=errorstline[0][0]**0.5
d=errorstline[1][1]**0.5
plt.plot(np.log(E),straightline(np.log(E),a+c,b-d),lw=0.5)
plt.plot(np.log(E),straightline(np.log(E),a-c,b+d),lw=0.5)
plt.show()
"""
plt.plot(logE,straightline(logE,a,b),lw=0.5)
print("straight line variables : slope-",a,'(+/-)',c,' intercept-',b,'(+/-)',d)
ene=[8.047,9.571,22.162,24.942,61.288,88]
print(channel)
del(channel[5:7])
del(channelerr[5:7])
print(channel)
print(channelerr)
p0=[]
vaars,errorstline=scipy.optimize.curve_fit(straightline,channel,ene,p0)
plt.errorbar(channel,ene,xerr=channelerr,fmt="none")
plt.plot(logE,straightline(logE,a,b),lw=0.5)
plt.plot(channel,ene,'--')
plt.show()
"""
'''

