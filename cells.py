from numba import jit
import numpy as np
from random import randint as rd
import matplotlib.pyplot as plt
import time
import pickle
#from numba.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
import time
import threading
 
#warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
#warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)
ydash=0
listorders=[]
for i in range(1000000):
    listorders.append(0)
class cell:
   m=1
   v=0
   h=0
   w=0
   e=100
   x=0
   y=0
   z=-5
   xr=0
   yr=0
   zr=0
   messageid=0
   messagekey=0
   message=0 
   settingm=1
   stateofaply=0
   mutationpercentage=10
   dna=np.array([[38,0.6],[[17,0,3,0],[1,0.05]],[4,0.1]])
   dna=np.delete(dna,1,0)
   scounter=1
   wall=150
   memory=0
   alpha=0
   beta=0
   indexme=0
   cid2=0.0
   pointx=0
   pointy=0
   pointz=0
   generation=0
   temperature=2
   round=0
   records=[]
   newsons=[]
   sons=[]
   parents=[]
   alive=True
   deathreason=0
   def record(self):
       lista=[self.m,self.wall,self.h,self.v,self.w,self.x,self.y,self.z,self.e,self.alpha,self.beta,self.memory,self.message,self.messageid,self.settingm,self.mutationpercentage,self.temperature]
       self.records.append(lista)
   def recordson(self,idtemp):
       self.sons.append(idtemp)
   def cr(self,a=0):
      iy=rd(0,5)
      #print(iy)
      if iy==1 and a==0:
         return [[117,rd(-1,16),rd(-1,5),rd(-10000000,10000000)/1000],self.cr(a=1)]
      elif iy==2 and a==0:
         return [[118,rd(-1,16),rd(-1,5),rd(-1,13)],self.cr(a=1)]
      elif a==1 :
        # ug=rd(1,4)
       #  print(ug)
         return [rd(1,41),rd(-10001,10000)/10000]
      else:
         return [rd(1,41),rd(-10001,10000)/10000]
   def aplyer(self,a):
      if self.e>abs(a):
         self.e-=abs(a)
         return a
      else :
       if self.stateofaply>=0:
         if a>0:
            a=self.e+0.1-0.1
            self.e=0
            return a
         else:
            a=self.e+0.1-0.1
            self.e=0
            return -1*a
       else:
           return 0;
   def mutat(self):
      #0 add 1 del 2 edit
      radiation=1
      if self.z>=0:
         radiation=3
      i=rd(-1,2+int(self.mutationpercentage/radiation))
      if i==0:
         #add
         indx=rd(-1,self.dna.shape[0]-1)
         
         gene= self.cr()
   
         return np.insert(self.dna,indx,gene,axis=0)
      elif i==1:
         #del
         indx=rd(-1,self.dna.shape[0]-1)
         return np.delete(self.dna,indx,0)
      else:
         #edit
         return self.dna
   def vv(self,a):
       #1
      # print("a")
       a=self.aplyer(a)
       #print(a)
       self.v+=a
      # print(self.v)
       
   def hh(self,a):
       #2
       a=self.aplyer(a)
       
       self.h+=a
   
   def mm(self,a):
    #3
    a=self.aplyer(a)
    self.m+=a/self.m
    if self.m==0:
       self.m=0.0001
   def setmessageid(self,n):
      #6
      a=self.aplyer(0.03)
      if a==0.03:
         if abs(n)<0.9:
            self.messageid=abs(n)
         else:
            self.setmessageid(n/100)
   def setmessagekey(self,n):
      #7
      a=self.aplyer(0.03)
      if a==0.03:
         if abs(n)<1:
            self.messagekey=abs(int(n*100))
         else:
            self.setmessagekey(n/100)
   def setmessage(self,n):
      #8
      a=self.aplyer(0.1)
      if a==0.1:
         listorders[int(abs(self.messageid)*1000000)+abs(self.messagekey)]=n
   def setqmessage(self,n):
      #9
      a=self.aplyer(0.1)
      if a==0.1:
         listorders[int(abs(self.messageid)*1000000)+abs(self.messagekey)]=n
   def setsettingm(self,n):
      #10
      a=self.aplyer(0.03)
      if a==0.03:
         self.settingm=n*100
   def mvv(self,n):
      #11
      self.vv(n*self.settingm*self.message)
   def mhh(self,n):
      #12
      self.hh(n*self.settingm*self.message)
   def mmm(self,n):
      #13
      self.mm(n*self.settingm*self.message)
  
   def copy(self,a):
    #4
    a=self.aplyer(a)
    if a>0:
       newc=cell()
       newc.dna=self.mutat()
       #newc.m=self.m+2
     #  print(xm)
       newc.xr=self.xr
       newc.yr=self.yr
       newc.zr=self.zr
       newc.h=self.h*0.8
       newc.v=self.v*0.8
       newc.w=self.w*0.8
       
       newc.x=self.x
       newc.y=self.y
       newc.z=self.z
       #rd(-1,1)*ydash
    #   print(ym)
       newc.e=a
       newc.v=self.v*0.8;newc.h=self.h*0.8;newc.w=self.w*0.8;
       newc.parents=self.parents.copy()
       newc.sons=[]
       newc.temperature=self.temperature
       newc.records=[]
       newc.parents.append(self.cid2)
       newc.cid2=rd(0,10**12)+0.0
       newc.generation=self.round
       self.recordson(newc.cid2)
       self.newsons.append(newc)
   def mcopy(self,n):
      #14
      self.copy(n*self.settingm*self.message)
   def setmessageidm(self,n):
      #16
      self.setmessageid(n*self.message)
   def setmessagekeym(self,n):
      #17
      self.setmessagekey(n*self.message)
   def setmessagem(self,n):
      #18
      self.setmessage(n*self.message*self.settingm)
   def setqmessagem(self,n):
       #19
      self.setqmessage(n*self.message*self.settingm)
   def setaplyer(self,n):
      #20
      if self.e>=0.03:
         self.e-=0.03
         self.stateofaply=n
   def setrefmessage(self,n):
      #21
      a=self.aplyer(0.1)
      self.e-=a
      n=int(abs(n)*13)
      if n==0:
         self.setmessage(self.x)
      elif n==1:
         self.setmessage(self.y)
      elif n==2:
         self.setmessage(self.e)
      elif n==3:
         self.setmessage(self.m)
      elif n==4:
         self.setmessage(self.v)
      elif n==5:
         self.setmessage(self.h)
      elif n==6:
         self.setmessage(self.memory)
      elif n==7:
         self.setmessage(self.z)
      elif n==8:
         self.setmessage(self.w)
      elif n==10:
         self.setmessage(self.alpha)
      elif n==11:
         self.setmessage(self.beta)
      elif n==12:
         self.setmessage(self.temperature)
   def setmutation(self,n):
      #22
      if self.e>=0.03:
         self.e-=0.03
         self.mutationpercentage=abs(int(n*1000))
   def setwall(self,n):
      #5
      self.wall+=(self.aplyer(n)/self.m)*80
   def msetwall(self,n):
      #24
      self.setwall(n*self.settingm*self.message)
 
   def setmemory(self,n):
      #25
      if self.e>=0.03:
         self.e-=0.03
         self.memory=n
   def setrefmemory(self,n):
    #33
      n=int(abs(n)*7)
      if n==0:
         self.setmemory(self.x)
      elif n==1:
         self.setmemory(self.y)
      elif n==2:
         self.setmemory(self.e)
      elif n==3:
         self.setmemory(self.m)
      elif n==4:
         self.setmemory(self.v)
      elif n==5:
         self.setmemory(self.h)
      elif n==6:
         self.setmemory(self.message)
   def mmsetsettingm(self,n):
      #26
      self.setsettingm(n*self.memory)
   def mmvv(self,n):
      #27
      self.vv(n*self.memory)
   def mmhh(self,n):
      #27
      self.hh(n*self.memory)
   def mmmm(self,n):
      #28
      self.mm(n*self.memory)
   def mmcopy(self,n):
      #29
      self.copy(n*self.memory)
   def mmsetmessageidm(self,n):
      #30
      self.setmessageid(n*self.memory)
   def mmsetmessagekeym(self,n):
      #31
      self.setmessagekey(n*self.memory)
   def mmsetwall(self,n):
      #32
      self.setwall(n*self.memory)
   def ww(self,a):
       #34
       a=self.aplyer(a)
       self.w+=a
   def mmww(self,n):
      #35
      self.ww(n*self.memory)
   def mww(self,n):
      #36
      self.ww(n*self.settingm*self.message)
   def alphaburn(self,n):
       #37
       if self.z>=0:
          if self.e>=0.01:
             self.e-=0.01
             self.e+=self.alpha*abs(n)*2*( 1/(1+2**(3-self.temperature) ))
             self.alpha-=self.alpha*abs(n)
             self.temperature+=0.1
       else :
           n=self.aplyer(0.3)
           self.e-=n;
           self.temperature+=0.7
   def betaburn(self,n):
       #38
      # print(self.cid2,self.z)
       if self.z<0:
       #   print(self.cid2,"True")
          if self.e>=0.07:
             self.e-=0.07
        #     print(self.cid2,"Trues")
             self.e+=self.beta*abs(n)*( 1/(1+2**(3-self.temperature) ))
             self.beta-=self.beta*abs(n)
             self.temperature+=0.1
       else :
           n=self.aplyer(0.3)
           self.e-=n;
           self.temperature+=0.7
   def settemperature(self,n):
       #39
       n=self.aplyer(n)
       self.e-=n
       self.temperature+=n
       if self.temperature<0:
          self.temperature=0
   def msettemperature(self,n):
       #40
       self.settemperature(n*self.settingm*self.message)
   def mmsettemperature(self,n):
       #41
       self.settemperature(n*self.memory)
   def ifst(self,a,b,c):
      if b==0:
         return a<c
      elif b==1:
         return a>c
      elif b==2:
         return a>=c
      elif b==3:
         return a<=c
      elif b==4:
         return a==c
      elif b==5:
         return a!=c
      else:
         return False
   def aplyvh(self,ydash):
   #   print(self.y,self.yr,self.v)
      self.xr+=self.h/2
      self.yr+=self.v/2
      self.zr+=self.w/2
      self.x=float(int(self.xr/ydash))*ydash
      self.y=float(int(self.yr/ydash))*ydash
      self.z=float(int(self.zr/ydash))*ydash
    #  print(ydash)
    #  print(self.y,self.yr,self.v)
      self.w-=0.3
      #Aply fraction . . .
      if self.z>=0:
         self.v*=0.95
         self.h*=0.95
         self.w*=0.95
      else:
         self.v*=0.85
         self.h*=0.85
         self.w*=0.85
   def test(self,gen):
     # print(type(1))
      if type(gen[0])==type(int(1)):
         print(gen)
         print(self.dna)
         return False
      else:
         return True
   
   def readdna(self,dna,dz=0):
   # try:
    scounter=0
    if dz==0:
       self.message=listorders[int(self.messageid*1000000)+self.messagekey]
       self.record()
  #  except:
     #  print(self.messageid,self.message,self.messagekey)
    a=0
    kd=1
    while a<len(dna) and self.e>0.02*1.4**a:
      gen=dna[a]
      #print(a,gen)
      self.e-=0.02*1.4**a
      a+=1
      if len(gen):
           if gen[0]==1:
              self.vv(gen[1])
              
           elif gen[0]==2:
              self.hh(gen[1])
              
           elif gen[0]==3:
              self.mm(gen[1])
              
           elif gen[0]==4:
              self.copy(gen[1])
           elif  gen[0]==23:
         #     print(gen[0])
              if dz==1 and self.e>=0.005*(2**kd):
                 self.e-=0.005*(2**kd)
                 kd+=1
                 a=int(len(dna)*abs( gen[1]))
           elif gen[0]==5:
                self.setwall(gen[1])
           elif gen[0]==6:
                self.setmessageid(gen[1])
           elif gen[0]==7:
                self.setmessagekey(gen[1])
           elif gen[0]==8:
                self.setmessage(gen[1])
           elif gen[0]==9:
                self.setqmessage(gen[1])
           elif gen[0]==10:
                self.setsettingm(gen[1])
           elif gen[0]==11:
                self.mvv(gen[1])
           elif gen[0]==12:
                self.mhh(gen[1])
           elif gen[0]==13:
                self.mmm(gen[1])
           elif gen[0]==14:
                self.mcopy(gen[1])
           elif  gen[0]==15:
                if dz==1 and self.e>=0.005*(2**kd) and self.message!=0:
                 self.e-=0.005*(2**kd)
                 kd+=1
                 a=int(len(dna)*abs( gen[1]*self.message)) 
           elif gen[0]==16:
                self.setmessageidm(gen[1])
           elif gen[0]==17:
                self.setmessagekeym(gen[1])
           elif gen[0]==18:
                self.setmessagem(gen[1])
           elif gen[0]==19:
                self.setqmessagem(gen[1])
           elif gen[0]==20:
                self.setaplyer(gen[1])
           elif gen[0]==21:
                self.setrefmessage(gen[1])
           elif gen[0]==22:
                self.setmutation(gen[1])
           elif gen[0]==24:
                self.msetwall(gen[1])
           elif gen[0]==25:
                self.setmemory(gen[1])
           elif gen[0]==26:
                self.mmsetsettingm(gen[1])
           elif gen[0]==27:
                self.mmvv(gen[1])
           elif gen[0]==28:
                self.mmhh(gen[1])
           elif gen[0]==29:
                self.mmmm(gen[1])
           elif gen[0]==30:
                self.mmsetmessageidm(gen[1])
           elif gen[0]==31:
                self.mmsetmessagekeym(gen[1])
           elif gen[0]==32:
                self.mmsetwall(gen[1])
           elif gen[0]==33:
                self.setrefmemory(gen[1])
           elif gen[0]==34:
                self.ww(gen[1])
           elif gen[0]==35:
                self.mmww(gen[1])
           elif gen[0]==36:
                self.mww(gen[1])
           elif gen[0]==37:
                self.alphaburn(gen[1])
           elif gen[0]==38:
                self.betaburn(gen[1])
           elif gen[0]==39:
                self.settemperature(gen[1])
           elif gen[0]==40:
                self.msettemperature(gen[1])
           elif gen[0]==41:
                self.mmsettemperature(gen[1])
           else:
              #  print(gen,self.cid2)
                dict=[self.v,self.h,self.m,self.e,self.x,self.y,self.messageid,self.messagekey,self.message,self.settingm,self.stateofaply,self.wall,self.w,self.z,self.alpha,self.beta,self.temperature]
                if gen[0][0]==117:
                   
                   if self.ifst(dict[gen[0][1]],gen[0][2],gen[0][3]):
                        self.readdna(np.array(gen[1:]),1)  
        
                elif gen[0][0]==118:  
                   
                   if self.ifst(dict[gen[0][1]],gen[0][2],dict[gen[0][3]]):
                        self.readdna(np.array(gen[1:]),1)
def createlist(sx,sy,sz,f=1):
        mylist=[None]*(sx*sy*sz)
        return mylist
def modify(sx,sy,sz,lista,add):
       # print(x,y,z)
        #print(lista)
#        tempcopy=lista[x][y][z].copy()
 #       print("t",tempcopy)
       # tempcopy.append(add)
       # print("tt",tempcopy,lista[x][y][z])
        #lista[x][y][z]=tempcopy
        if lista[sx*sy*sz]!=None:
           lista[sx*sy*sz].append(add)
        else:
            lista[sx*sy*sz]=[add]
        #print("Lista",lista[sx*sy*sz])
        return lista
def addonebyone(list1,list2):
    for l in list2:
        list1.append(l)
    return list1
def aplythread1(listtempa,ep,alpha,beta,xmax,ymax,zmax,cl,round,airt,watert,aplydna=False):
    listnew=[]
    listofxsyszs=[]
    #print("listtempa ",listtempa)
    for ilist in listtempa:
        ilist.newsons=[]
        ilist.aplyvh(cl)
        if ilist.z>=0:
           ilist.temperature+=(airt-ilist.temperature)*0.1
        else:
           ilist.temperature+=(watert-ilist.temperature)*0.1
        if ilist.x>xmax or ilist.x<-xmax or ilist.y>ymax or ilist.y<-ymax or ilist.z>zmax or ilist.z<-zmax:
           ilist.alive=False;
           ilist.deathreason=0
        else:
           if aplydna:
              ilist.round=round
              ilist.readdna(ilist.dna)
           sx=(xmax+ilist.x)/cl-(xmax+ilist.x)%cl/cl
           sy=(ymax+ilist.y)/cl-(ymax+ilist.y)%cl/cl
           sz=(zmax+ilist.z)/cl-(zmax+ilist.z)%cl/cl
           listofxsyszs.append([int(sx),int(sy),int(sz),ilist])
           if ilist.newsons!=[]:
              listnew=addonebyone(listnew,ilist.newsons)
    for ilist  in listnew:
        sx=(xmax+ilist.x)/cl-(xmax+ilist.x)%cl/cl
        sy=(ymax+ilist.y)/cl-(ymax+ilist.y)%cl/cl
        sz=(zmax+ilist.z)/cl-(zmax+ilist.z)%cl/cl
        ilist.e+=ep
        ilist.alpha+=alpha
        ilist.beta+=beta
        
        listofxsyszs.append([int(sx),int(sy),int(sz),ilist])
    return listofxsyszs , listnew
def aplythread2(listarranged):
  listalive=[]
  listx=[]
  listy=[]
  listz=[]
  #print(listarranged)
  i=0
  for liu in listarranged:
      listalive1=[]
      if liu!=None:
        for liu2 in liu:
                h=liu2.wall
                liu2.wall-=1*(liu2.temperature+1)*len(liu)
                
                if liu2.wall>0:
                   if listalive1==[]:
                       #print("rrrrrrrr")
                       listalive1.append(liu2)
                       listx.append(liu2.x)
                       listy.append(liu2.y)
                       listz.append(liu2.z)
                   else:
                       if liu2.m>listalive1[0].m:
                          for li in listalive1:
                              liu2.alpha+=li.alpha*0.5
                              liu2.e+=li.e*0.3
                              liu2.beta+=li.beta*0.5
                              li.alive=False;
                              li.deathreason=2 #ate
                          listalive1=[liu2]
                       elif liu2.m<listalive1[0].m:
                            liu2.alive=False
                            liu2.deathreason=2 #ate
                            listalive1[len(listalive1)-1].alpha=liu2.alpha*0.5
                            listalive1[len(listalive1)-1].beta=liu2.beta*0.5
                            listalive1[len(listalive1)-1].e=liu2.e*0.3
                       else :
                            listalive1.append(liu2)
                else:
                    liu2.alive=False
                    liu2.deathreason=1 #wall break
                    if len(listalive1)>0:  
                       listalive1[len(listalive1)-1].alpha=liu2.alpha*0.5
                       listalive1[len(listalive1)-1].beta=liu2.beta*0.5
                       listalive1[len(listalive1)-1].e=liu2.e*0.3
        listalive=addonebyone(listalive,listalive1)
  return listalive,listx,listy,listz
print("Valuable cells V2.0 ")
def dome(n_threads):
    iu=int(input("press zero to loop  1 to load saved data: "))
    if iu==1:
       strfile=input("Enter file name : ")
       f=open(strfile,'rb')
       listcells,listall,watert,airt,factor1,factor2,a=pickle.load(f)
    else:
        listall=[]
        listcells=[]
        airt=0
        watert=0
        factor1=1
        factor2=1
        a=0
    ep=float(input("enter energy per time : "))
    alpha=float(input("enter alpha per time : "))
    beta=float(input("enter beta per time : "))
    maxairt=float(input("enter Maximum air temperature : "))
    maxwatert=float(input("enter maximum water temperature : "))
    ymax=int(input("enter Y axis : "))
    xmax=int(input("enter X axis : "))
    zmax=int(input("enter Z axis : "))
    cl=float(input("enter crash limit : "))
    times = int(input("enter number of times : "))
    maxi=0
    maxi=int(input("enter Maximum number of cells : "))
    c=cell()
    listcells.append(c)
    listall.append(c)
    dxmax=2*xmax
    sx=2*xmax/cl-(2*xmax)%cl/cl
    sy=2*ymax/cl-(2*ymax)%cl/cl
    sz=2*zmax/cl-(2*zmax)%cl/cl
    sx=int(sx)
    sy=int(sy)
    sz=int(sz)
    listarranged2=createlist(sx,sy,sz,f=1)
    listarranged=listarranged2.copy()
    
 
    #print(len(listarranged),len(listarranged[0]),len(listarranged[0][0]))
    while True:
      time.sleep(0.2)
      a+=1
      #print("1",listcells)
      airt+=0.2*factor1
      watert+=0.05*factor2
      if watert>=maxwatert or watert<0:
         factor2*=-1
      if airt>=maxairt or airt<0:
          factor1*=-1
      listofxsyszs,listnew=aplythread1(listcells,ep,alpha,beta,xmax,ymax,zmax,cl,a/3-a%3,airt,watert,aplydna=a%3==0)
      #
      #print("xsl",listofxsyszs,listnew)
      listall=addonebyone(listall,listnew)
      #print(listall)
      del(listnew)
      i=0
      listarranged=listarranged2.copy()
      #print(listarranged)
      while i< len(listofxsyszs):
       #   print("begin")
          li=listofxsyszs[i]
          #print(li[0],li[1],li[2])
        #  print(listarranged[0][0][0])
         # print(listarranged[li[0]][li[1]][li[2]])
        #  print(listarranged)
          listarranged=modify(li[0],li[1],li[2],listarranged,li[3])
          #print("dddddddddddd  ====",listarranged2[li[0]*li[1]*li[2]])
         # print(listarranged)
          #print(listarranged[li[0]][li[1]][li[2]])
          #print(listarranged[0][0][0])
          i+=1
      #print("Done")
      #del(listarranged)
      del(listofxsyszs)
      listcells,listx,listy,listz=aplythread2(listarranged)
      del(listarranged)
      #print("cells",listcells)
      print("cells : ",len(listcells))
      print("Time : ",a)
      print("WaterT : ",watert)
      print("AirT : ",airt)
      #fig2 = plt.figure(figsize=(7,7))
      #ax2 = fig2.add_subplot(111, projection='3d')
      #ax2.scatter(listx,listy,listz)
      xf="";yf="";zf="";
      try:
       file1=open("x.txt","w+")
       file2=open("y.txt","w+")
       file3=open("z.txt","w+")
       for i in range(len(listx)):
         xf+=(str(listx[i]))
         yf+=(str(listy[i]))
         zf+=(str(listz[i]))
         if i<len(listx)-1:
            xf+="\n"
            yf+="\n"
            zf+="\n"
       file1.writelines(xf)
       file2.writelines(yf)
       file3.writelines(zf)
       file1.close();
       file2.close();
       file3.close();
      except :
         print("\nerror\n")
      del(xf)
      del(yf)
      del(zf)
      del(listx)
      del(listy)
      del(listz)
      #ax2.axis([-xmax,xmax,-ymax,ymax])
      #plt.show()
      if a==times or len(listcells)>=maxi or len(listcells)==0 :
       print("variation ended")
      #print(len(listcells)-len(ist2),len(listcells),len(ist2),id(ist2))
       while True:
         iu=int(input("press zero to loop again 1 to save the data: "))
         if iu==0:
            ep=float(input("enter energy per time : "))
            alpha=float(input("enter alpha per time : "))
            beta=float(input("enter beta per time : "))
            maxairt=float(input("enter Maximum air temperature : "))
            maxwatert=float(input("enter maximum water temperature : "))
            ymax=int(input("enter Y axis : "))
            xmax=int(input("enter X axis : "))
            zmax=int(input("enter Z axis : "))
            cl=float(input("enter crash limit : "))
            times += int(input("enter number of more times : "))
            maxi= int(input("enter maximum number of cells : "))
            dxmax=2*xmax
            sx=(2*xmax+cl)/cl-(2*xmax+cl)%cl/cl
            sy=(2*ymax+cl)/cl-(2*ymax+cl)%cl/cl
            sz=(2*zmax+cl)/cl-(2*zmax+cl)%cl/cl
            sx=int(sx)
            sy=int(sy)
            sz=int(sz)
            listarranged2=createlist(sx,sy,sz,f=1)
            listarranged=listarranged2.copy()
            break;
         elif iu==1:
          iku=0
          r=0
          strfile=input("Enter file name : ")
          f=open(strfile,'wb')
          pickle.dump([listcells,listall,watert,airt,factor1,factor2,a],f)
          f.close()
dome(1)