from djinn import *
import random,math,itertools
pi=math.pi
def midpoint(x1, y1,z1, x2, y2,z2):
    return ((x1 + x2)/2, (y1 + y2)/2, (z1 + z2)/2)
def calc_dist(p1,p2):
    return math.sqrt((p2[0] - p1[0]) ** 2 +
                     (p2[1] - p1[1]) ** 2 +
                    (p2[2] - p1[2]) ** 2)

     
class proton():
    def __init__(self):
        self.fnet=0
        self.Q=1.60217656535*10**-19
        self.r=  8.768*10**-16
        self.A=4*pi*self.r**2
        self.m=1.6726219*10**-27
        self.v=(4/3)*pi*self.r**3
        self.p=2*pi*self.r
        self.S=Sphere(1,random.randrange(-5,5),random.randrange(-5,5),random.randrange(-5,5),(1,0,0))
        self.pos=[self.S._x,self.S._y,self.S._z]
class neutron():
    def __init__(self):
        self.Q=0.00000000001*10**-30
        self.r=8.768*10**-16
        self.A=4*pi*self.r**2
        self.m=1.67492747121*10**-27
        self.v=(4/3)*pi*self.r**3
        self.p=2*pi*self.r
        self.S=Sphere(1,random.randrange(-5,5),random.randrange(-5,5),random.randrange(-5,5),(0,1,0))
        self.pos=[self.S._x,self.S._y,self.S._z]
class electron():
    def __init__(self):
        self.Q=-1.60217656535*10**-19
        self.r=2.817940326727 * 10**-15
        self.A=4*pi*self.r**2
        self.m= 9.10938291 * 10**-31
        self.v=(4/3)*pi*self.r**3
        self.p=2*pi*self.r
        self.M=9.10938291 * 10**-31
        self.S=Sphere(1,random.randrange(-5,5),random.randrange(-5,5),random.randrange(-5,5),(0,0,1))
        self.pos=[self.S._x,self.S._y,self.S._z]
class atom():
    def __init__(self):
        self.N=0
        self.ol=[]
        self.Q=self.N*-1.60217656535*10**-19
    
    def M(self):
        tm=0
        for i in self.ol:
            tm+=i.m
        return tm
    def com(self):
        X=0
        Y=0
        Z=0
        mt=0
        mgx=0
        mgy=0
        mgz=0
        ll=[]
        ml=[]
        ML=0
        t1=0
        popr=0
        for i in self.ol:
                ll.append(i.pos)
        for i in ll:
            try:
                ll.pop(popr)
                ll.pop(popr)
                ll.pop(popr)
                popr+=1
            except IndexError:
                pass
        for i in self.ol:
            ml.append(i.m)
            ML+=i.m
            
        for i in ml:
            try:
                mgx+=i*(ll[t1][0])
                mgy+=i*(ll[t1][1])
                mgz+=i*(ll[t1][2])
                t1+=1
            except IndexError:
                pass
        X+=mgx/ML
        Y+=mgy/ML
        Z+=mgz/ML
        return [X,Y,Z]    
    def bs(self,n):
        v=n
        while v!=0:
            exec('self.p'+str(self.N)+'=proton()')
            exec('self.n'+str(self.N)+'=neutron()')
            exec('self.e'+str(self.N)+'=electron()')
            self.ol.append(eval('self.p'+str(self.N)))
            self.ol.append(eval('self.n'+str(self.N)))
            self.ol.append(eval('self.e'+str(self.N)))
            self.N+=1
            v-=1
    def ds(self):
        for i in self.ol:
            i.S.build()
class atmcnt():
    def __init__(self):
        self.c=0
        self.la=[]
        self.lM=[]
atm=atmcnt()
def spwn(a):
    globals()["atm"+str(atm.c)]=atom()
    atm.la.append(eval("atm"+str(atm.c)))
    exec('atm'+str(atm.c)+'.bs('+str(a)+')')
    atm.c+=1




def gea(m):
    k=8.9875517873681764*10**9
    G=6.671281903963040991511534289*10**-11
    cl=itertools.combinations(m.ol,2)    
    for i in cl:
        try:
            mp=midpoint(i[0].S._x,i[0].S._y,i[0].S._z,i[1].S._x,i[1].S._y,i[1].S._z)
            R=(calc_dist(i[0].pos,i[1].pos)/2)*10**-13           
            F=k*(i[0].Q*i[1].Q/R**2)
            gF=G*(i[0].m*i[1].m/R**2)    
            fstr=str(F)
            gstr=str(gF)
            pf=fstr.partition('e')
            pgf=gstr.partition('e')
            pvf=(int(pf[2]))-(2*int(pf[2]))#+9.5
            pvgf=(int(pgf[2]))-(2*int(pgf[2]))#+9.5           
            #print R
            if i[0].S._x>mp[0]:
                i[0].S._x=i[0].S._x-(gF*10**pvgf)
            if i[1].S._x>mp[0]:
                i[1].S._x=i[1].S._x-(gF*10**pvgf)
            if i[0].S._x<mp[0]:
                i[0].S._x=i[0].S._x+(gF*10**pvgf)
            if i[1].S._x<mp[0]:
                i[1].S._x=i[1].S._x+(gF*10**pvgf)

            if i[0].S._y>mp[1]:
                i[0].S._y=i[0].S._y-(gF*10**pvgf)
            if i[1].S._y>mp[1]:
                i[1].S._y=i[1].S._y-(gF*10**pvgf)
            if i[0].S._y<mp[1]:
                i[0].S._y=i[0].S._y+(gF*10**pvgf)
            if i[1].S._y<mp[1]:
                i[1].S._y=i[1].S._y+(gF*10**pvgf)

            if i[0].S._z>mp[2]:
                i[0].S._z=i[0].S._z-(gF*10**pvgf)
            if i[1].S._z>mp[2]:
                i[1].S._z=i[1].S._z-(gF*10**pvgf)
            if i[0].S._z<mp[2]:
                i[0].S._z=i[0].S._z+(gF*10**pvgf)
            if i[1].S._z<mp[2]:
                i[1].S._z=i[1].S._z+(gF*10**pvgf)
                
            if i[0].S._x<mp[0]:
                i[0].S._x=i[0].S._x-(F*10**pvf)
            if i[1].S._x<mp[0]:
                i[1].S._x=i[1].S._x-(F*10**pvf)
            if i[0].S._x>mp[0]:
                i[0].S._x=i[0].S._x+(F*10**pvf)
            if i[1].S._x>mp[0]:
                i[1].S._x=i[1].S._x+(F*10**pvf)

            if i[0].S._y<mp[1]:
                i[0].S._y=i[0].S._y-(F*10**pvf)
            if i[1].S._y<mp[1]:
                i[1].S._y=i[1].S._y-(F*10**pvf)
            if i[0].S._y>mp[1]:
                i[0].S._y=i[0].S._y+(F*10**pvf)
            if i[1].S._y>mp[1]:
                i[1].S._y=i[1].S._y+(F*10**pvf)

            if i[0].S._z<mp[2]:
                i[0].S._z=i[0].S._z-(F*10**pvf)
            if i[1].S._z<mp[2]:
                i[1].S._z=i[1].S._z-(F*10**pvf)
            if i[0].S._z>mp[2]:
                i[0].S._z=i[0].S._z+(F*10**pvf)
            if i[1].S._z>mp[2]:
                i[1].S._z=i[1].S._z+(F*10**pvf)

        except ValueError:
            pass
def gem(m):
    k=8.9875517873681764*10**9
    G=6.671281903963040991511534289*10**-11
    cl=itertools.combinations(m.al,2)    
    for i in cl:
        try:
            mp=midpoint(i[0].com()[0],i[0].com()[1],i[0].com()[2],i[1].com()[0],i[1].com()[1],i[1].com()[2])
            R=(calc_dist(i[0].com(),i[1].com())/2)*10**-13           
            F=k*(i[0].Q*i[1].Q/R**2)
            gF=G*(i[0].M()*i[1].M()/R**2)    
            fstr=str(F)
            gstr=str(gF)
            pf=fstr.partition('e')
            pgf=gstr.partition('e')
            pvf=(int(pf[2]))-(2*int(pf[2]))#+9.5
            pvgf=(int(pgf[2]))-(2*int(pgf[2]))#+9.5           
            print R
            if i[0].S._x>mp[0]:
                i[0].S._x=i[0].S._x-(gF*10**pvgf)
            if i[1].S._x>mp[0]:
                i[1].S._x=i[1].S._x-(gF*10**pvgf)
            if i[0].S._x<mp[0]:
                i[0].S._x=i[0].S._x+(gF*10**pvgf)
            if i[1].S._x<mp[0]:
                i[1].S._x=i[1].S._x+(gF*10**pvgf)

            if i[0].S._y>mp[1]:
                i[0].S._y=i[0].S._y-(gF*10**pvgf)
            if i[1].S._y>mp[1]:
                i[1].S._y=i[1].S._y-(gF*10**pvgf)
            if i[0].S._y<mp[1]:
                i[0].S._y=i[0].S._y+(gF*10**pvgf)
            if i[1].S._y<mp[1]:
                i[1].S._y=i[1].S._y+(gF*10**pvgf)

            if i[0].S._z>mp[2]:
                i[0].S._z=i[0].S._z-(gF*10**pvgf)
            if i[1].S._z>mp[2]:
                i[1].S._z=i[1].S._z-(gF*10**pvgf)
            if i[0].S._z<mp[2]:
                i[0].S._z=i[0].S._z+(gF*10**pvgf)
            if i[1].S._z<mp[2]:
                i[1].S._z=i[1].S._z+(gF*10**pvgf)
                
            if i[0].S._x<mp[0]:
                i[0].S._x=i[0].S._x-(F*10**pvf)
            if i[1].S._x<mp[0]:
                i[1].S._x=i[1].S._x-(F*10**pvf)
            if i[0].S._x>mp[0]:
                i[0].S._x=i[0].S._x+(F*10**pvf)
            if i[1].S._x>mp[0]:
                i[1].S._x=i[1].S._x+(F*10**pvf)

            if i[0].S._y<mp[1]:
                i[0].S._y=i[0].S._y-(F*10**pvf)
            if i[1].S._y<mp[1]:
                i[1].S._y=i[1].S._y-(F*10**pvf)
            if i[0].S._y>mp[1]:
                i[0].S._y=i[0].S._y+(F*10**pvf)
            if i[1].S._y>mp[1]:
                i[1].S._y=i[1].S._y+(F*10**pvf)

            if i[0].S._z<mp[2]:
                i[0].S._z=i[0].S._z-(F*10**pvf)
            if i[1].S._z<mp[2]:
                i[1].S._z=i[1].S._z-(F*10**pvf)
            if i[0].S._z>mp[2]:
                i[0].S._z=i[0].S._z+(F*10**pvf)
            if i[1].S._z>mp[2]:
                i[1].S._z=i[1].S._z+(F*10**pvf)

        except ValueError:
            pass
class molecule():
    def __init__(self):
        self.al=[]
        #self.mf=self.com()
    def com(self):
        X=0
        Y=0
        Z=0
        mt=0
        mgx=0
        mgy=0
        mgz=0
        ll=[]
        ml=[]
        t1=0
        popr=0
        for i in self.al:
                ll.append(i.com())
        for i in ll:
            try:
                ll.pop(popr)
                ll.pop(popr)
                ll.pop(popr)
                popr+=1
            except IndexError:
                pass
        for i in self.al:
            ml.append(i.M())
            
        for i in ml:
            try:
                mgx+=i*(ll[t1][0])
                mgy+=i*(ll[t1][1])
                mgz+=i*(ll[t1][2])
                t1+=1
            except IndexError:
                pass
        
        X+=mgx/sum(ml)
        Y+=mgy/sum(ml)
        Z+=mgz/sum(ml)
        return [X,Y,Z]
    
    def sets(self,N):
        
        n=N
        y=0
        while y==0:
            if n>0:
                spwn(n)
                self.al.append(eval("atm"+str(atm.c-1)))
                n-=1
            elif n==0:
                y+=1
        
#spwn(1)
#atm0.ol[0].m=5.9723*10**24
#spwn(1)
mol1=molecule()
mol1.sets(3)
atm.lM.append(mol1)
mol2=molecule()
mol2.sets(1)
atm.lM.append(mol2)
mol3=molecule()
mol3.sets(3)
atm.lM.append(mol3)
mol4=molecule()
mol4.sets(1)
atm.lM.append(mol4)
def geA():
    for i in atm.la:
        i.ds()
        gea(i)
def geM():
    for i in atm.lM:
        gem(i)

#he=atom()
#he.bs(2)
if __name__=="__main__":
    window = Window((800,600))
    window.start(70)
    window.caption("Atomix")  
    #window.icon('djinn-iv-logo.bmp')
    play = Player(0,0,-10.0)
    
    
    
    
    play.setOrigin(0,0,0)
    light0 = Light(-20,5,-5,[1,1,1,1],1)
    light1 = Light(0,1,15,[1,1,1,1],1)
    light0.bake(GL_LIGHT0)
    light1.bake(GL_LIGHT1)
    moveList = [0,0,0]
    keymap = {'w':[0,0,25],'s':[0,0,-25],'a':[25,0,0],'d':[-25,0,0],'up':[0,-25,0],'down':[0,25,0]}
    x=0
    while True:
        KeyboardEvent(moveList,keymap)
        window.clear()
        play.move(moveList[0],moveList[1],moveList[2])
        geA()
        geM()
        window.update()
