import finesse
_loss = 5e-6
_etmT = 15e-6
def gmi(power=125,wlen=1064,Larm=4e3,phi1=0.03,phi2=89.97,phi3=0,mass=40,phi0=0,lossy=1,nsr=1,pr=0,prmT=0.03,dphi=180,sr=0,pendula=1):
    nsr = 'True' if nsr else 'False'
    loss = 0
    etmT = 0
    if lossy:
        loss = _loss
        etmT = _etmT
    kat = finesse.Model()
    kat.parse(f'''
    lambda({wlen}n)
    var Larm {Larm}
    var Mtm  {mass}
    var phi2 {phi2}
    var phi1 {phi1}
    var etmT {etmT}
    var l0 {loss}
    var phi0 {phi0}
    var phi3 {phi3}
    var dphi {dphi}
    

    laser l1 {power}
    
    s s0 l1.p1 prm.p1 L=1
    s prc prm.p2 BS1.p1 L=1
    m prm R=0 L=0 phi=0
    
    bs BS1 R=0.5 T=0.5 L=0 phi=0 alpha=45
    bs BS2 R=0.5 T=0.5 L=0 phi=0 alpha=45
    s s1 BS1.p2 m1.p1 L=1
    s s2 BS2.p4 m2.p1 L=1
    s bridge BS1.p3 BS2.p1 L=1
    s syarm BS2.p2 ETMy.p1 L=Larm
    s sxarm BS2.p3 ETMx.p1 L=Larm
    m ETMy T=etmT L=l0 phi=phi2
    m ETMx T=etmT L=l0 phi=phi1
    m m1 T=etmT L=l0 phi=phi0 # GMI homodyne mirror
    m m2 T=etmT L=l0 phi=phi3 # GMI feedback mirror

    # Apply a signal to each arm with 180 phase difference between them
    fsig(1)
    sgen darmx syarm.h phase=dphi
    sgen darmy sxarm.h
    
    # output the upper sideband
    ad signal BS1.p4.o f=fsig
    pd output BS1.p4.o
    qshot NSR_without_RP BS1.p4.o nsr={nsr}
    qnoised NSR_with_RP BS1.p4.o nsr={nsr}
    ''')
    if pr:
        kat.set('prm.T',prmT)
        kat.set('prm.phi', 90)
        if lossy:
            kat.set('prm.L',loss)
    if sr:
        kat.parse('''
                m srm T=0.2 L=l0 phi=-90
                s src BS1.p4 srm.p1 L=1
                ''')
        [kat.remove(i) for i in ['signal','output','NSR_without_RP','NSR_with_RP']]
        kat.parse(f'''
        ad signal srm.p2.o f=fsig
        pd output srm.p2.o
        qshot NSR_without_RP srm.p2.o nsr={nsr}
        qnoised NSR_with_RP srm.p2.o nsr={nsr}
            ''')
    if pendula:
        kat.parse('''
        pendulum etmx_sus ETMx.mech mass=Mtm fz=1 Qz=1M
        pendulum etmy_sus ETMy.mech mass=Mtm fz=1 Qz=1M
            ''')
    return kat

def ligo(power=125,wlen=1064,Larm=4e3,lossy=1,dphi=180,pendula=1,nsr=1):
    nsr = 'True' if nsr else 'False'
    etmT = 0
    loss = 0
    if lossy:
        loss = _loss
        etmT = _etmT
    kat = finesse.Model()
    kat.parse(f'''
        lambda({wlen}n)
        ###########################################################################
        ###   Variables
        ###########################################################################
        var Larm {Larm}
        var Mtm  40
        var itmT 0.014
        var lmichx 4.5
        var lmichy 4.45
        var l0 {loss}
        var etmT {etmT}
        var phi1 0
        var dphi {dphi}
        

        ###########################################################################
        ###   Input optics
        ###########################################################################
        l L0 {power}
        s l_in L0.p1 prm.p1
        # Power recycling mirror
        m prm T=0.03 L=l0 phi=90
        s prc prm.p2 bs.p1 L=53


        # Central beamsplitter
        bs bs T=0.5 L=0 alpha=45

        ###########################################################################
        ###   X arm
        ###########################################################################
        s lx bs.p3 itmx.p1 L=lmichx
        m itmx T=itmT L=l0 phi=90+phi1
        s LX itmx.p2 etmx.p1 L=Larm
        m etmx T=etmT L=l0 phi=89.999875

        ###########################################################################
        ###   Y arm
        ###########################################################################
        s ly bs.p2 itmy.p1 L=lmichy
        m itmy T=itmT L=l0 phi=0
        s LY itmy.p2 etmy.p1 L=Larm
        m etmy T=etmT L=l0 phi=0.000125  
        

        ###########################################################################
        ###   Output and squeezing
        ###########################################################################
        s src bs.p4 srm.p1 L=50.525
        m srm T=0.2 L=l0 phi=-90

        # A squeezed source could be injected into the dark port
        sq sq1 db=0 angle=90
        s lsqz sq1.p1 srm.p2

        # Differentially modulate the arm lengths
        fsig(1)
        sgen darmx LX.h
        sgen darmy LY.h phase=dphi

        # Output the full quantum noise limited sensitivity
        qnoised NSR_with_RP srm.p2.o nsr={nsr}
        # Output just the shot noise limited sensitivity
        qshot NSR_without_RP srm.p2.o nsr={nsr}

        # We could also display the quantum noise and the signal
        # separately by uncommenting these two lines.
        # qnoised noise srm.p2.o
        ad signal srm.p2.o f=fsig
        pd output srm.p2.o
        '''
    )
    if pendula:
        kat.parse('''
        pendulum itmx_sus itmx.mech mass=Mtm fz=1 Qz=1M
        pendulum etmx_sus etmx.mech mass=Mtm fz=1 Qz=1M
        pendulum itmy_sus itmy.mech mass=Mtm fz=1 Qz=1M
        pendulum etmy_sus etmy.mech mass=Mtm fz=1 Qz=1M
        ''')
    return kat

def mi(power=125,Larm=4e3,lossy=1,pendula=0,nsr=1):
    etmT = 0
    loss = 0
    if lossy:
        loss = _loss
        etmT = _etmT
    kat = finesse.Model()
    kat.parse(f'''
    var etmT {etmT}
    var Larm {Larm}
    var l0 {loss}
    
    laser l1 {power}
    
    s s0 l1.p1 BS.p1 L=1
    bs BS R=0.5 T=0.5 L=0 alpha=45
    s syarm BS.p2 ETMy.p1 L=Larm
    m ETMy T=etmT L=l0 phi=0
    s sxarm BS.p3 ETMx.p1 L=Larm
    m ETMx T=etmT L=l0 phi=90
    # Apply a signal each arm with 180 phase difference between them
    fsig(1)
    sgen darmx syarm.h phase=180
    sgen darmy sxarm.h 
    # output the upper sideband
    ad signal BS.p4.o f=fsig
    pd output BS.p4.o
    qnoised NSR_with_RP BS.p4.o nsr=False
    qshot NSR_without_RP BS.p4.o nsr=False
    ''')
    if pendula:
        kat.parse('''
        pendulum etmx_sus ETMx.mech mass=Mtm fz=1 Qz=1M
        pendulum etmy_sus ETMy.mech mass=Mtm fz=1 Qz=1M
        ''')
    if nsr:
        kat.set('NSR_with_RP.nsr', True)
        kat.set('NSR_without_RP.nsr', True)
    return kat
