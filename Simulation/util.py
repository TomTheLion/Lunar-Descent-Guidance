import numpy as np
import matplotlib.pyplot as plt
from math import factorial


#rodrigues vector rotation formula
def rodrigues(vector, axis, angle):
    axis = axis / np.linalg.norm(axis)
    return vector * np.cos(angle) + np.cross(axis, vector) * np.sin(angle) + axis * np.dot(axis, vector) * (1 - np.cos(angle))


#fourth order runge kutta
def RK4(r0, v0, a, h, gm):
    def f(r):
        g = -gm / np.linalg.norm(r) ** 3 * r
        return a + g
    kv1 = f(r0)
    kr1 = v0
    kv2 = f(r0 + kr1 * h / 2)
    kr2 = v0 + kv1 * h / 2
    kv3 = f(r0 + kr2 * h / 2)
    kr3 = v0 + kv2 * h / 2
    kv4 = f(r0 + kr3 * h)
    kr4 = v0 + kv3 * h   
    v1 = v0 + h / 6 * (kv1 + 2 * kv2 + 2 * kv3 + kv4)
    r1 = r0 + h / 6 * (kr1 + 2 * kr2 + 2 * kr3 + kr4)   
    return r1, v1


#state transition matrix
def STM(t1, t0):
    dt = t1 - t0
    xlist = []
    ylist = []
    for i in range(5):
        for j in range(5):
            if j >= i:
                yval = j - i
                xval = dt ** yval / factorial(yval)
            else:
                xval = 0
            xlist.append(xval)
            xval = 0
        ylist.append(xlist)
        xlist = []
    return np.vstack(ylist)


def aoplan(inputs):
    i = inputs['inclination']
    psi = inputs['psi']

    aop = np.arcsin(abs(np.sin(np.pi / 2 - psi)) / np.sin(i))
    if np.pi / 2 - psi < 0:
        aop += np.pi
    lan = np.arcsin(abs(np.tan(np.pi / 2 - psi)) / np.tan(i))
    print('AOP: ' + str(aop * 180 / np.pi) + ', LAN: ' + str(lan * 180 / np.pi))
    return 0

    
#calculates state vector in geocentric coordinates
def statevector(inputs):
    #calculate orbital parameters
    gm = inputs['gravparam']
    ra = inputs['bodyradius'] + inputs['apoapsis']
    rp = inputs['bodyradius'] + inputs['periapsis']
    i = inputs['inclination']
    ta = inputs['trueanomaly']
    theta = inputs['theta']
    psi = inputs['psi']
    a = (ra + rp) / 2
    e = (ra - rp) / (ra + rp)
    aop = np.arcsin(abs(np.sin(np.pi / 2 - psi)) / np.sin(i))
    lan = theta - np.arcsin(abs(np.tan(np.pi / 2 - psi)) / abs(np.tan(i)))
    #only works for inclination <90
    if i < np.pi / 2:
            if np.pi / 2 - psi < 0:
                aop = 2 * np.pi - aop
                lan = 2 * np.pi - lan
                
    else:
        if np.pi / 2 - psi < 0:
            aop = np.pi + aop
            lan = np.pi - lan
        else:
            aop = np.pi - aop
            lan = np.pi + lan

    h = rp * np.sqrt(gm * (2 / rp - 1 / a))
    r = h ** 2 / gm / (1 + e * np.cos(ta)) * np.array([np.cos(ta), np.sin(ta), 0])
    v = gm / h * np.array([-np.sin(ta), e + np.cos(ta), 0])
    
    #calculate transformation matrix from perifocal frame to geocentric frame
    qXx1 = np.array([[np.cos(aop), np.sin(aop), 0], [-np.sin(aop), np.cos(aop), 0], [0, 0, 1]])
    qXx2 = np.array([[1, 0, 0], [0, np.cos(i), np.sin(i)], [0, -np.sin(i), np.cos(i)]])
    qXx3 = np.array([[np.cos(lan), np.sin(lan), 0], [-np.sin(lan), np.cos(lan), 0], [0, 0, 1]])
    qXx = np.matmul(np.matmul(qXx1, qXx2), qXx3)
    qxX = np.transpose(qXx)
    
    #calculate state vector in geocentric coordinates
    rG = np.matmul(qxX, r)
    vG = np.matmul(qxX, v)
        
    return rG, vG


#p63, p64 guidance algorithm
def guide(T, TG, rT, rG, vG, k, inputs):
    #generate guidance coordinate axes
    xaxis = rT / np.linalg.norm(rT)
    yaxis = (np.cross(rT, rG - k * (vG - np.cross([0, 0, inputs['bodyrotation']], rG)) * T / 4)) / np.linalg.norm(np.cross(rT, rG - k * (vG - np.cross([0, 0, inputs['bodyrotation']], rG)) * T / 4))
    zaxis = np.cross(xaxis, yaxis) / np.linalg.norm(np.cross(xaxis, yaxis))     
    q = np.array([xaxis, yaxis, zaxis])
    
    #calculation state vectors in guidance coordinates
    RG = np.matmul(q, rG - rT)
    VG = np.matmul(q, vG - np.cross([0, 0, inputs['bodyrotation']], rG))
  
    #calculation T
    a = TG[3][2]
    b = 6 * TG[2][2]
    c = 18 * TG[1][2] + 6 * VG[2]
    d = 24 * (TG[0][2] - RG[2])
    T += inputs['dt']
    while True:
        dT = -(a * T ** 3 + b * T ** 2 + c * T + d) / (3 * a * T ** 2 + 2 * b * T + c)
        T = T + dT
        if abs(dT) <  1 / 128:
            break
 
    #calculate required acceleration
    TP = T + inputs['LEADTIME']
    a = (3 * (TP / T) ** 2 - 2 * (TP / T)) * 12 / T ** 2
    b = (4 * (TP / T) ** 2 - 3 * (TP / T))  * 6 / T
    c = (2 * (TP / T) ** 2 - (TP / T)) * 6 / T
    d = (6 * (TP / T) ** 2 - 6 * (TP / T) + 1)
    ACG1 = a * (TG[0] - RG)
    ACG2 = b * TG[1]
    ACG3 = c * VG
    ACG4 = d * TG[2]
    ACG = ACG1 + ACG2 + ACG3 + ACG4
    GP = -inputs['gravparam'] / np.linalg.norm(rG) ** 3 * rG
    AFCP = np.matmul(np.transpose(q), ACG) - GP  
    
    return AFCP, T, RG, VG

   
def target(inputs):   
    #equation 21.2
    def EQ2102(TIF, TMF, m3):
        m1row1 = np.array([TMF  ** 2 / 2, TMF  ** 3 / 6, TMF  ** 4 / 24])
        m1row2 = np.array([TMF, TMF  ** 2 / 2, TMF  ** 3 / 6])
        m1row3 = np.array([TIF  ** 2 / 2, TIF  ** 3 / 6, TIF  ** 4 / 24])
        m1 = np.array([m1row1, m1row2, m1row3])
        m1 = np.linalg.inv(m1)
        m2row1 = np.array([-1, -TMF, 1, 0, 0])
        m2row2 = np.array([0, -1, 0, 1, 0])
        m2row3 = np.array([-1, -TIF, 0, 0, 1])
        m2 = np.array([m2row1, m2row2, m2row3])      
        APFGX = np.matmul(np.matmul(m1, m2), m3)
        APFGX = np.array([m3[0], m3[1], APFGX[0], APFGX[1], APFGX[2]])
        return APFGX
    
    #equations 21.4 - 21.6
    def EQ2104(TIF, TMF, TAU, m2):
        m1row1 = np.array([TAU ** 2 - TAU * TMF + TMF ** 2 / 2, TMF ** 3 / 6, TMF ** 4 / 24])
        m1row2 = np.array([-TAU + TMF, TMF ** 2 / 2, TMF ** 3 / 6])
        m1row3 = np.array([TAU ** 2 - TAU * TIF + TIF ** 2 / 2, TIF ** 3 / 6, TIF ** 4 / 24])
        m1 = np.array([m1row1, m1row2, m1row3])
        m1 = np.linalg.inv(m1)
        APFGZ  = np.matmul(m1, m2)
        APFGZ = np.array([APFGZ[0] * TAU ** 2, -APFGZ[0] * TAU, APFGZ[0], APFGZ[1], APFGZ[2]])
        return APFGZ
    
    #calculation of APTG and APIG
    #equation 21.2
    TMF = inputs['TAPM'] - inputs['TAPF']
    TIF = inputs['TAPI'] - inputs['TAPF']
    APFGX = EQ2102(TIF, TMF, [inputs['RAPFGX'], inputs['VAPFGX'], inputs['RAPMG'][0], inputs['VAPMG'][0], inputs['RAPIG'][0]])
    APFGZ = EQ2104(TIF, TMF, inputs['TAU'], [inputs['RAPMG'][2], inputs['VAPMG'][2], inputs['RAPIG'][2]])
    APFG = []
    for i in range(5):
        APFG.append([APFGX[i], 0, APFGZ[i]])
    APFG = np.vstack(APFG)
    
    #equation 21.7
    APTG = np.matmul(STM(0, inputs['TAPF']), APFG)
    APIG = np.matmul(STM(inputs['TAPI'], 0), APTG)
    
    #calculation of BRTG, BRFG, and BRIG
    #equations 22.1-22.6
    RBRFG = APIG[0]
    VBRFG = APIG[1]
    ABRFG = [np.cos(inputs['PBRF']) * inputs['FBRF'] / inputs['MBRF'] - inputs['surfacegravity'], 0, -np.sin(inputs['PBRF']) * inputs['FBRF'] / inputs['MBRF']]
    MDBRF = -inputs['FBRF'] / inputs['exhaustvelocity']
    #if not first pass use calculated values from equations 22.9 - 22.14
    if inputs['JBRFGFLAG'] == 0:
        JBRFG = [0, 0, -inputs['KJ'] * ABRFG[2] * MDBRF / inputs['MBRF']]
        SBRFG = [0, 0, 0]
    else:
        JBRFG = [inputs['JBRFG'][0], 0, -inputs['KJ'] * ABRFG[2] * MDBRF / inputs['MBRF']]
        SBRFG = [inputs['SBRFG'][0], 0, inputs['SBRFG'][2]]

    BRFG = [RBRFG, VBRFG, ABRFG, JBRFG, SBRFG]
    BRTG = np.matmul(STM(0, inputs['TBRF']), BRFG)
    BRIG = np.matmul(STM(inputs['TBRI'], 0), BRTG)
    
    #if not first pass user calculated value from equation 21.8
    if inputs['RBRIGZ'] != None:
        BRIG[0][2] = inputs['RBRIGZ']   
       
    return APIG, APTG, BRIG, BRTG, BRFG


def initialize(inputs, initialstate):
    #generate initial state and targets
    state = np.copy(initialstate)
    rTi = inputs['bodyradius'] * np.array([np.cos(inputs['theta']) * np.sin(inputs['psi']), np.sin(inputs['theta']) * np.sin(inputs['psi']), np.cos(inputs['psi'])])
    rT = rodrigues(rTi, np.array([0, 0, 1]), inputs['bodyrotation'] * inputs['t'])
    APIG, APTG, BRIG, BRTG, BRFG = target(inputs)

   #intigrate trajectory until T is greater than TBRI     
    while True:
        AFCP, T, RG, VG = guide(inputs['GUIDETIME'], BRTG, rT, state[0], state[1], 1, inputs)
        if T > inputs['TBRI']:
            break
        else:
            inputs['t'] += inputs['dt']
            state = RK4(state[0], state[1], 0, inputs['dt'], inputs['gravparam'])
            rT = rodrigues(rTi, np.array([0, 0, 1]), inputs['bodyrotation'] * inputs['t'])
        
    #p63 ignition algorithm                    
    while True:
        #equations 11.6 - 11.10
        for i in range(4):
            RP, VP = cse(initialstate[0], initialstate[1], inputs['t'], inputs['gravparam'])
            UNFCP = AFCP / np.linalg.norm(AFCP)
            VP += UNFCP * inputs['AFTRIM'] * inputs['TTRIM']
            rT = rodrigues(rTi, np.array([0, 0, 1]), inputs['bodyrotation'] * inputs['t'])
            AFCP, T, RG, VG = guide(inputs['GUIDETIME'], BRTG, rT, RP, VP, 1, inputs)
        
        #equations 11.11 - 11.12
        KXYV = inputs['KXYV']                
        DRGX = KXYV[0] * (RG[0] - BRIG[0][0])
        DRGY = KXYV[1] * RG[1] ** 2
        DRGZ = (RG[2] - BRIG[0][2])
        DRGV = KXYV[2] * (np.linalg.norm(VG) - np.linalg.norm(BRIG[1]))
        DEN = VG[2] + KXYV[0] * VG[0]  
        DGUIDETIME = (DRGX + DRGY + DRGZ + DRGV) / DEN     
        inputs['t'] -= DGUIDETIME
        inputs['GUIDETIME'] -= DGUIDETIME
        if abs(DGUIDETIME) < 1 / 128:
            break
    #equations 11.13
    inputs['GUIDETIME'] -= inputs['TTRIM']
    inputs['t'] -= inputs['TTRIM']
    inputs['ti'] = inputs['t']
    return APTG, BRFG, BRTG, BRIG


def simulation(inputs):
    #equation 21.9
    def EQ2109(T, m2):
        m1row1 = np.array([-24 / T ** 3, -18 / T ** 2, - 6 / T, 24 / T ** 3, -6 / T ** 2])
        m1row2 = np.array([72 / T ** 4, 48 / T ** 3, 12 / T ** 2, -72 / T ** 4, 24 / T ** 3])
        m1 = np.array([m1row1, m1row2])
        return np.matmul(m1, m2)

    #generate initial state and calculate targets
    initialstate = statevector(inputs)
    APTG, BRFG, BRTG, BRIG = initialize(inputs, initialstate)
    state = cse(initialstate[0], initialstate[1], inputs['t'], inputs['gravparam'])  
    rTi = inputs['bodyradius'] * np.array([np.cos(inputs['theta']) * np.sin(inputs['psi']), np.sin(inputs['theta']) * np.sin(inputs['psi']), np.cos(inputs['psi'])])

    #set simulation variables
    T = inputs['GUIDETIME']
    TG = BRTG
    TTHROTAFLAG = 0
    TTHROTA = 0
    dtinc = 3
    
    #main simulation loop    
    while T < inputs['TBRF']:
        #reduce dt when approaching TBRF
        if inputs['TBRF'] - T <= 2 * inputs['dt'] and dtinc > 0:
            inputs['dt'] /= 10
            dtinc -= 1
            
        #calculate body rotation, desired acceleration, and maximum acceleration
        rT = rodrigues(rTi, np.array([0, 0, 1]), inputs['bodyrotation'] * inputs['t'])      
        AFCP, T, RG, VG = guide(T, TG, rT, state[0], state[1], 1, inputs)
        maxAFCP = inputs['maxthrust'] / inputs['mass']
        
        #limit acceleration
        if inputs['t'] - inputs['ti'] < 26:
            maxAFCP = inputs['AFTRIM']
        if np.linalg.norm(AFCP) > maxAFCP * inputs['minthrott']:
            AFCP *= maxAFCP / np.linalg.norm(AFCP)
        #record time of throttle recovery
        elif TTHROTAFLAG == 0:
            TTHROTA = T
            TTHROTAFLAG = 1
        #calculate state after dt
        state = RK4(state[0], state[1], AFCP, inputs['dt'], inputs['gravparam'])
        FCP = AFCP * inputs['mass']
        inputs['mass'] -= np.linalg.norm(FCP) / inputs['exhaustvelocity'] * inputs['dt']
        inputs['t'] += inputs['dt']

    #calculate iterative values
    RBRIGZ = BRIG[0][2] + inputs['KIG'] * (TTHROTA - inputs['TTHROT'])
    MASSRATIO = np.exp((np.linalg.norm(BRFG[1]) - np.linalg.norm(VG)) / inputs['exhaustvelocity'])
    PACK = EQ2109(T,  [BRTG[0], BRTG[1], BRTG[2], RG, VG])
    JBRTGA = PACK[0]
    SBRTGA = PACK[1]       
    BRFGA = np.matmul(STM(inputs['TBRF'], 0), [BRTG[0], BRTG[1], BRTG[2], [JBRTGA[0], JBRTGA[1], BRTG[3][2]], SBRTGA])
     
    return BRFG, BRFGA, RBRIGZ, MASSRATIO, TTHROTA
    

def fullsimulation(inputs):
    #generate initial state and calculate targets
    initialstate = statevector(inputs)
    APTG, BRFG, BRTG, BRIG = initialize(inputs, initialstate)
    state = cse(initialstate[0], initialstate[1], inputs['t'], inputs['gravparam'])
    rTi = inputs['bodyradius'] * np.array([np.cos(inputs['theta']) * np.sin(inputs['psi']), np.sin(inputs['theta']) * np.sin(inputs['psi']), np.cos(inputs['psi'])])
    
    #set simulation variables
    T = inputs['GUIDETIME']
    TG = BRTG
    k = 1
    APFLAG = 0
    statelist = []
       
    #main simulation loop
    while T < inputs['TAPF']:
        #switch to p64
        if T > inputs['TBRF'] and APFLAG == 0:
            TG = APTG
            T = inputs['TAPI']
            k = 0
            APFLAG = 1
        
        #calculate body rotation, desired acceleration, and maximum acceleration
        rT = rodrigues(rTi, np.array([0, 0, 1]), inputs['bodyrotation'] * inputs['t'])     
        AFCP, T, RG, VG = guide(T, TG, rT, state[0], state[1], k, inputs)
        maxAFCP = inputs['maxthrust'] / inputs['mass']
         
        #limit acceleration
        if inputs['t'] - inputs['ti'] < 26:
            maxAFCP = inputs['AFTRIM']
        if np.linalg.norm(AFCP) > maxAFCP * inputs['minthrott']:
            AFCP *= maxAFCP / np.linalg.norm(AFCP)
        
        #calculate state after dt
        state = RK4(state[0], state[1], AFCP, inputs['dt'], inputs['gravparam'])
        FCP = AFCP * inputs['mass']
        inputs['mass'] -= np.linalg.norm(FCP) / inputs['exhaustvelocity'] * inputs['dt']
        inputs['t'] += inputs['dt']
        
        #output loop 1
        outputlist = [[T], rT.tolist(), AFCP.tolist(), RG.tolist(), VG.tolist(), state[0].tolist(), state[1].tolist(),[inputs['mass']]]
        flatoutputlist = []
        for sublist in outputlist:
            for item in sublist:
                flatoutputlist.append(item)
        statelist.append(flatoutputlist)
        
    #output loop 2
    file = open('output.txt', 'w')
    for i in range(len(statelist)):
        file.write(str(statelist[i]).strip("[]") + '\n')      
    file.close

    return 0


#conic state extrapolation
def cse(r0, v0, dt, gm):
    def scfunc(z):
        az = abs(z)
        if az < 1e-4:
            return (1 - z * (0.05 - z / 840)) / 6, 0.5 - z * (1 - z / 30) / 24
        else:
            saz = np.sqrt(az)
            if z > 0:
                x = saz
                return (saz - np.sin(x)) / (saz * az), (1 - np.cos(x)) / az
            else:
                x = np.exp(saz)
                return (0.5 * (x - 1 / x) - saz) / (saz * az), (0.5 * (x + 1 / x) - 1) / az

    rscale = np.linalg.norm(r0)
    vscale = np.sqrt(gm / rscale)
    r0s = r0 / rscale
    v0s = v0 / vscale
    dts = dt * vscale / rscale
    v2s = np.linalg.norm(v0) ** 2 * rscale / gm
    alpha = 2 - v2s
    armd1 = v2s - 1
    rvr0s = np.dot(r0, v0) / np.sqrt(gm * rscale)

    x = 0
    ratio = 1
    x2 = x * x
    z = alpha * x2
    s, c = scfunc(z)
    x2c = x2 * c
    f = 0
    df = 0

    while abs(ratio) > 5e-9:
        f = x + rvr0s * x2c + armd1 * x * x2 * s - dts
        df = x * rvr0s * (1 - z * s) + armd1 * x2c + 1
        ratio = f / df
        x = x - ratio
        x2 = x * x
        z = alpha * x2
        s, c = scfunc(z)
        x2c = x2 * c

    lf = 1 - x2c
    lg = dts - x2 * x * s

    r1 = lf * r0s + lg * v0s
    ir1 = 1 / np.linalg.norm(r1)
    lfdot = ir1 * x * (z * s - 1)
    lgdot = 1 - x2c * ir1

    v1 = lfdot * r0s + lgdot * v0s

    return r1 * rscale, v1 * vscale


#plotting function
def plot(APTG, BRTG, TAPI, TBRI):
    def STM(t1, t0):
        dt = t1 - t0
        xlist = []
        ylist = []
        for i in range((-TAPI + 1)):
            for j in range(5):
                if j >= i:
                    yval = j - i
                    xval = dt ** yval / factorial(yval)
                else:
                    xval = 0
                xlist.append(xval)
                xval = 0
            ylist.append(xlist)
            xlist = []
        return np.vstack(ylist)
     
    plot1x = []
    plot1y = []
    plot2x = []
    plot2y = []
      
    val1 = 1
    val2 = 2
    
    for i in range(160):
        t1 = -i + 0
        x = np.matmul(STM(t1, 0), APTG)
        plot1x.append(x[0][2])
        plot1y.append(x[val1][val2])
    
    
#    for i in range(-TBRI + 1):
    for i in range((-TAPI + 1)  +60):
        t1 = -i
        x = np.matmul(STM(t1, 0), BRTG)
        plot2x.append(x[0][2])
        plot2y.append(x[val1][val2])
      
    f = plt.figure(1)
    plt.scatter(plot1x,plot1y, label='-', color='k', s=25, marker="o")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('p63')
    plt.legend()
        
    f = plt.figure(2)
    
    plt.scatter(plot2x,plot2y, label='-', color='r', s=25, marker="o")
    plt.scatter(plot1x,plot1y, label='-', color='k', s=25, marker="o")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('p64')
    plt.legend()
    plt.show()
