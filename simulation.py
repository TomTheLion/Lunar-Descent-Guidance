import numpy as np
import util as util

#inputs_reserved = {
#        #body parameters
#        'gravparam': 4.90389580165072e12,
#        'surfacegravity': 1.622782,
#        'bodyradius': 1744363,
#        'bodyrotation': 2.66627e-6,
#        #target parameters
#        'theta': 0.409680,
#        'psi': 1.559044,
#        #orbit parameters
#        'apoapsis': 110000,
#        'periapsis': 15000,
#        'inclination': 0.017628,
#        'trueanomaly': 5.95,
#        'argumentofperiapsis': None,
#        'longitudeofascendingnode': None,
#        #ship parameters
#        'mass': 5764,
#        'exhaustvelocity': 320 * 9.807,
#        'maxthrust': 0.60 * 30000,
#        'AFTRIM': 0.10 * 0.60 * 30000 / 5764,
#        'TTRIM': 26,
#        #descent parameters
#        'GUIDETIME': -664.4,
#        'LEADTIME': 2,
#        'TAPF': -10,
#        'TAPM': -60,
#        'TAPI': -160,
#        'TBRF': -60,
#        'TBRI': -660,
#        'TAU': 8,
#        'RAPFGX': 30,
#        'VAPFGX': -1,
#        'RAPMG': [150, 0, -150],
#        'VAPMG': [-7, 0, 7],
#        'RAPIG': [2200, 0, -7500],
#        'FBRF': 0.5 * 0.60 * 30000,
#        'MBRF': 3400,
#        'PBRF': 0.958,
#        'KJ': 1.2,
#        'KXYV': [0.2, 1, 1],
#        'TTHROT': -630,
#        'KIG': -500,
#        'JBRFGFLAG' : 0,
#        'RBRIGZ': None,
#        'JBRFG': None,
#        'SBRFG': [0, 0, 0],
#        #simulation parameters
#        'dt': 1,
#        't': 0,
#        'ti': 0
#        }

inputs_reserved = {
        #body parameters
        'gravparam': 4.90389580165072e12,
        'surfacegravity': 1.626,
        'bodyradius': 1744363,
        'bodyrotation': 2.66627e-6,
        #target parameters
        'theta': 0.409680,
        'psi': 1.559044,
        #orbit parameters
        'apoapsis': 110000,
        'periapsis': 15000,
        'inclination': np.pi - 0.017628,
        'trueanomaly': 5.95,
        'argumentofperiapsis': None,
        'longitudeofascendingnode': None,
        #ship parameters
        'mass': 14861,
        'exhaustvelocity': 311 * 9.802,
        'maxthrust': 43900,
        'AFTRIM': 0.111 * 43900 / 14861,
        'minthrott': 0.57,
        'TTRIM': 26,
        #descent parameters
        'GUIDETIME': -664.4,
        'LEADTIME': 2,
        'TAPF': -10,
        'TAPM': -60,
        'TAPI': -170,
        'TBRF': -60,
        'TBRI': -760,
        'TAU': 8,
        'RAPFGX': 30,
        'VAPFGX': -1,
        'RAPMG': [150, 0, -523],
        'VAPMG': [-5, 0, 17.4],
        'RAPIG': [2200, 0, -7500],
        'FBRF': 0.50 * 43900,
        'MBRF': 8545,
        'PBRF': 1.05,
        'KJ': 1.2,
        'KXYV': [0.2, 0, 0],
        'TTHROT': -280,
        'KIG': -100,
        'JBRFGFLAG' : 0,
        'RBRIGZ': None,
        'JBRFG': None,
        'SBRFG': [0, 0, 0],
        #simulation parameters
        'dt': 1,
        't': 0,
        'ti': 0
}

for run in range(100):
    inputs = inputs_reserved.copy()
    BRFG, BRFGA, RBRIGZ, MASSRATIO, TTHROTA = util.simulation(inputs)

    inputs_reserved['RBRIGZ'] = RBRIGZ
    inputs_reserved['JBRFG'] = BRFGA[3]
    inputs_reserved['SBRFG'] = np.array(BRFGA[4]) * 0.5 + np.array(inputs_reserved['SBRFG']) * 0.5
    inputs_reserved['MBRF'] = 0.5 * inputs_reserved['MBRF'] * MASSRATIO + 0.5 * inputs_reserved['MBRF']
    inputs_reserved['JBRFGFLAG'] = 1

    if abs(TTHROTA - inputs['TTHROT']) < 2:
        inputs_reserved['KIG'] = -10

    if abs(TTHROTA - inputs['TTHROT']) < 1:
        inputs_reserved['KIG'] = -1

    error1 = '{:.2e}'.format(abs(BRFGA[2][2] - BRFG[2][2]) / np.linalg.norm(BRFG[2]))
    error2 = '{:.2e}'.format(abs(BRFGA[2][0] - BRFG[2][0]) / np.linalg.norm(BRFG[2]))
    error3 = '{:.2e}'.format(abs(BRFGA[1][0] - BRFG[1][0]) / np.linalg.norm(BRFG[1]))
    error4 = '{:.2e}'.format(abs(TTHROTA - inputs['TTHROT']))
    error5 = '{:.2e}'.format(abs(inputs['mass'] - inputs['MBRF']))

    print(error1, error2, round(BRFGA[1][0], 2), round(BRFG[1][0], 2), round(TTHROTA, 1), round(inputs['TTHROT'], 1),
          inputs['mass'], inputs['MBRF'], BRFGA[3][2], RBRIGZ, BRFGA[4])

APIG, APTG, BRIG, BRTG, BRFG = util.target(inputs_reserved)

print('RBRIGZ: ' + str(inputs_reserved['RBRIGZ']))
print('JBRFG: ' + str(inputs_reserved['JBRFG']))
print('SBRFG: ' + str(inputs_reserved['SBRFG']))
print('MBRF: ' + str(inputs_reserved['MBRF']))
util.aoplan(inputs_reserved)

util.plot(APTG, BRTG, inputs_reserved['TAPI'], inputs_reserved['TBRI'])

inputs = inputs_reserved.copy()
util.fullsimulation(inputs)  














