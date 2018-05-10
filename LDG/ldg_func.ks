//P64 Approach phase targeting routine
FUNCTION APTR {

	//Equation 21.1
	LOCAL TMF IS inputs["TAPM"] - inputs["TAPF"].
	LOCAL TIF IS inputs["TAPI"] - inputs["TAPF"].
		
	//Equation 21.2: Terminal x components
	LOCAL APFGX IS EQ2102(TMF, TIF, LIST(inputs["RAPFGX"], inputs["VAPFGX"], inputs["RAPMG"][0], inputs["VAPMG"][0], inputs["RAPIG"][0])).
	
	//Equation 21.4 - 21.6: Terminal z components
	LOCAL APFGZ IS EQ2104(inputs["TAU"], TMF, TIF, LIST(inputs["RAPMG"][2], inputs["VAPMG"][2], inputs["RAPIG"][2])).

	//Equation 21.7: Compute target and initial states
	LOCAL APFG IS LIST().
	FROM { LOCAL i IS 0. } UNTIL i = 5 STEP { SET i TO i + 1. } DO {
		APFG:ADD(LIST(APFGX[i], 0, APFGZ[i])).
	}			
	LOCAL APTG IS MMUL(STM(0, inputs["TAPF"]), APFG).
	LOCAL APIG IS MMUL(STM(inputs["TAPI"], 0), APTG).
		
	RETURN LIST(APTG, APIG).
}

//P63 Braking phase targeting routine
FUNCTION BRTR {
	
	LOCAL RBRFG IS APIG[0].
	LOCAL VBRFG IS APIG[1].
	LOCAL UNFBRFG IS LIST(COS(inputs["PBRF"]), 0, -SIN(inputs["PBRF"])).
	LOCAL ABRFG IS LIST((inputs["FBRF"] / inputs["MBRF"]) * UNFBRFG[0] - inputs["GM"], 0, (inputs["FBRF"] / inputs["MBRF"]) * UNFBRFG[2]).
	LOCAL MDBRF IS inputs["FBRF"] / inputs["VEXBRF"].
	LOCAL JBRFG IS LIST(inputs["JBRFG"][0], 0, inputs["KJ"] * ABRFG[2] * MDBRF / inputs["MBRF"]).
	LOCAL SBRFG IS LIST(inputs["SBRFG"][0], 0, inputs["SBRFG"][2]).
	
	LOCAL BRFG IS LIST(RBRFG, VBRFG, ABRFG, JBRFG, SBRFG).		
	LOCAL BRTG IS MMUL(STM(0, inputs["TBRF"]), BRFG).
	LOCAL BRIG IS MMUL(STM(inputs["TBRI"], 0), BRTG).
	
	SET BRIG[0][2] TO inputs["RBRIGZ"].
		
	RETURN LIST(BRTG, BRIG).
}

//P63 ignition algorithm
FUNCTION P63I {
	
	LOCAL T IS inputs["TBRI"].
	LOCAL RTP IS rodrigues(inputs["RTPI"], V(0, 1, 0), CONSTANT:RADTODEG * BODY:ANGULARVEL:Y * (clocktime - starttime)).
	LOCAL RP IS inputs["RPI"].
	LOCAL VP IS inputs["VPI"].

	UNTIL 0 {
		LOCAL PACK IS GUIDE(T, BRTG, RTP, RP, VP, 1, 1).
		SET T TO PACK[3].
		IF T > inputs["TBRI"] {
			BREAK.
		} ELSE {
			SET clocktime TO clocktime + 10.
			SET PACK TO RK4(RP, VP, V(0, 0, 0), 10).
			SET RP TO PACK[0].
			SET VP TO PACK[1].
			SET RTP TO rodrigues(inputs["RTPI"], V(0, 1, 0), CONSTANT:RADTODEG * BODY:ANGULARVEL:Y * (clocktime - starttime)).
			}
		}
		
	UNTIL 0 {
		LOCAL RG IS 0.
		LOCAL VG IS 0.
		LOCAL RP IS inputs["RPI"].
		LOCAL VP IS inputs["VPI"].
		LOCAL PACK IS cse(RP, VP, clocktime - starttime).
		SET RP0 TO PACK[0].
		SET VP0 TO PACK[1].
		FROM { LOCAL I IS 0. } UNTIL I = 4 STEP { SET I TO I + 1. } DO {

			SET RTP TO rodrigues(inputs["RTPI"], V(0, 1, 0), CONSTANT:RADTODEG * BODY:ANGULARVEL:Y * (clocktime - starttime)).
			LOCAL PACK IS GUIDE(T, BRTG, RTP, RP0, VP0, 1, 0).
			LOCAL UNFCP IS PACK[0]:NORMALIZED.
			SET RG TO PACK[1].
			SET VG TO PACK[2].
			LOCAL T IS PACK[3].
			SET VP TO VP0 + UNFCP * inputs["AFTRIM"] * inputs["TTRIM"] / MASS / 1000.
		}
		LOCAL KXYV IS inputs["KXYV"].              
        LOCAL DRGX IS KXYV[0] * (RG[0] - BRIG[0][0]).
        LOCAL DRGY IS KXYV[1] * RG[1] ^ 2.
        LOCAL DRGZ IS (RG[2] - BRIG[0][2]).
        LOCAL DRGV IS KXYV[2] * (MMAG(VG) - MMAG(BRIG[1])).
        LOCAL DEN IS VG[2] + KXYV[0] * VG[0].  
        LOCAL DGUIDETIME IS (DRGX + DRGY + DRGZ + DRGV) / DEN.
        SET clocktime TO clocktime - DGUIDETIME.
        SET T TO T - DGUIDETIME.
        IF ABS(DGUIDETIME) < 1 / 128 {
            BREAK.
		}	
	}
	
	RETURN T.
}

//P63, P64 Guidance algorithm
FUNCTION GUIDE{
	PARAMETER T.
	PARAMETER TG.
	PARAMETER RTP.
	PARAMETER RP.
	PARAMETER VP.
	PARAMETER k.
	PARAMETER deltaclocktime.
	
	//Creating guidance coordinate axes
	LOCAL XAXIS IS RTP:NORMALIZED.
	LOCAL YAXIS IS VCRS(RTP, RP - k * (VP - VCRS(BODY:ANGULARVEL, RP)) * T / 4):NORMALIZED.
	LOCAL ZAXIS IS VCRS(XAXIS, YAXIS).
	
	//Creating coordinate transform matrices
	LOCAL qrow1 IS LIST(XAXIS:X, XAXIS:Z, XAXIS:Y).
	LOCAL qrow2 IS LIST(YAXIS:X, YAXIS:Z, YAXIS:Y).
	LOCAL qrow3 IS LIST(ZAXIS:X, ZAXIS:Z, ZAXIS:Y).	
	LOCAL qmatrix IS LIST(qrow1, qrow2, qrow3).
	LOCAL Tqmatrix IS MTRNS(qmatrix).

	//Calculation of position and velocity in guidance coordinates
	LOCAL RG IS RP - RTP.
	LOCAL VG IS VP - VCRS(BODY:ANGULARVEL, RP).
	LOCAL RG IS MMUL(qmatrix, LIST(RG:X, RG:Z, RG:Y)).
	LOCAL VG IS MMUL(qmatrix, LIST(VG:X, VG:Z, VG:Y)).

	//Equation 6.10 - 6.13: Calculation of T
	SET T TO T + deltaclocktime.
	LOCAL a IS TG[3][2].
	LOCAL b IS 6 * TG[2][2].
	LOCAL c IS 18 * TG[1][2] + 6 * VG[2].
	LOCAL d IS 24 * (TG[0][2] - RG[2]).
	LOCAL deltaT IS 1.
	UNTIL ABS(deltaT) <  1 / 128 {
		SET deltaT TO -(a * T ^ 3 + b * T ^ 2 + c * T + d) / (3 * a * T ^ 2 + 2 * b * T + c).
		SET T TO T + deltaT.
	}

	//Equation 6.14 - 6.15: calculation of ACG
	LOCAL TP IS T + inputs["LEADTIME"].
	LOCAL a IS (3 * (TP / T) ^ 2 - 2 * (TP / T)) * 12 / T ^ 2.
	LOCAL b IS (4 * (TP / T) ^ 2 - 3 * (TP / T))  * 6 / T.
	LOCAL c IS (2 * (TP / T) ^ 2 - (TP / T)) * 6 / T.
	LOCAL d IS (6 * (TP / T) ^ 2 - 6 * (TP / T) + 1).
	LOCAL ACG1 IS SMUL(a, MSUB(TG[0], RG)).
	LOCAL ACG2 IS SMUL(b, TG[1]).
	LOCAL ACG3 IS SMUL(c, VG).
	LOCAL ACG4 IS SMUL(d, TG[2]).
	LOCAL ACG IS MADD(MADD(MADD(ACG1, ACG2), ACG3), ACG4).

	//Calculation of FCP
	LOCAL GP IS BODY:MU / BODY:POSITION:MAG ^ 2 * BODY:POSITION:NORMALIZED.
	LOCAL AFCP IS MMUL(Tqmatrix, ACG).
	SET AFCP TO V(AFCP[0], AFCP[2], AFCP[1]) - GP.
	LOCAL FCP IS AFCP * SHIP:MASS * 1000.

	RETURN LIST(FCP, RG, VG, T).
}

//P66 Terminal Descent
FUNCTION TGUIDE {
	UNTIL ALT:RADAR < inputs["offset"] + 1.7 {
		SET GP TO BODY:MU / BODY:POSITION:MAG ^ 2 * BODY:POSITION:NORMALIZED.
		SET AFYZ TO VXCL(UP:VECTOR, -VELOCITY:SURFACE) / 10.
		IF AFYZ:MAG > 0.35 * GP:MAG { SET AFYZ TO 0.35 * GP:MAG * AFYZ:NORMALIZED. }		
		SET AFX TO (((-2 - VERTICALSPEED) / 1.5) * UP:VECTOR - GP) / COS(VANG(UP:VECTOR, FACING:VECTOR)).
		SET AFCP TO AFYZ + AFX.

		LOCK STEERING TO AFCP:NORMALIZED.
		LOCK THROTTLE TO THROTT(AFCP:MAG * MASS * 1000).
		WAIT 0.001.
	}
}

//Equation 21.2
FUNCTION EQ2102 {
	PARAMETER TMF.
	PARAMETER TIF.
	PARAMETER m3.
	
	LOCAL m1row1 IS LIST(TMF ^ 2 / 2, TMF ^ 3 / 6, TMF ^ 4 / 24).
	LOCAL m1row2 IS LIST(TMF, TMF ^ 2 / 2, TMF ^ 3 / 6).
	LOCAL m1row3 IS LIST(TIF ^ 2 / 2, TIF ^ 3 / 6, TIF ^ 4 / 24).
	LOCAL m1 IS LIST(m1row1, m1row2, m1row3).
	LOCAL m1 IS MINV(m1).
	
	LOCAL m2row1 IS LIST(-1, -TMF, 1, 0, 0).
	LOCAL m2row2 IS LIST(0, -1, 0, 1, 0).
	LOCAL m2row3 IS LIST(-1, -TIF, 0, 0, 1).
	LOCAL m2 IS LIST(m2row1, m2row2, m2row3).
		
	LOCAL APFGX IS MMUL(MMUL(m1, m2), m3).
	APFGX:INSERT(0, inputs["VAPFGX"]).
	APFGX:INSERT(0, inputs["RAPFGX"]).
	
	RETURN APFGX.
}

//Equation 21.4
FUNCTION EQ2104 {
	PARAMETER TAU.
	PARAMETER TMF.
	PARAMETER TIF.
	PARAMETER m2.

	LOCAL m1row1 IS LIST(TAU ^ 2 - TAU * TMF + TMF ^ 2 / 2, TMF ^ 3 / 6, TMF ^ 4 / 24).
	LOCAL m1row2 IS LIST(-TAU + TMF, TMF ^ 2 / 2, TMF ^ 3 / 6).
	LOCAL m1row3 IS LIST(TAU ^ 2 - TAU * TIF + TIF ^ 2 / 2, TIF ^ 3 / 6, TIF ^ 4 / 24).
	LOCAL m1 IS LIST(m1row1, m1row2, m1row3).
	LOCAL m1 IS MINV(m1).
	
	LOCAL APFGZ IS MMUL(m1, m2).
	APFGZ:INSERT(0, -APFGZ[0] * TAU).
	APFGZ:INSERT(0, APFGZ[1] * TAU ^ 2).
	
	RETURN APFGZ.
}

//State transition matrix
FUNCTION STM {
	PARAMETER t1.
	PARAMETER t0.
	
	LOCAL dt IS t1 - t0.
	LOCAL xlist IS LIST().
	LOCAL ylist IS LIST().
	
	FROM { LOCAL i IS 0. } UNTIL i = 5 STEP { SET i TO i + 1. } DO {
		FROM { LOCAL j IS 0. } UNTIL j = 5 STEP { SET j TO j + 1. } DO {
			IF j >= i {
				SET yval TO j - i.
				SET xval TO dt ^ yval / FACT(yval).
			} ELSE {
				SET xval TO 0.
			}
			xlist:ADD(xval).
			SET xval TO 0.
		}
		ylist:ADD(xlist).
		SET xlist TO LIST().
	}

	RETURN ylist.
}

//State extrapolation for P63 ignition algorithm, runge kutta
FUNCTION RK4 {
	PARAMETER r0.
	PARAMETER v0.
	PARAMETER a.
	PARAMETER h.
	
	FUNCTION f {
		PARAMETER r.		
		LOCAL g IS -BODY:MU / r:MAG ^ 3 * r.
		RETURN a + g.
	}

	LOCAL kv1 IS f(r0).
    LOCAL kr1 IS v0.
    LOCAL kv2 IS f(r0 + kr1 * h / 2).
    LOCAL kr2 IS v0 + kv1 * h / 2.
    LOCAL kv3 IS f(r0 + kr2 * h / 2).
    LOCAL kr3 IS v0 + kv2 * h / 2.
    LOCAL kv4 IS f(r0 + kr3 * h).
    LOCAL kr4 IS v0 + kv3 * h.   
    LOCAL v1 IS v0 + h / 6 * (kv1 + 2 * kv2 + 2 * kv3 + kv4).
    LOCAL r1 IS r0 + h / 6 * (kr1 + 2 * kr2 + 2 * kr3 + kr4).
    RETURN LIST(r1, v1).	
}

FUNCTION THROTT {
	PARAMETER REQTHRUST.
	
	LOCAL SETTHROTT IS 1.
	
	IF THROTTLECONTROL = 1 {
		IF REQTHRUST > inputs["FMID"] { SET REQTHRUST TO inputs["FMID"]. }
		SET SETTHROTT TO (REQTHRUST - inputs["FMIN"]) / (inputs["FMAX"] - inputs["FMIN"]).
	}
			
	RETURN SETTHROTT.
}







