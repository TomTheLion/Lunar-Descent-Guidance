CLEARSCREEN.
SET CONFIG:IPU TO 1000.

GLOBAL starttime IS TIME:SECONDS.
GLOBAL clocktime IS starttime.
GLOBAL clocktimeOLD IS starttime.

GLOBAL mission IS LEXICON(
					"latitude", 0.67337,
					"longitude", 23.47293,
					"inclination", 1.01					
).

GLOBAL LS IS LATLNG(mission["latitude"], mission["longitude"]).

GLOBAL inputs IS LEXICON(
					"GM", 1.626,
					"LEADTIME", 2,
					"OFFSET", 5,
					"TAPF", -10,
					"TAPM", -60, //-60
					"TAPI", -170,
					"TBRF", -60,
					"TBRI", -760,										
					"TAU", 8,
					"RAPFGX", 30,
					"VAPFGX", -1,		
					"RAPMG", LIST(150, 0, -523),
					"VAPMG", LIST(-5, 0, 17.4),	//7				
					"RAPIG", LIST(2200, 0, -7500),
					"FMAX", 43900,
					"FMID", 25023,
					"FMIN", 4670,
					"FBRF", 0.5 * 43900,
					"VEXBRF", 311 * 9.802,
					"AFTRIM", 0.111 * 43900,
					"TTRIM", 26,
					"MBRF", 8540,
					"PBRF", 60.13,
					"KJ", 1.2,
					"KXYV", LIST(0.2, 0.0000, 0.00),
					"RBRIGZ", -453398.8776200049,
					"JBRFG", LIST(-6.09621070e-03, 0, 0),
					"SBRFG", LIST(-3.60810506e-05, 0, -1.16655679e-05),
					"RTPI", -BODY:POSITION + LS:POSITION,
					"RPI", -BODY:POSITION,
					"VPI", ORBIT:VELOCITY:ORBIT
					
).

RUN ldg_math.
RUN ldg_cser.
RUN ldg_func.

STEERINGMANAGER:RESETTODEFAULT().
SET STEERINGMANAGER:MAXSTOPPINGTIME TO 20.
SET STEERINGMANAGER:ROLLTS TO 5.
SET STEERINGMANAGER:PITCHTS TO 5.
SET STEERINGMANAGER:YAWTS TO 5.

LOCAL PACK TO APTR().
GLOBAL APTG IS PACK[0].
GLOBAL APIG IS PACK[1].

LOCAL PACK TO BRTR().
GLOBAL BRTG IS PACK[0].
GLOBAL BRIG IS PACK[1].

LOCAL T IS P63I().
LOCAL ignitiontime IS clocktime - inputs["TTRIM"].

LOCK STEERING TO RETROGRADE.

UNTIL TIME:SECONDS > ignitiontime - 85 {
	CLEARSCREEN.
	PRINT "T-" + ROUND(ignitiontime - TIME:SECONDS).
	WAIT 0.1.
}
KUNIVERSE:timewarp:cancelwarp.

LOCAL TG IS BRTG.
SET k TO 1.
SET TDRAW TO 0.
GLOBAL THROTTLECONTROL IS 1.
CLEARSCREEN.

UNTIL TIME:SECONDS > ignitiontime - 10 {
	
	LOCAL RTP IS -BODY:POSITION + LS:POSITION.
	LOCAL RP IS -BODY:POSITION.
	LOCAL VP IS ORBIT:VELOCITY:ORBIT.
	
	LOCAL deltaclocktime IS TIME:SECONDS - clocktimeOLD.
	SET clocktimeOLD TO TIME:SECONDS.

	LOCAL PACK IS GUIDE(T, TG, RTP, RP, VP, k, deltaclocktime).
	SET FCP TO PACK[0].
	SET UNFCP TO FCP:NORMALIZED.
	SET T TO PACK[3].
	
	LOCK STEERING TO UNFCP.
	
	IF TDRAW + deltaclocktime > 1 {
		CLEARSCREEN.
		PRINT "T-" + ROUND(ignitiontime - TIME:SECONDS).
		SET TDRAW TO 0.
	} ELSE {
		SET TDRAW TO TDRAW + deltaclocktime.
	}
	
	WAIT 0.1.
}

UNTIL TIME:SECONDS > ignitiontime {
	
	SET SHIP:CONTROL:FORE TO 1.0.
	LOCAL deltaclocktime IS TIME:SECONDS - clocktimeOLD.
	SET clocktimeOLD TO TIME:SECONDS.
	
	IF TDRAW + deltaclocktime > 1 {
		CLEARSCREEN.
		PRINT "T-" + ROUND(ignitiontime - TIME:SECONDS).
		SET TDRAW TO 0.
	} ELSE {
		SET TDRAW TO TDRAW + deltaclocktime.
	}
	
	WAIT 0.1.
}



UNTIL TIME:SECONDS > ignitiontime + inputs["TTRIM"] {

	IF TIME:SECONDS > ignitiontime + 2 { SET SHIP:CONTROL:FORE TO 0. }
	
	LOCAL RTP IS -BODY:POSITION + LS:POSITION.
	LOCAL RP IS -BODY:POSITION.
	LOCAL VP IS ORBIT:VELOCITY:ORBIT.
	
	LOCAL deltaclocktime IS TIME:SECONDS - clocktimeOLD.
	SET clocktimeOLD TO TIME:SECONDS.
	
	LOCAL PACK IS GUIDE(T, TG, RTP, RP, VP, k, deltaclocktime).
	SET FCP TO PACK[0].
	SET UNFCP TO FCP:NORMALIZED.
	SET T TO PACK[3].
	
	LOCK STEERING TO UNFCP.
	LOCK THROTTLE TO THROTT(inputs["AFTRIM"]).
	
	IF TDRAW + deltaclocktime > 1 {
		CLEARSCREEN.
		PRINT "CURRENT THROTTLE: " + ROUND(THROTT(inputs["AFTRIM"]), 3).
		PRINT "T: " + ROUND(T).
		SET TDRAW TO 0.
	} ELSE {
		SET TDRAW TO TDRAW + deltaclocktime.
	}
	
	WAIT 0.001.
}
SET THROTTLECONTROL TO 0.
UNTIL T > inputs["TBRF"] {
	LOCAL RTP IS -BODY:POSITION + LS:POSITION.
	LOCAL RP IS -BODY:POSITION.
	LOCAL VP IS ORBIT:VELOCITY:ORBIT.
	
	LOCAL deltaclocktime IS TIME:SECONDS - clocktimeOLD.
	SET clocktimeOLD TO TIME:SECONDS.
	
	LOCAL PACK IS GUIDE(T, TG, RTP, RP, VP, k, deltaclocktime).
	SET FCP TO PACK[0].
	LOCAL RG IS PACK[1].
	LOCAL VG IS PACK[2].
	SET UNFCP TO FCP:NORMALIZED.
	SET T TO PACK[3].
	
	LOCK THROTTLE TO THROTT(FCP:MAG).

	IF FCP:MAG < inputs["FMID"] { SET THROTTLECONTROL TO 1. }
	
	IF TDRAW + deltaclocktime > 1 {
		CLEARSCREEN.
		PRINT "CURRENT THROTTLE: " + ROUND(THROTT(FCP:MAG), 3).
		PRINT "CURRENT THRUST: " + ROUND(FCP:MAG, 3).
		PRINT "T: " + ROUND(T).
		PRINT RG.
		PRINT VG.
		SET TDRAW TO 0.
	} ELSE {
		SET TDRAW TO TDRAW + deltaclocktime.
	}
	
	WAIT 0.001.
}

SET T TO inputs["TAPI"].
SET TG TO APTG.
SET k TO 0.

UNTIL T > inputs["TAPF"] {
	LOCAL RTP IS -BODY:POSITION + LS:POSITION.
	LOCAL RP IS -BODY:POSITION.
	LOCAL VP IS ORBIT:VELOCITY:ORBIT.
	
	LOCAL deltaclocktime IS TIME:SECONDS - clocktimeOLD.
	SET clocktimeOLD TO TIME:SECONDS.
	
	LOCAL PACK IS GUIDE(T, TG, RTP, RP, VP, k, deltaclocktime).
	SET FCP TO PACK[0].
	LOCAL RG IS PACK[1].
	LOCAL VG IS PACK[2].
	SET UNFCP TO FCP:NORMALIZED.
	SET T TO PACK[3].
	
	LOCK STEERING TO UNFCP.
	LOCK THROTTLE TO THROTT(FCP:MAG).
	
	IF TDRAW + deltaclocktime > 1 {
		CLEARSCREEN.
		PRINT "CURRENT THROTTLE: " + ROUND(THROTT(FCP:MAG), 3).
		PRINT "CURRENT THRUST: " + ROUND(FCP:MAG, 3).
		PRINT "T: " + ROUND(T).
		PRINT RG.
		PRINT VG.
		SET TDRAW TO 0.
	} ELSE {
		SET TDRAW TO TDRAW + deltaclocktime.
	}
	
	WAIT 0.001.
}

TGUIDE().






