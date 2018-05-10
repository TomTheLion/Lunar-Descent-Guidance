//Matrix inversion, inverts a 3x3 matrix
FUNCTION MINV {
	PARAMETER m1.
	
	LOCAL A IS  (m1[1][1] * m1[2][2] - m1[1][2] * m1[2][1]).
	LOCAL B IS -(m1[1][0] * m1[2][2] - m1[1][2] * m1[2][0]).
	LOCAL C IS  (m1[1][0] * m1[2][1] - m1[1][1] * m1[2][0]).
	
	LOCAL D IS -(m1[0][1] * m1[2][2] - m1[0][2] * m1[2][1]).
	LOCAL E IS  (m1[0][0] * m1[2][2] - m1[0][2] * m1[2][0]).
	LOCAL F IS -(m1[0][0] * m1[2][1] - m1[0][1] * m1[2][0]).
	
	LOCAL G IS  (m1[0][1] * m1[1][2] - m1[0][2] * m1[1][1]).
	LOCAL H IS -(m1[0][0] * m1[1][2] - m1[0][2] * m1[1][0]).
	LOCAL I IS  (m1[0][0] * m1[1][1] - m1[0][1] * m1[1][0]).
	
	LOCAL det IS m1[0][0] * A + m1[0][1] * B + m1[0][2] * C.
		
	LOCAL row1 IS LIST(A / det, D / det, G / det).
	LOCAL row2 IS LIST(B / det, E / det, H / det).
	LOCAL row3 IS LIST(C / det, F / det, I / det).
	
	RETURN LIST(row1, row2, row3).
}

//Matrix transpose, tranposes a 2d matrix
FUNCTION MTRNS {
	PARAMETER m1.
	
	LOCAL n IS m1:LENGTH.
	LOCAL xlist IS LIST().
	LOCAL ylist IS LIST().
	
	FROM { LOCAL i IS 0. } UNTIL i = n STEP { SET i TO i + 1. } DO {				
		FROM { LOCAL j IS 0. } UNTIL j = n STEP { SET j TO j + 1. } DO {
			xlist:ADD(m1[j][i]).
		}
		ylist:ADD(xlist).
		SET xlist TO LIST().
	}
	
	RETURN ylist.
}

//Matrix magnitude
FUNCTION MMAG {
	PARAMETER m1.
	
	RETURN (m1[0] ^ 2 + m1[1] ^ 2 + m1[1] ^ 2) ^ 0.5.
}

//Matrix multiplication, multiplies a 2d matrix by a 1d or 2d matrix
FUNCTION MMUL {
	PARAMETER m1.
	PARAMETER m2.
	
	LOCAL n IS m1:LENGTH.
	LOCAL m IS m1[0]:LENGTH.
	LOCAL xval IS 0.
	LOCAL ylist IS LIST().	
	
	IF m2[0]:TYPENAME() = "List" {
		LOCAL p IS m2[0]:LENGTH.
		LOCAL xlist IS LIST().

		FROM { LOCAL i IS 0. } UNTIL i = n STEP { SET i TO i + 1. } DO {
			FROM { LOCAL j IS 0. } UNTIL j = p STEP { SET j TO j + 1. } DO {				
				FROM { LOCAL k IS 0. } UNTIL k = m STEP { SET k TO k + 1. } DO {
					SET xval TO xval + m1[i][k] * m2[k][j].
				}
				xlist:ADD(xval).
				SET xval TO 0.
			}
			ylist:ADD(xlist).
			SET xlist TO LIST().
		}
	} ELSE {	
		FROM { LOCAL i IS 0. } UNTIL i = n STEP { SET i TO i + 1. } DO {				
			FROM { LOCAL k IS 0. } UNTIL k = m STEP { SET k TO k + 1. } DO {
				SET xval TO xval + m1[i][k] * m2[k].
			}
			ylist:ADD(xval).
			SET xval TO 0.	
		}
	}
	
	RETURN ylist.
}

//Scalar multiplication, multiplies a 1d or 2d matrix by a scalar
FUNCTION SMUL {
	PARAMETER s.
	PARAMETER m1.
	
	LOCAL n IS m1:LENGTH.
	LOCAL xval IS 0.
	LOCAL ylist IS LIST().	
	
	IF m1[0]:TYPENAME() = "List" {
		LOCAL m IS m1[0]:LENGTH.
		LOCAL xlist IS LIST().

		FROM { LOCAL i IS 0. } UNTIL i = n STEP { SET i TO i + 1. } DO {
			FROM { LOCAL j IS 0. } UNTIL j = m STEP { SET j TO j + 1. } DO {
				xlist:ADD(s * m1[i][j]).
			}
			ylist:ADD(xlist).
			SET xlist TO LIST().
		}
	} ELSE {	
		FROM { LOCAL i IS 0. } UNTIL i = n STEP { SET i TO i + 1. } DO {				
			ylist:ADD(s * m1[i]).
		}
	}
	
	RETURN ylist.
}

//Matrix addition, adds two 1d or 2d matrices
FUNCTION MADD {
	PARAMETER m1.
	PARAMETER m2.

	LOCAL n IS m1:LENGTH.
	LOCAL ylist IS LIST().
	
	IF m1[0]:TYPENAME() = "List" {
		LOCAL m IS m1[0]:LENGTH.
		LOCAL xlist IS LIST().

		FROM { LOCAL i IS 0. } UNTIL i = n STEP { SET i TO i + 1. } DO {
			FROM { LOCAL j IS 0. } UNTIL j = m STEP { SET j TO j + 1. } DO {
				xlist:ADD(m1[i][j] + m2[i][j]).
			}
		ylist:ADD(xlist).
		SET xlist TO LIST().
		}
	} ELSE {	
		FROM { LOCAL i IS 0. } UNTIL i = n STEP { SET i TO i + 1. } DO {
			ylist:ADD(m1[i] + m2[i]).	
		}
	}
	
	RETURN ylist.
}

//Matrix subtraction, subtracts a 1d or 2d matrix from another
FUNCTION MSUB {
	PARAMETER m1.
	PARAMETER m2.
	
	LOCAL ylist IS LIST().
	
	IF m1[0]:TYPENAME() = "List" {
		LOCAL n IS m1:LENGTH.
		LOCAL m IS m1[0]:LENGTH.
		LOCAL xlist IS LIST().
		FROM { LOCAL i IS 0. } UNTIL i = n STEP { SET i TO i + 1. } DO {
			FROM { LOCAL j IS 0. } UNTIL j = m STEP { SET j TO j + 1. } DO {
				xlist:ADD(m1[i][j] - m2[i][j]).
			}
		ylist:ADD(xlist).
		SET xlist TO LIST().
		}
	} ELSE {
		LOCAL n IS m1:LENGTH.
		
		FROM { LOCAL i IS 0. } UNTIL i = n STEP { SET i TO i + 1. } DO {
			ylist:ADD(m1[i] - m2[i]).	
		}
	}
	
	RETURN ylist.
}

//Factorial, calculates a factorial
FUNCTION FACT {
	PARAMETER x.
	
	LOCAL y IS 1.
	
	UNTIL x = 1 OR x = 0 {
		SET y TO y * x.
		SET x TO x - 1.	
	}

	RETURN y.
}

//Rodrigues vector rotation
FUNCTION rodrigues {
	PARAMETER vector.
	PARAMETER axis.
	PARAMETER angle.
	
	SET axis TO axis:NORMALIZED.
	
	RETURN vector * COS(angle) + VCRS(axis, vector) * SIN(angle) + axis * VDOT(axis, vector) * (1 - COS(angle)).
}

















