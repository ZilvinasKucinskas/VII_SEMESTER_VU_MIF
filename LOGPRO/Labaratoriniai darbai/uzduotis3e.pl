/*

3. Faktais uþraðyta logikos algebros konjunkcijos bei neiginio reikðmiø
lentelë. Apibrëþti predikatus
       e) "X,Y priklauso vienai tiesei A*X + B = Y.
	   
*/

/* ---- Þilvinas Kuèinskas 4 kursas 1 grupë PS ----*/

neiginys('true', 'false').
neiginys('false', 'true').

konjunkcija('true', 'true', 'true').
konjunkcija('true', 'false', 'false').
konjunkcija('false', 'true', 'false').
konjunkcija('false', 'false', 'false').

/*
     A - konstanta
	 B - konstanta
	 X - kintamieji
	 Y - kintamieji
	 
	 logine daugyba - konjunkcija - and
	 logine sudetis - disjunkcija - or
		
		a and b
		true true -> true
		false true -> false
		true false -> false
		false false -> false
		
		a or b
		true true -> true
		false true -> true
		true false -> true
		false false -> false
		
		ne(ne(a) ir ne(b)) -> logine sudetis per neigima ir konjunkcija
		true true -> true
		true false -> true
		false true -> true
		false false -> false
*/

priklauso_tiesei(A, X, B, Y) :- konjunkcija(A, X, Z), 
								neiginys(Z, NeZ), 
								neiginys(B, NeB),
								konjunkcija(NeB, NeZ, TempKonjunkcija), 
								neiginys(TempKonjunkcija, RezultatasY), 
								RezultatasY=Y.
	
/*
	1 variantas:
	priklauso_tiesei('true', 'true', 'true', 'true').

	true and true -> true
	true or true -> true
	true == true.
	
	Expected value = yes.
*/
	
	
/*
	2 variantas:
	priklauso_tiesei('false', 'true', 'false', 'false').

	false and true -> false
	false or false -> false
	false == false.
	
	Expected value = yes.
*/

/*
	3 variantas:
	priklauso_tiesei('false', 'true', 'false', 'true').

	false and true -> false
	false or false -> false
	false == true.
	
	Expected value = no.
*/
