/* Duom. bazëje saugomi duomenys apie asmenis ir jø giminystës ryðius faktais:
	asmuo(Vardas, Lytis, Amþius, Pomëgis);
	mama (Mama,  Vaikas).
	pora (Vyras, Þmona);

Apibrëþti ðiuos predikatus:
	5) brolis_ir_sesuo(X,Y);
	---
	33) vidurinis(X) : " vaikas, ne vyriausias ir ne jauniausias ðeimoje"; */
	
/* ---- Þilvinas Kuèinskas 4 kursas 1 grupë PS ----*/

/*
									   Aleksas-Ona
											|
											|
											|
			-----------------------------------------------------------------------
			|						        |									  |
	     Egidijus-Nijolë			    Liutauras-Cikita						Audrius-Egle
			|						        |									  |
	  		|	                      _____________________				   ________________________
	        |			  			  |		   |		  |				   |		  |			  |
	     Þilvinas				    Anton	  Petr		Ramune			Aurimas		Marius		Raminta
*/	

/* asmuo(vardas, lytis, amzius, pomegis) */

asmuo('Aleksas', 'vyras', 70, 'verslas'  ).
asmuo('Ona', 'moteris', 68, 'ukis'       ).
asmuo('Egidijus', 'vyras', 50, 'verslas' ).
asmuo('Nijole', 'moteris', 49, 'darbas'  ).
asmuo('Liutauras', 'vyras', 48, 'verslas').
asmuo('Cikita', 'moteris', 25, 'mada'    ).
asmuo('Audrius', 'vyras', 52, 'darbas'   ).
asmuo('Egle', 'moteris', 51, 'mezgimas'  ).
asmuo('Zilvinas', 'vyras', 22, 'pokeris' ).
asmuo('Anton', 'vyras', 12, 'kartingai'  ).
asmuo('Petr', 'vyras', 14, 'matematika'  ).
asmuo('Ramune', 'moteris', 17, 'menas'    ).
asmuo('Aurimas', 'vyras', 23, 'grojimas' ).
asmuo('Marius', 'vyras', 25, 'zaidimai'  ).
asmuo('Raminta', 'moteris', 22, 'verimas').

/* mama(mama, vaikas) */

mama('Ona', 'Egidijus').
mama('Ona', 'Liutauras').
mama('Ona', 'Audrius').

mama('Nijole', 'Zilvinas').

mama('Cikita', 'Anton').
mama('Cikita', 'Petr').
mama('Cikita', 'Ramune').

mama('Egle', 'Aurimas').
mama('Egle', 'Marius').
mama('Egle', 'Raminta').

/* pora(vyras, zmona) */

pora('Aleksas', 'Ona').
pora('Egidijus', 'Nijole').
pora('Liutauras', 'Cikita').
pora('Audrius', 'Egle').

/* 5) brolis_ir_sesuo(X,Y); */
brolis_ir_sesuo(X,Y) :- mama(Mama, X), mama(Mama, Y), asmuo(X, LV,_,_), asmuo(Y, LM, _, _), LV \= LM.

/* 33) vidurinis(X) : " vaikas, ne vyriausias ir ne jauniausias ðeimoje"; */
	 
vidurinis(X) :- asmuo(X,_,M1,_),
				asmuo(Y,_,M2,_), 
				asmuo(Z,_,M3,_), 
				M1<M2, 
				M1>M3, 
				mama(Mama, X), 
				mama(Mama, Y), 
				mama(Mama, Z).
/*
bendraamziai('Cikita', 'Marius').
bendraamziai('Cikita', 'Egle').
vidurinis('Anton').
vidurinis('Petr').
*/
