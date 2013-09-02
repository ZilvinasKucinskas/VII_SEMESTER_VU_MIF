/*

4.5. skaic_skaitm(Skaic,Sar):
Sar - duotojo deðimtainio skaièiaus Skaic skaitmenø sàraðas.
	Pvz. goal: skaic_skaitm(1456,Sar).
	Sar = [1,4,5,6].

*/

/* ---- Þilvinas Kuèinskas 4 kursas 1 grupë PS ----*/

skaic_skaitm(Skaic, Sar) :- paversti(Skaic, [], Sar).

paversti(0, TSar, Sar) :- Sar = TSar.
paversti(Skaic, TSar, Sar) :- (Skaic < 10 -> TEMP is Skaic,
									prideti_prie_pradzios(TEMP, TSar, NAUJAS),
									paversti(0, NAUJAS, Sar);
									
									LIEKANAMOD is Skaic mod 10,
									LIEKANADIV is Skaic div 10,
									prideti_prie_pradzios(LIEKANAMOD, TSar, NAUJAS),
									paversti(LIEKANADIV, NAUJAS, Sar)
						).

						
/* 
	prideti_prie_pradzios(X,L,L1).
	Prideda i saraso pradzia elementa
	ir grazina L1
*/	
prideti_prie_pradzios(X, L, [X|L]).



/*

	TEST CASES:
	
	1.
	skaic_skaitm(1456,Sar).
	Expected : Sar = [1,4,5,6].
	2.
	skaic_skaitm(1,Sar).
	Expected : Sar = [1].
	3.
	skaic_skaitm(12,Sar).
	Expected : Sar = [1,2].
	4.
	skaic_skaitm(54321,Sar).
	Expected : Sar = [5,4,3,2,1].

*/