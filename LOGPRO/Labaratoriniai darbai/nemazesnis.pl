/*

Tikrinam ar pirmasis skaicius yra nemazesnis uz pirmaji :)

*/

/* ---- Þilvinas Kuèinskas 4 kursas 1 grupë PS ----*/

nemazesnis(_, nul).
nemazesnis(s(X), s(Y)) :- nemazesnis(X, Y).

/*

	nemazesnis(nul, nul). -> yes Expected : yes
	
	nemazesnis(s(nul), nul). -> yes Expected : yes
	
	nemazesnis(s(nul), s(nul)). -> yes Expected : yes
	
	nemazesnis(s(s(nul)), s(nul)). -> yes Expected : yes 
	
	nemazesnis(s(nul), s(s(nul))). -> no Expected : no
	
	5 TESTS CORRECT
	0 TESTS FAIL
*/