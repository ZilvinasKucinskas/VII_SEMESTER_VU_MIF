/*

1.2. maz(S):
	duoto sveikøjø skaièiø sàraðo S elementai iðdëstyti maþejimo tvarka.
	Pvz. goal: maz([41,8,2]).
	true.
	   
*/

/* ---- Þilvinas Kuèinskas 4 kursas 1 grupë PS ----*/

maz([]). /* jeigu sarasas tuscias tai cia true viskas tvarkoj */
maz([_]). /* jeigu sarase vienas elementas tai cia true viskas tvarkoj */
maz([X,Y|T]):- X >= Y, 
			   maz([Y|T]).


/* 
	Testavimo atvejai
	
	maz([]). -> Expected: true
	maz([2]). -> Expected : true
	maz([4,3,2]). -> Expected: true
	maz([5,3,7]). -> Expected : false
	maz([5,5,5]). -> Expected : true


*/