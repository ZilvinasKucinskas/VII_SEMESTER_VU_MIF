/*

2.4. tsuma(S,K):
K - skaièiø sàraðo S  teigiamø elementø suma.
	Pvz. goal:t_suma([5,-1,3],K).
	K = 8.

*/

/* ---- Þilvinas Kuèinskas 4 kursas 1 grupë PS ----*/

t_suma(S, K) :- suma(S, 0, K).

suma([], R, K) :- K is R.
suma([Y|T], R, K) :- (Y > 0 -> TEMP is R+Y; 
							 TEMP is R), 
						suma(T, TEMP, K).

/*

    Test atvejai
	t_suma([3,1,-2,4,5], K). -> Expected : K = 13
	t_suma([3,1,2,4,5], K). -> Expected : K = 15
	t_suma([], K). -> Expected : K = 0
	
*/