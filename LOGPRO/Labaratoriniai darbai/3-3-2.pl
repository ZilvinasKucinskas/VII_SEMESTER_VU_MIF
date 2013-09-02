/*


3.2. iterpti(S,K,R):
sàraðas R gautas ið duotojo skaièiø sàraðo S, áterpus á S duotajá skaièiø K,
kad K 'kaimynas' ið kairës sàraðe R bûtø maþesnis, o ið deðinës - didesnis
negu K.
	Pvz. goal: iterpti([10,2,14,8,1], 13, R).
	R = [10,2,13,14,8,1].

	
*/

/* ---- Þilvinas Kuèinskas 4 kursas 1 grupë PS ----*/

iterpti(S, K, R) :- iterpti2(S, K, [], R).

iterpti2([X,Y|T], K, Temp, R) :- (X < K -> 
							(Y > K ->
								prideti_prie_saraso(X, Temp, New),
								prideti_prie_saraso(K, New, New2),
								iterpti_iki_galo([Y|T], New2, R)
								;
								prideti_prie_saraso(X, Temp, New),
								iterpti2([Y|T], K, New, R)
							)
							;
							prideti_prie_saraso(X, Temp, New),
							iterpti2([Y|T], K, New, R)
						   ).

						   /*
it() :- X < K, Y < K, !, prideti_prie_saraso(X, Temp, New),
								prideti_prie_saraso(K, New, New2),
								iterpti_iki_galo([Y|T], New2, R).
it() :- 							prideti_prie_saraso(X, Temp, New),
							iterpti2([Y|T], K, New, R).
						

test(petras).
test(jonas) :- fail, !.
test(antanas).				

*/

		
/* Grazina saraso dydi */
size([],0).
size([H|T],N) :- size(T,N1), N is N1+1.
						   
/* 
	prideti_prie_saraso(X,L,L1).
	Prideda i saraso gala elementa
	ir grazina L1
*/			
   
prideti_prie_saraso(X,[],[X]).
prideti_prie_saraso(X,[A|L],[A|L1]):- prideti_prie_saraso(X,L,L1).

/*
	saraso elementus iterpia prie saraso List ir issaugo reza Final
*/
iterpti_iki_galo([], List, List).
iterpti_iki_galo([Y|T], List, Final) :- prideti_prie_saraso(Y, List, New), iterpti_iki_galo(T, New, Final).