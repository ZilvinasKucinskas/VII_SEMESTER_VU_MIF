/*

sqrt_of_power_of_two(X) - predicate, which returns 
		true	if exist such A for pow(A, 2) == X

*/

/* ---- TooHighToPlay ----*/

sqrt_of_power_of_two(X) :- Current_value is X - 2, sqrt_of_power_of_two_helper(X, 2, Current_value).

sqrt_of_power_of_two_helper(Number,Current_multiplicator,0).
sqrt_of_power_of_two_helper(1,Current_multiplicator, Current_value).
sqrt_of_power_of_two_helper(Number, Current_multiplicator, Current_value) :-
	New_multiplicator is Current_multiplicator * 2,
	New_value is Current_value - Current_multiplicator,
	Current_value >= 0,
	sqrt_of_power_of_two_helper(Number, New_multiplicator, New_value).
	
/*

        TEST CASES:
        
        1.
        sqrt_of_power_of_two(0).
        true.
        2.
        sqrt_of_power_of_two(16).
        true.
        3.
        sqrt_of_power_of_two(512).
		true.
        4.
        sqrt_of_power_of_two(15).
		false.
*/