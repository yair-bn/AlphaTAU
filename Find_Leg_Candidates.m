nums = [1:1:70];
P_nums_idx = isprime(nums);
P_nums = nums(P_nums_idx);

leg_cond = mod(P_nums,4)==3
P_nums(leg_cond)
