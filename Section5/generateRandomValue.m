function new_value =generateRandomValue(current_value, nlevels)
    all_possible_values = 1:nlevels;
    all_possible_values(all_possible_values == current_value) = []; 
    new_value_idx = randi(length(all_possible_values)); 
    new_value = all_possible_values(new_value_idx);
end
