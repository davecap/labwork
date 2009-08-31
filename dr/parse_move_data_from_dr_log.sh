#!/bin/sh


grep -A 5 "will be incremented" t2.log | grep -A 1 -B 4 "Move rejected" > rejected_moves
echo "system change,cancellation_change,DRPE change,total change" > rejected_move_data
cat rejected_moves  | grep "cancellation_change" | cut -d ":" -f 2 | sed -e 's/^ //g' | sed -e 's/ /,/g' > rejected_move_data

grep -A 5 "will be incremented" t2.log | grep -A 1 -B 4 "Move accepted" > accepted_moves 
echo "system change,cancellation_change,DRPE change,total change" > accepted_move_data
cat accepted_moves  | grep "cancellation_change" | cut -d ":" -f 2 | sed -e 's/^ //g' | sed -e 's/ /,/g' > accepted_move_data

