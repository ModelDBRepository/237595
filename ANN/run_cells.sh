EPOCHS=300
SEED=1

for SEED in 1; do
	(
	for CELL in 1 2 3 4 5 6 7 8; do
		python inn_cell.py  --epochs ${EPOCHS} --act linear --seed $SEED  --cell $CELL
		python inn_cell.py  --epochs ${EPOCHS} --act sub --seed $SEED  --cell    $CELL
		python inn_cell.py  --epochs ${EPOCHS} --act supra --seed $SEED   --cell $CELL
		python inn_cell.py  --epochs ${EPOCHS} --act mixed  --seed $SEED  --cell $CELL
	done
	) & 
done

