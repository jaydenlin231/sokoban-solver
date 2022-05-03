import itertools
import time
import search 
import sokoban
import mySokobanSolver

if __name__ == "__main__":
    test = ["./warehouses/warehouse_07.txt", 
            "./warehouses/warehouse_09.txt",
            "./warehouses/warehouse_47.txt",
            "./warehouses/warehouse_81.txt",
            "./warehouses/warehouse_147.txt",
            "./warehouses/warehouse_5n.txt"]

    for i in range(len(test)):
        warehouse = sokoban.Warehouse()
        warehouse.load_warehouse(test[i])
        print("-------------------------------------")
        print(test[i])
        t0 = time.time()
        solution, total_cost = mySokobanSolver.solve_weighted_sokoban(warehouse)
        t1 = time.time()
        print (f'\nAnalysis took {t1-t0:.6f} seconds\n')
        if solution == 'impossible':
            print('\nNo solution found!\n')
        else:
            print(f"\nSolution found with a cost of {total_cost} \n", solution, '\n')

        
        
