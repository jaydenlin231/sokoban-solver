import itertools
import time
import search 
import sokoban
import mySokobanSolver

if __name__ == "__main__":
    # warehouse = sokoban.Warehouse()
    # warehouse.load_warehouse("./warehouses/warehouse_07.txt")
    # print(mySokobanSolver.taboo_cells(warehouse))

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
        print(mySokobanSolver.taboo_cells(warehouse))
        
