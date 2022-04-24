
'''

    Sokoban assignment


The functions and classes defined in this module will be called by a marker script. 
You should complete the functions and classes according to their specified interfaces.

No partial marks will be awarded for functions that do not meet the specifications
of the interfaces.

You are NOT allowed to change the defined interfaces.
In other words, you must fully adhere to the specifications of the 
functions, their arguments and returned values.
Changing the interfacce of a function will likely result in a fail 
for the test of your code. This is not negotiable! 

You have to make sure that your code works with the files provided 
(search.py and sokoban.py) as your code will be tested 
with the original copies of these files. 

Last modified by 2022-03-27  by f.maire@qut.edu.au
- clarifiy some comments, rename some functions
  (and hopefully didn't introduce any bug!)

'''

# You have to make sure that your code works with 
# the files provided (search.py and sokoban.py) as your code will be tested 
# with these files
import itertools
import search 
import sokoban


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def my_team():
    '''
    Return the list of the team members of this assignment submission as a list
    of triplet of the form (student_number, first_name, last_name)
    
    '''
#    return [ (1234567, 'Ada', 'Lovelace'), (1234568, 'Grace', 'Hopper'), (1234569, 'Eva', 'Tardos') ]
    raise NotImplementedError()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

class Floor(search.Problem):
    def __init__(self, warehouse):
        self.warehouse = sokoban.Warehouse.copy(warehouse)
        self.initial = warehouse.worker
        self.goal = None

    def actions(self, state):
        """
        Return the list of actions that can be executed in the given state.
        
        """
        # list of legal actions
        L = []
        # get current state worker coords
        [worker_x, worker_y] = state 
        # check if action results in the worker on top of a wall, if not, append to legal action
        # target and box coords are ignored for this problem of finding inside wall area.
        if(worker_x + 1, worker_y) not in self.warehouse.walls:
            L.append("Right")
        if(worker_x - 1, worker_y) not in self.warehouse.walls:
            L.append("Left")
        if(worker_x, worker_y + 1) not in self.warehouse.walls:
            L.append("Down")
        if(worker_x, worker_y - 1) not in self.warehouse.walls:
            L.append("Up")
        
        return L

    def result(self, state, action):
        """
        Return the state that results from executing the given
        action in the given state. The action must be one of
        self.actions(state).
        """

        next_worker_coords = list(state)
        
        assert action in self.actions(state)  # defensive programming!
        if action == 'Up':
            next_worker_coords[1] = next_worker_coords[1] - 1
        if action == 'Down':
            next_worker_coords[1] = next_worker_coords[1] + 1
        if action == 'Left':
            next_worker_coords[0] = next_worker_coords[0] - 1
        if action == 'Right':
            next_worker_coords[0] = next_worker_coords[0] + 1

        # self.warehouse = sokoban.Warehouse.copy(self.warehouse, worker=next_worker_coords)

        return tuple(next_worker_coords)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def floor_breadth_first_graph_search(problem):
    "Graph search version of BFS.  [Fig. 3.11]"
    return floor_graph_search(problem, search.FIFOQueue())

def floor_graph_search(problem, frontier):
    """
    Search through the successors of a problem to find a goal.
    The argument frontier should be an empty queue.
    If two paths reach a state, only use the first one. [Fig. 3.7]
    Return
        the node of the first goal state found
        or None is no goal state is found
    """
    assert isinstance(problem, search.Problem)
    frontier.append(search.Node(problem.initial))
    explored = set() # initial empty set of explored states
    while frontier:
        node = frontier.pop()
        if problem.goal_test(node.state):
            return node
        explored.add(node.state)
        # Python note: next line uses of a generator
        frontier.extend(child for child in node.expand(problem)
                        if child.state not in explored
                        and child not in frontier)
    return explored

def get_floor_cells(floor):
    return floor_breadth_first_graph_search(floor)
    
def has_up_neighbour(cell, allCells):
    return (cell[0], cell[1] - 1) in allCells

def has_down_neighbour(cell, allCells):
    return (cell[0], cell[1] + 1) in allCells

def has_right_neighbour(cell, allCells):
    return (cell[0] + 1, cell[1]) in allCells

def has_left_neighbour(cell, allCells):
    return (cell[0] - 1, cell[1]) in allCells

def has_any_neighbour(cell, allCells):
    return has_up_neighbour(cell, allCells) \
           or has_down_neighbour(cell, allCells)\
           or has_right_neighbour(cell, allCells)\
           or has_left_neighbour(cell, allCells)

def is_bot_left_corner(cell, allCells):
    return has_up_neighbour(cell, allCells) \
           and has_right_neighbour(cell, allCells) \
           and not has_down_neighbour(cell, allCells) \
           and not has_left_neighbour(cell, allCells) 
           
def is_bot_right_corner(cell, allCells):
    return has_up_neighbour(cell, allCells) \
           and has_left_neighbour(cell, allCells) \
           and not has_down_neighbour(cell, allCells) \
           and not has_right_neighbour(cell, allCells) 
           
def is_top_left_corner(cell, allCells):
    return has_right_neighbour(cell, allCells) \
           and has_down_neighbour(cell, allCells) \
           and not has_left_neighbour(cell, allCells) \
           and not has_up_neighbour(cell, allCells)

def is_top_right_corner(cell, allCells):
    return has_left_neighbour(cell, allCells) \
           and has_down_neighbour(cell, allCells) \
           and not has_right_neighbour(cell, allCells) \
           and not has_up_neighbour(cell, allCells)

def is_corner(cell, allCells):
    return is_bot_left_corner(cell, allCells) \
           or is_bot_right_corner(cell, allCells) \
           or is_top_right_corner(cell, allCells) \
           or is_top_left_corner(cell, allCells) 

# def get_floor_area_corners(floor_cells):
#     return {cell for cell in floor_cells if is_corner(cell, floor_cells)}

get_floor_area_corners = lambda floor_cells : {cell for cell in floor_cells if is_corner(cell, floor_cells)}

are_inline_vert = lambda cellA, cellB : cellA[0] == cellB[0]
are_inline_hor = lambda cellA, cellB : cellA[1] == cellB[1]

are_in_line = lambda cellA, cellB : are_inline_hor(cellA, cellB) or are_inline_vert(cellA, cellB)

def space_in_line(cellA, cellB):
    if are_inline_hor(cellA, cellB):
        y = cellA[1]
        x_start = min(cellA[0], cellB[0]) + 1
        x_end = max(cellA[0], cellB[0])
        if(x_start < x_end and x_start > 0):
            return {(x, y) for x in range(x_start, x_end)}
    elif are_inline_vert(cellA, cellB):
        x = cellA[0]
        y_start = min(cellA[1], cellB[1]) + 1
        y_end = max(cellA[1], cellB[1])
        if(y_start < y_end and y_start > 0):
            return {(x, y) for y in range(y_start, y_end)}


def taboo_cells(warehouse):
    '''  
    Identify the taboo cells of a warehouse. A "taboo cell" is by definition
    a cell inside a warehouse such that whenever a box get pushed on such 
    a cell then the puzzle becomes unsolvable. 
    
    Cells outside the warehouse are not taboo. It is a fail to tag an 
    outside cell as taboo.
    
    When determining the taboo cells, you must ignore all the existing boxes, 
    only consider the walls and the target  cells.  
    Use only the following rules to determine the taboo cells;
     Rule 1: if a cell is a corner and not a target, then it is a taboo cell.
     Rule 2: all the cells between two corners along a wall are taboo if none of 
             these cells is a target.
    
    @param warehouse: 
        a Warehouse object with the worker inside the warehouse

    @return
       A string representing the warehouse with only the wall cells marked with 
       a '#' and the taboo cells marked with a 'X'.  
       The returned string should NOT have marks for the worker, the targets,
       and the boxes.  
    '''
    floor_area = get_floor_cells(Floor(warehouse))
    rule_1 = {cell for cell in get_floor_area_corners(floor_area) if cell not in warehouse.targets}
    # Prepare results
    X,Y = zip(*warehouse.walls) 
    x_size, y_size = 1+max(X), 1+max(Y)
    
    vis = [[" "] * x_size for y in range(y_size)]

    
    rule_2 = set()
    for combination in itertools.combinations(rule_1, 2):
        cellA = combination[0]
        cellB = combination[1]
        if not are_in_line(cellA, cellB):
            continue

        if(not space_in_line(cellA, cellB)):
            continue

        print("---------------------")
        print(f"{cellA} and {cellB}")
        print("space in line")
        print(space_in_line(cellA, cellB))
        not_taboo_blocks = False
        for space in space_in_line(cellA, cellB):
            if(not(space in floor_area and has_any_neighbour(space, warehouse.walls))):
                not_taboo_blocks = True
                print("not taboo")
                break   

        if not_taboo_blocks:
            continue

        print("taboo")
        for space in space_in_line(cellA, cellB):
            rule_2.add(space)    
            
    for (x,y) in warehouse.walls:
        vis[y][x] = "#"
    for (x,y) in rule_1.union(rule_2):
        vis[y][x] = "X"

    # print(get_floor_area_corners(floor_area))
    return "\n".join(["".join(line) for line in vis])

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


class SokobanPuzzle(search.Problem):
    '''
    An instance of the class 'SokobanPuzzle' represents a Sokoban puzzle.
    An instance contains information about the walls, the targets, the boxes
    and the worker.

    Your implementation should be fully compatible with the search functions of 
    the provided module 'search.py'. 
    
    '''
    
    #
    #         "INSERT YOUR CODE HERE"
    #
    #     Revisit the sliding puzzle and the pancake puzzle for inspiration!
    #
    #     Note that you will need to add several functions to 
    #     complete this class. For example, a 'result' method is needed
    #     to satisfy the interface of 'search.Problem'.
    #
    #     You are allowed (and encouraged) to use auxiliary functions and classes

    
    def __init__(self, warehouse):
        raise NotImplementedError()

    def actions(self, state):
        """
        Return the list of actions that can be executed in the given state.
        
        """
        raise NotImplementedError

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def check_elem_action_seq(warehouse, action_seq):
    '''
    
    Determine if the sequence of actions listed in 'action_seq' is legal or not.
    
    Important notes:
      - a legal sequence of actions does not necessarily solve the puzzle.
      - an action is legal even if it pushes a box onto a taboo cell.
        
    @param warehouse: a valid Warehouse object

    @param action_seq: a sequence of legal actions.
           For example, ['Left', 'Down', Down','Right', 'Up', 'Down']
           
    @return
        The string 'Impossible', if one of the action was not valid.
           For example, if the agent tries to push two boxes at the same time,
                        or push a box into a wall.
        Otherwise, if all actions were successful, return                 
               A string representing the state of the puzzle after applying
               the sequence of actions.  This must be the same string as the
               string returned by the method  Warehouse.__str__()
    '''
    
    ##         "INSERT YOUR CODE HERE"
    
    raise NotImplementedError()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def solve_weighted_sokoban(warehouse):
    '''
    This function analyses the given warehouse.
    It returns the two items. The first item is an action sequence solution. 
    The second item is the total cost of this action sequence.
    
    @param 
     warehouse: a valid Warehouse object

    @return
    
        If puzzle cannot be solved 
            return 'Impossible', None
        
        If a solution was found, 
            return S, C 
            where S is a list of actions that solves
            the given puzzle coded with 'Left', 'Right', 'Up', 'Down'
            For example, ['Left', 'Down', Down','Right', 'Up', 'Down']
            If the puzzle is already in a goal state, simply return []
            C is the total cost of the action sequence C

    '''
    
    raise NotImplementedError()

if __name__ == "__main__":
    wh = sokoban.Warehouse()
    wh.load_warehouse("./warehouses/warehouse_03.txt")
    
    print(wh)
    # print("\n")
    print(taboo_cells(wh))
    # for i in space_in_line((1, 2), (1, 10)):
    #     if(i):
    #         print(type(i))
    #         print(i)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

