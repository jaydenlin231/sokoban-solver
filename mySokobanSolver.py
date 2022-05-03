"""

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

"""

import itertools
import search
import sokoban

# Code has been autoformatted using 'black'
# (a PEP 8 compliant opinionated formatter)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Team Declaration
def my_team():
    """
    Return the list of the team members of this assignment submission as a list
    of triplet of the form (student_number, first_name, last_name)

    """
    return [(10664599, "Wei-Chung", "Lin"), (10794565, "Alexander", "Kim")]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Coordinate Utilities


def move_up(coords):
    """
    Return the resulting (x,y) coordinates after applying an up move on the
    specified input coordinates. This function does not consider any rules
    and simply decrement the y coordinate by 1. The y coordinate corresponds to
    the row number in the sokoban, hence it decreases if moving up.

    @param coords:
        A tuple representing (x,y) coordinates

    @return:
        A tuple representing (x,y) coordinates after an up move has been
        applied to the input tuple

    """
    return (coords[0], coords[1] - 1)


def move_down(coords):
    """
    Return the resulting (x,y) coordinates after applying a down move on the
    specified input coordinates. This function does not consider any rules
    and simply increment the y coordinate by 1. The y coordinate corresponds to
    the row number in the sokoban, hence it increases if moving down.

    @param coords:
        A tuple representing (x,y) coordinates

    @return:
        A tuple representing (x,y) coordinates after a down move has been
        applied to the input tuple

    """
    return (coords[0], coords[1] + 1)


def move_right(coords):
    """
    Return the resulting (x,y) coordinates after applying a right move on the
    specified input coordinates. This function does not consider any rules
    and simply increment the x coordinate by 1. The x coordinate corresponds to
    the column number in the sokoban, hence it increases if moving right.


    @param coords:
        A tuple representing (x,y) coordinates

    @return:
        A tuple representing the (x,y) coordinates after a right move has been
        applied to the input tuple

    """
    return (coords[0] + 1, coords[1])


def move_left(coords):
    """
    Return the resulting (x,y) coordinates after applying a left move on the
    specified input coordinates. This function does not consider any rules
    and simply decrement the x coordinate by 1. The x coordinate corresponds to
    the column number in the sokoban, hence it decreases if moving left.

    @param coords:
        A tuple representing (x,y) coordinates

    @return:
        A tuple representing the (x,y) coordinates after a left move has been
        applied to the input tuple

    """
    return (coords[0] - 1, coords[1])


def has_up_neighbour(coords, all_coords):
    """
    Determine if a tuple representing (x,y) coordinates has an up neighbour
    in a list of (x,y) coordinates tuples.

    @param coords:
        A tuple representing (x,y) coordinates

    @param all_coords:
        A list of (x,y) coordinates tuples to check for neighbours

    @return:
       True if the input (x,y) coordinates tuple has an up neighbour in
       the list of (x, y) coordinates tuples to check.
       False otherswise.

    """
    return move_up(coords) in all_coords


def has_down_neighbour(coords, all_coords):
    """
    Determine if a tuple representing (x,y) coordinates has a down neighbour
    in a list of (x,y) coordinates tuples.

    @param coords:
        A tuple representing (x,y) coordinates

    @param all_coords:
        A list of (x,y) coordinates tuples to check for neighbours

    @return:
       True if the input (x,y) coordinates tuple has a down neighbour in
       the list of (x, y) coordinates tuples to check.
       False otherswise.

    """
    return move_down(coords) in all_coords


def has_right_neighbour(coords, all_coords):
    """
    Determine if a tuple representing (x,y) coordinates has a right neighbour
    in a list of (x,y) coordinates tuples.

    @param coords:
        A tuple representing (x,y) coordinates

    @param all_coords:
        A list of (x,y) coordinates tuples to check for neighbours

    @return:
       True if the input (x,y) coordinates tuple has a right neighbour in
       the list of (x, y) coordinates tuples to check.
       False otherswise.

    """
    return move_right(coords) in all_coords


def has_left_neighbour(coords, all_coords):
    """
    Determine if a tuple representing (x,y) coordinates has a left neighbour
    in a list of (x,y) coordinates tuples.

    @param coords:
        A tuple representing (x,y) coordinates

    @param all_coords:
        A list of (x,y) coordinates tuples to check for neighbours

    @return:
       True if the input (x,y) coordinates tuple has a left neighbour in
       the list of (x, y) coordinates tuples to check.
       False otherswise.

    """
    return move_left(coords) in all_coords


def has_any_neighbour(coords, all_coords):
    """
    Determine if a tuple representing (x,y) coordinates has any neighbour
    in a list of (x,y) coordinates tuples.

    @param coords:
        A tuple representing (x,y) coordinates

    @param all_coords:
        A list of (x,y) coordinates tuples to check for neighbours

    @return:
       True if the input (x,y) coordinates tuple has any neighbour in
       the list of (x, y) coordinates tuples to check.
       False otherswise.

    """
    return (
        has_up_neighbour(coords, all_coords)
        or has_down_neighbour(coords, all_coords)
        or has_right_neighbour(coords, all_coords)
        or has_left_neighbour(coords, all_coords)
    )


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Distance Utilities


def get_distance(point_a, point_b):
    """
    Calculate the Absolute Distance between two 1D points.

    @param point_a:
        A integer representing a 1D point to calculate the Absolute Distance

    @param point_b:
        A integer representing another 1D point to calculate the Absolute Distance

    @return:
        An integer representing the Absolute Distance between point_a and
        point_b

    """
    return abs(point_a - point_b)


def get_manhattan_distance(coords_a, coords_b):
    """
    Calculate the Manhattan Distance between a (x,y) cell coordinates
    tuple to another (x,y) cell coordinates tuple.

    @param coords_a:
        A tuple representing (x,y) coordinates to calculate the Manhanttan Distance

    @param coords_b:
        Another tuple representing (x,y) coordinates to calculate the Manhanttan Distance

    @return:
        An integer representing the Manhanttan Distance between coords_a and
        coords_b

    """
    return get_distance(coords_a[0], coords_b[0]) + get_distance(
        coords_a[1], coords_b[1]
    )


def get_min_manhattan_distance(coords, all_coords):
    """
    Calculate the minimum Manhattan Distance between a (x,y) cell coordinates
    tuple to any cell in the list of (x,y) coordinates tuples.

    @param coords:
       A tuple representing (x,y) coordinates which is the start of the Manhanttan Distance

    @param all_coords:
        A list of (x,y) coordinates tuples to be the ends of the Manhanttan Distance

    @return:
       An integer representing the minimun Manhanttan Distance between the input coordinates
       tuple to another coordinates tuple in the list.

    """
    return min(
        [
            get_manhattan_distance(coords, another_coords)
            for another_coords in all_coords
        ]
    )


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Taboo Cell Utilities


def is_bot_left_corner(floor_coords, all_floor_coords):
    """
    Determine if a tuple representing (x,y) coordinates of a floor cell is a
    bottom left corner in the Warehouse.

    For example:
        #  \n
        #X \n
        ###\n

    @param floor_coords:
        A tuple representing (x,y) coordinates of a floor cell in the Warehouse

    @param all_floor_coords:
        A list of (x, y) tuples representing coordinates of all the inner floor
        cells of the Warehouse

    @return:
        True if the specified (x,y) coordinates tuple is a bottom left corner cell
        in the Warehouse.
        False otherswise.

    """
    return (
        has_up_neighbour(floor_coords, all_floor_coords)
        and has_right_neighbour(floor_coords, all_floor_coords)
        and not has_down_neighbour(floor_coords, all_floor_coords)
        and not has_left_neighbour(floor_coords, all_floor_coords)
    )


def is_bot_right_corner(floor_coords, all_floor_coords):
    """
    Determine if a tuple representing (x,y) coordinates of a floor cell is a
    bottom right corner in the Warehouse.

    For example:
          #\n
         X#\n
        ###\n

    @param floor_coords:
        A tuple representing (x,y) coordinates of a floor cell in the Warehouse

    @param all_floor_coords:
        A list of (x, y) tuples representing coordinates of all the inner floor
        cells of the Warehouse

    @return:
        True if the specified (x,y) coordinates tuple is a bottom right corner cell
        in the Warehouse.
        False otherswise.

    """
    return (
        has_up_neighbour(floor_coords, all_floor_coords)
        and has_left_neighbour(floor_coords, all_floor_coords)
        and not has_down_neighbour(floor_coords, all_floor_coords)
        and not has_right_neighbour(floor_coords, all_floor_coords)
    )


def is_top_left_corner(floor_coords, all_floor_coords):
    """
    Determine if a tuple representing (x,y) coordinates of a floor cell is a
    top left corner in the Warehouse.

    For example:
        ###\n
        #X \n
        #  \n

    @param floor_coords:
        A tuple representing (x,y) coordinates of a floor cell in the Warehouse

    @param all_floor_coords:
        A list of (x, y) tuples representing coordinates of all the inner floor
        cells of the Warehouse

    @return:
        True if the specified (x,y) coordinates tuple is a top left corner cell
        in the Warehouse.
        False otherswise.

    """
    return (
        has_right_neighbour(floor_coords, all_floor_coords)
        and has_down_neighbour(floor_coords, all_floor_coords)
        and not has_left_neighbour(floor_coords, all_floor_coords)
        and not has_up_neighbour(floor_coords, all_floor_coords)
    )


def is_top_right_corner(floor_coords, all_floor_coords):
    """
    Determine if a tuple representing (x,y) coordinates of a floor cell is a
    top right corner in the Warehouse.

    For example:
        ###\n
         X#\n
         ##\n

    @param floor_coords:
        A tuple representing (x,y) coordinates of a floor cell in the Warehouse

    @param all_floor_coords:
        A list of (x, y) tuples representing coordinates of all the inner floor
        cells of the Warehouse

    @return:
        True if the specified (x,y) coordinates tuple is a top right corner cell
        in the Warehouse.
        False otherswise.

    """
    return (
        has_left_neighbour(floor_coords, all_floor_coords)
        and has_down_neighbour(floor_coords, all_floor_coords)
        and not has_right_neighbour(floor_coords, all_floor_coords)
        and not has_up_neighbour(floor_coords, all_floor_coords)
    )


def is_caved_corner(floor_coords, all_floor_coords):
    """
    Determine if a tuple representing (x,y) coordinates of a floor cell is a
    caved corner in the Warehouse.

    For example:

        ###\n
        #X\n
        ###\n

    @param floor_coords:
        A tuple representing (x,y) coordinates of a floor cell in the Warehouse

    @param all_floor_coords:
        A list of (x, y) tuples representing coordinates of all the inner floor
        cells of the Warehouse

    @return:
        True if the specified (x,y) coordinates tuple is a caved corner in
        the Warehouse.
        False otherswise.

    """
    floor_neighbour_count = 0
    if has_down_neighbour(floor_coords, all_floor_coords):
        floor_neighbour_count += 1
    if has_left_neighbour(floor_coords, all_floor_coords):
        floor_neighbour_count += 1
    if has_up_neighbour(floor_coords, all_floor_coords):
        floor_neighbour_count += 1
    if has_right_neighbour(floor_coords, all_floor_coords):
        floor_neighbour_count += 1

    return floor_neighbour_count == 1


def is_corner(floor_coords, all_floor_coords):
    """
    Determine if a tuple representing (x,y) coordinates of a floor cell is a
    corner in the Warehouse.

    @param floor_coords:
        A tuple representing (x,y) coordinates of a floor cell in the Warehouse

    @param all_floor_coords:
        A list of (x, y) tuples representing coordinates of all the inner floor
        cells of the Warehouse

    @return:
        True if the specified (x,y) coordinates tuple is a corner floor cell
        of the Warehouse.
        False otherswise.

    """
    return (
        is_bot_left_corner(floor_coords, all_floor_coords)
        or is_bot_right_corner(floor_coords, all_floor_coords)
        or is_top_right_corner(floor_coords, all_floor_coords)
        or is_top_left_corner(floor_coords, all_floor_coords)
        or is_caved_corner(floor_coords, all_floor_coords)
    )


def get_floor_area_corners(all_floor_coords):
    """
    Return the set of (x, y) tuples which contains pairs of coordinates of all the corners in
    the inner floor area of the warehouse

    @param all_floor_coords:
        A list of (x, y) tuples representing coordinates of all the inner floor
        cells of the Warehouse

    @return:
        A set of (x, y) tuples containing coordinates of all the corners in
        the inner floor area of the warehouse

    """
    return {
        floor_coords
        for floor_coords in all_floor_coords
        if is_corner(floor_coords, all_floor_coords)
    }


def are_inline_vert(coords_a, coords_b):
    """
    Determine if a tuple representing (x,y) coordinates is vertically in line with
    another (x,y) coordinates tuple

    @param coords_a:
        A tuple representing (x,y) coordinates

    @param coords_b:
        Another tuple representing (x,y) coordinates

    @return:
       True if the two specified coordinates tuples are vertically in line,
       meaning the x coordinates are the same between the tuples.
       False otherswise.

    """
    return coords_a[0] == coords_b[0]


def are_inline_hor(coords_a, coords_b):
    """
    Determine if a tuple representing (x,y) coordinates is horizontally in line with
    another (x,y) coordinates tuple.

    @param coords_a:
        A tuple representing (x,y) coordinates

    @param coords_b:
        Another tuple representing (x,y) coordinates

    @return:
       True if the two specified coordinates tuples are horizontally in line,
       meaning the y coordinates are the same between the tuples.
       False otherswise.

    """
    return coords_a[1] == coords_b[1]


def space_in_line(coords_a, coords_b):
    """
    Returns the set of (x, y) tuples which contains pairs of coordinates in between
    the two specified (x, y) coordinates tuples

    @param coords_a:
        A tuple representing (x,y) coordinates

    @param coords_b:
        Another tuple representing (x,y) coordinates

    @return:
        A set containing (x, y) tuples which represents pairs of coordinates in between
        the two specified (x, y) coordinates tuples, if they are either
        vertically or horizontally in line.

        None, if coords_a and coords_b are not in line, or have no space in between
        (i.e. directly neighbours).

    """
    if are_inline_hor(coords_a, coords_b):
        # y same if aligned horizontally
        y = coords_a[1]
        # Inclusive
        x_start = min(coords_a[0], coords_b[0]) + 1
        # Exclusive
        x_end = max(coords_a[0], coords_b[0])
        if x_start < x_end:
            return {(x, y) for x in range(x_start, x_end)}

    elif are_inline_vert(coords_a, coords_b):
        # x same if aligned vertically
        x = coords_a[0]
        # Inclusive
        y_start = min(coords_a[1], coords_b[1]) + 1
        # Exclusive
        y_end = max(coords_a[1], coords_b[1])
        if y_start < y_end:
            return {(x, y) for y in range(y_start, y_end)}

    # Not in line, or have no space in between
    return None


def taboo_cells_coords(warehouse):
    """
    Returns the set of (x, y) tuples which contains pairs of coordinates of all taboo cells
    in the Warehouse object specified.

    Use only the following rules to determine the taboo cells;
        Rule 1: if a cell is a corner and not a target, then it is a taboo cell.
        Rule 2: all the cells between two corners along a wall are taboo if none of
        these cells is a target.

    @param warehouse:
        a Warehouse object with the worker inside the warehouse

    @return:
       A set of (x, y) tuples containing coordinates of the taboo cells
       in the Warehouse object specified.

    """
    floor_coords = get_floor_cells(Floor(warehouse))

    # Rule 1 coordinates set
    rule_1 = {
        cell
        for cell in get_floor_area_corners(floor_coords)
        if cell not in warehouse.targets
    }

    rule_2 = set()

    # Test every corner combination in rule 1
    for corner_combination in itertools.combinations(rule_1, 2):
        cell_a = corner_combination[0]
        cell_b = corner_combination[1]

        # Corners not aligned
        if not (are_inline_hor(cell_a, cell_b) or are_inline_vert(cell_a, cell_b)):
            continue

        # Corners are directly neighbours
        if not space_in_line(cell_a, cell_b):
            continue

        # Flag to determine whether a cell in between cell_a and cell_b does not comply with Rule 2.
        # If true, don't need to check the rest of the cells in between
        not_taboo_blocks = False
        for space in space_in_line(cell_a, cell_b):
            # Check Rule 2 fails
            if (space in warehouse.targets) or not (
                space in floor_coords and has_any_neighbour(space, warehouse.walls)
            ):
                not_taboo_blocks = True
                break

        # Continue to the next combination of corners if flag is true
        if not_taboo_blocks:
            continue

        # Else, the spaces all comply with Rule 2, hence add to set
        for space in space_in_line(cell_a, cell_b):
            rule_2.add(space)

    return rule_1 | rule_2


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Warehouse Updating Utility


def warehouse_result(warehouse, action, taboo_coords):
    """
    Return the resulting Warehouse object after applying the specified action
    on the specified Warehouse object. This function takes into consideration
    the illegal actions and taboo cells.

    Important notes:
        If an illegal action has been applied
            Or
        If the taboo cell coordinates specified is not an empty list, and
        the action results a box to be located on a taboo cell

        Then the original Warehouse object is returned

        Illegal actions have been defined as
        - if the agent tries to go into a wall
        - if the agent tries to push a box into a wall
        - if the agent tries to push two boxes at the same time

    @param warehouse:
        A valid Warehouse object

    @param action:
        A string of a single action.
        For example, one of the following: 'Up', 'Down', 'Left','Right'

    @param taboo_coords:
        A set of (x,y) tuples containing pairs of coordinates of the taboo cells
        in the corresponding Warehouse

    @return:
        A valid Warehouse object after applying the specified action
        on the input Warehouse object.

    """
    current_warehouse = warehouse

    worker_coords = current_warehouse.worker
    next_worker_coords = worker_coords

    # Apply specified action on worker coords accordingly
    if action == "Up":
        next_worker_coords = move_up(worker_coords)
    if action == "Down":
        next_worker_coords = move_down(worker_coords)
    if action == "Left":
        next_worker_coords = move_left(worker_coords)
    if action == "Right":
        next_worker_coords = move_right(worker_coords)

    # Worker attempts to go into wall, return warehouse unchanged
    if next_worker_coords in current_warehouse.walls:
        return warehouse

    boxes_coords = list(current_warehouse.boxes)

    # Worker attempts push a box
    if next_worker_coords in boxes_coords:
        # Identify the index of the box being pushed
        idx = boxes_coords.index(next_worker_coords)

        # Apply speficied actions on a box coords accordingly
        if action == "Up":
            boxes_coords[idx] = move_up(boxes_coords[idx])
        if action == "Down":
            boxes_coords[idx] = move_down(boxes_coords[idx])
        if action == "Left":
            boxes_coords[idx] = move_left(boxes_coords[idx])
        if action == "Right":
            boxes_coords[idx] = move_right(boxes_coords[idx])

        # The box attempts to go into wall, or into taboo cells, return warehouse unchanged
        if (boxes_coords[idx] in current_warehouse.walls) or (
            boxes_coords[idx] in taboo_coords
        ):
            return warehouse

        new_boxes_coords_set = set(boxes_coords)

        # Two boxes attempts to be pushed at the same time, return warehouse unchanged
        if len(new_boxes_coords_set) != len(boxes_coords):
            return warehouse

    # The actions is successful, copy and update the Warehouse object accordingly
    current_warehouse = current_warehouse.copy(
        worker=next_worker_coords, boxes=tuple(boxes_coords)
    )

    return current_warehouse


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Floor Problem


class Floor(search.Problem):
    """
    An instance of the class 'Floor' represents the inner floor searching/
    flooding Problem

    An instance contains information about the warehouse, but only considering
    the worker and the walls
    """

    def __init__(self, warehouse):
        self.warehouse = sokoban.Warehouse.copy(warehouse)
        self.initial = warehouse.worker
        # No goal to find, the problem is used to explored the inner floor
        self.goal = None

    def actions(self, state):
        """
        Return the list of actions that can be executed in the given state.
        """
        available_actions = []

        # Check if action results in the worker on top of a wall.
        # Target and box coords are ignored for this problem of
        # finding the internal floor area.
        if move_right(state) not in self.warehouse.walls:
            available_actions.append("Right")
        if move_left(state) not in self.warehouse.walls:
            available_actions.append("Left")
        if move_down(state) not in self.warehouse.walls:
            available_actions.append("Down")
        if move_up(state) not in self.warehouse.walls:
            available_actions.append("Up")

        return available_actions

    def result(self, state, action):
        """
        Return the state that results from executing the given
        action in the given state. The action must be one of
        self.actions(state).
        """
        # defensive programming
        assert action in self.actions(state)

        # Apply the actions accordingly
        if action == "Up":
            next_state = move_up(state)
        elif action == "Down":
            next_state = move_down(state)
        elif action == "Left":
            next_state = move_left(state)
        elif action == "Right":
            next_state = move_right(state)

        return next_state


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Floor Search


def floor_graph_search(problem, frontier):
    """
    A modified version of the search.graph_search. The return statement after
    the frontier has been exhausted is modified to the explored set rather than
    None. Hence, the graph search returns the coordinates of all internal floor
    cells to a warehouse.

    Assuming the state of the problem is the coordinates of the worker position,
    and that goal_state has been defined as None.

    """

    assert isinstance(problem, search.Problem)
    frontier.append(search.Node(problem.initial))
    explored = set()
    while frontier:
        node = frontier.pop()
        if problem.goal_test(node.state):
            return node
        explored.add(node.state)
        frontier.extend(
            child
            for child in node.expand(problem)
            if child.state not in explored and child not in frontier
        )
    return explored


def get_floor_cells(floor):
    """
    Returns the set of (x, y) tuples which contains pairs of coordinates of all internal
    floor cells to a warehouse using a modified version of Graph BFS

    @param floor:
        A Floor object (inherits search.Problem), representing the subject to this search problem.

    @return:
        A set of (x, y) tuples representing pairs of coordinates of all internal
        floor cells to a warehouse.

    """
    return floor_graph_search(floor, search.FIFOQueue())  # Graph search version of BFS


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SokobanPuzzel Problem
class SokobanPuzzle(search.Problem):
    """
    An instance of the class 'SokobanPuzzle' represents a Sokoban puzzle.
    An instance contains information about the walls, the targets, the boxes
    and the worker.

    """

    def __init__(self, warehouse):
        self.warehouse = warehouse
        # State represented by worker and boxes coordinates
        self.initial = (warehouse.worker, tuple(warehouse.boxes))
        # Store the box to weight mapping
        self.box_weight_map = dict(zip(warehouse.boxes, warehouse.weights))
        # Identify taboo cells
        self.taboo = taboo_cells_coords(warehouse)

    def actions(self, state):
        """
        Return the list of actions that can be executed in the given state.
        """
        available_actions = []
        self.warehouse = self.warehouse.copy(worker=state[0], boxes=list(state[1]))

        # If the action doesn't result in the same warehouse (using warehose_result).
        # (Checking for illegal actions and moving box into taboo cells).
        if str(self.warehouse) != str(
            warehouse_result(self.warehouse, "Up", self.taboo)
        ):
            available_actions.append("Up")
        if str(self.warehouse) != str(
            warehouse_result(self.warehouse, "Down", self.taboo)
        ):
            available_actions.append("Down")
        if str(self.warehouse) != str(
            warehouse_result(self.warehouse, "Left", self.taboo)
        ):
            available_actions.append("Left")
        if str(self.warehouse) != str(
            warehouse_result(self.warehouse, "Right", self.taboo)
        ):
            available_actions.append("Right")

        return available_actions

    def result(self, state, action):
        """
        Return the state that results from executing the given
        action in the given state. The action must be one of
        self.actions(state).
        """
        # defensive programming
        assert action in self.actions(state)

        # Update the Warehouse object
        self.warehouse = warehouse_result(self.warehouse, action, self.taboo)
        # Update the box to weight mapping
        self.box_weight_map = dict(zip(self.warehouse.boxes, self.warehouse.weights))

        return (self.warehouse.worker, tuple(self.warehouse.boxes))

    def goal_test(self, state):
        """
        Return True if the state is a goal. Compares the equality
        of the set of boxes coordinates in the current state, to the
        set of target coordinates in the Warehouse object.
        """
        return set(state[1]) == set((self.warehouse.targets))

    def path_cost(self, c, state1, action, state2):
        """
        Return the cost of a solution path that arrives at state2 from
        state1 via action, assuming cost c to get up to state1. The cost
        of an action is 1 + weight of the box being pushed (if any).
        """
        state1_worker, state1_boxes = state1
        state2_worker, state2_boxes = state2

        # Worker coordinates unchanged, return original cost
        if state1_worker == state2_worker:
            return c

        # Worker has moved, add 1 for worker term
        c += 1

        moved_box = set(state2_boxes) - set(state1_boxes)
        # Box has not been moved, return worker term only
        if len(moved_box) <= 0:
            return c

        # Box has been moved, account for box weight
        search_key = tuple(moved_box)
        moved_box_weight = self.box_weight_map[search_key[0]]
        c += moved_box_weight

        return c

    def h(self, node):
        """
        Heuristic function h(n) for the SokobanPuzzel

        Consists of a worker term and a boxes term

        The worker term is determined by the minimum Manhanttan Distance from the
        worker to the closest box.

        The boxes term is determined by the minimum summed Manhanttan Distance
        between permutations of box to target pairs, taking into account the weight
        of the boxes.
        """
        warehouse = self.warehouse
        worker, boxes, targets = warehouse.worker, warehouse.boxes, warehouse.targets

        worker_term = get_min_manhattan_distance(worker, boxes)

        # Set initial value to large number in order to calculate minimum
        boxes_term = 100000000

        for perm in itertools.permutations([i for i in range(len(boxes))]):
            # Counters (same range as i for i in range(len(boxes)))
            # to use in list comprehension
            counter = itertools.count(0)
            counter_2 = itertools.count(0)
            # Calculate the current permutation's box term
            perm_boxes_term = sum(
                [
                    get_manhattan_distance(boxes[next(counter)], targets[j])
                    + (self.box_weight_map[boxes[next(counter_2)]])
                    for j in perm
                ]
            )
            boxes_term = min(boxes_term, perm_boxes_term)

        return worker_term + boxes_term


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sanity Check


def taboo_cells(warehouse):
    """
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

    @return:
       A string representing the warehouse with only the wall cells marked with
       a '#' and the taboo cells marked with a 'X'.
       The returned string should NOT have marks for the worker, the targets,
       and the boxes.
    """
    floor_coords = get_floor_cells(Floor(warehouse))

    # Rule 1 coordinates set
    rule_1 = {
        cell
        for cell in get_floor_area_corners(floor_coords)
        if cell not in warehouse.targets
    }

    rule_2 = set()

    # Test every corner combination in rule 1
    for corner_combination in itertools.combinations(rule_1, 2):
        cell_a = corner_combination[0]
        cell_b = corner_combination[1]

        # Corners not aligned
        if not (are_inline_hor(cell_a, cell_b) or are_inline_vert(cell_a, cell_b)):
            continue

        # Corners are directly neighbours
        if not space_in_line(cell_a, cell_b):
            continue

        # Flag to determine whether a cell in between cell_a and cell_b does not comply with Rule 2.
        # If true, don't need to check the rest of the cells in between
        not_taboo_blocks = False
        for space in space_in_line(cell_a, cell_b):
            # Check Rule 2 fails
            if (space in warehouse.targets) or not (
                space in floor_coords and has_any_neighbour(space, warehouse.walls)
            ):
                not_taboo_blocks = True
                break

        # Continue to the next combination of corners if flag is true
        if not_taboo_blocks:
            continue

        # Else, the spaces all comply with Rule 2, hence add to set
        for space in space_in_line(cell_a, cell_b):
            rule_2.add(space)

    # Prepare to string results
    X, Y = zip(*warehouse.walls)
    x_size, y_size = 1 + max(X), 1 + max(Y)

    vis = [[" "] * x_size for y in range(y_size)]
    for (x, y) in warehouse.walls:
        vis[y][x] = "#"
    for (x, y) in rule_1 | rule_2:
        vis[y][x] = "X"

    return "\n".join(["".join(line) for line in vis])


def check_elem_action_seq(warehouse, action_seq):
    """

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
    """

    current_warehouse = sokoban.Warehouse.copy(warehouse)

    for action in action_seq:
        # If the action doesn't result in the same warehouse (using warehose_result).
        # (Checking for illegal actions only, but allows boxes to be pushed on taboo cells).
        if str(current_warehouse) == str(
            warehouse_result(current_warehouse, action, [])
        ):
            return "Impossible"

        current_warehouse = warehouse_result(current_warehouse, action, [])

    return str(current_warehouse)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sokoban Solver


def trace_actions(goal_node):
    """
    Return the list of actions required to transform the initial state
    to the goal state specified in the input Node object

    The action string must be in the following: 'Up', 'Down', 'Left','Right'.

    @param
        goal_node: a Node object

    @return:
        A list of action strings to be applied in order to get from the initial
        state to the state in the goal_node

    """
    path = goal_node.path()
    return [node.action for node in path if node.action]


def solve_weighted_sokoban(warehouse):
    """
    This function analyses the given warehouse.
    It returns the two items. The first item is an action sequence solution.
    The second item is the total cost of this action sequence.

    @param
     warehouse: a valid Warehouse object

    @return:

        If puzzle cannot be solved
            return 'Impossible', None

        If a solution was found,
            return S, C
            where S is a list of actions that solves
            the given puzzle coded with 'Left', 'Right', 'Up', 'Down'
            For example, ['Left', 'Down', Down','Right', 'Up', 'Down']
            If the puzzle is already in a goal state, simply return []
            C is the total cost of the action sequence C

    """

    sokoban_puzzle = SokobanPuzzle(warehouse)

    goal_node = search.astar_graph_search(sokoban_puzzle, sokoban_puzzle.h)

    # Puzzle is already in a goal state
    if sokoban_puzzle.goal_test(sokoban_puzzle.initial):
        return []
    # Puzzel cannot be solved
    elif (
        goal_node is None
        or check_elem_action_seq(warehouse, trace_actions(goal_node)) == "Impossible"
    ):
        return ["Impossible", None]
    # Solution is found
    else:
        return [trace_actions(goal_node), goal_node.path_cost]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
