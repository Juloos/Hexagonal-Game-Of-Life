from math import dist, sqrt
from statistics import mean
from typing import NamedTuple, Tuple, Union, Set, Dict
from copy import deepcopy

HEX_X_RATIO = sqrt(3) / 2
HEX_Y_RATIO = 1 / 2
UNCERTAINTY = 1 / 1000


class GeometricError(Exception):
    pass


class HexError(Exception):
    pass


class Coords(NamedTuple):
    x: float
    y: float


class SixDoubleUplet(NamedTuple):
    e1: Coords
    e2: Coords
    e3: Coords
    e4: Coords
    e5: Coords
    e6: Coords


class HEXVector(NamedTuple):
    H: int  # increasing to the right
    E: int  # increasing to the lower left
    X: int  # increasing to the lower right

    def __add__(self, other):
        if other.__class__.__name__ != self.__class__.__name__:
            err = "Cannot add two vectors that are not `HEXVector`"
            raise HexError(err)
        return HEXVector(self.H + other.H, self.E + other.E, self.X + other.E)

    def __str__(self):
        return f"({self.H}, {self.E}, {self.X})"


class HEXCell:
    """Hexagonal representation of a cell with its coordinates in a plan and its living state"""

    def __init__(self, HEX: HEXVector, coords: SixDoubleUplet, alive: bool = 0):
        self.HEX = HEX
        self.coords = coords
        self.alive = alive  # Boolean representation of the cell, 0 means dead

        # Checking that the cell's coords are canonical
        if HEX.H + HEX.E + HEX.X != 0:
            err = "HEXCell's coordinates should be canonical (H + E + X == 0) to ensure unicity"
            raise HexError(err)

        # Checking coords are coherent relatively to a regular hexagone geometry
        last_side_length = dist((coords[0].x, coords[0].y), (coords[1].x, coords[1].y))
        self.mean_coords = (mean([x[0] for x in coords]), mean([y[1] for y in coords]))
        for i in range(6):
            x1, y1 = coords[i].x, coords[i].y
            x2, y2 = coords[(i + 1) % 6].x, coords[(i + 1) % 6].y
            side_length = dist((x1, y1), (x2, y2))
            dist_to_mean_point = dist(self.mean_coords, (x1, y1))

            if abs(last_side_length - side_length) >= UNCERTAINTY:
                err = "Not a regular hexagone, all sides are not of the same length"
                raise GeometricError(err)

            if abs(dist_to_mean_point - side_length) >= UNCERTAINTY:
                err = "Not a regular hexagone, all points have not the same distance to the geometric mean"
                raise GeometricError(err)

    def kill(self):
        """Kills the current cell, if alive"""
        self.alive = False

    def revive(self):
        """Revives the current cell, if dead"""
        self.alive = True

    def flip(self):
        """Invert the living state version of the current cell"""
        self.alive = not self.alive

    def dead(self):
        """Returns a dead version of the current cell"""
        return HEXCell(self.HEX, self.coords, False)

    def living(self):
        """Returns a living version of the current cell"""
        return HEXCell(self.HEX, self.coords, True)

    def flipped(self):
        """Returns the inversed living state version of the current cell"""
        return HEXCell(self.HEX, self.coords, not self.alive)


class HEXEngine:
    """Simulates Conway's Game Of Life with hexagonal cells !

    Here are the rules that apply to all cells at each round (from Conway's Game Of Life, remanipulated):
        The neighbours are calculated with a tier system : tier 1 neighbours are contact cells, tier 2 neighbours are
        cells which points forms a David star. The "number" of alive cells is the sum of alive tier 1 neighbours + the
        sum of alive tier 2 neighbours weighted 0.3;
    - Births: Each dead cell with between 2 excluded and 3 included alive neighbors will become live in the next generation.
    - Death by isolation: Each live cell with strictly less that 2 alive neighbors will die in the next generation.
    - Death by overcrowding: Each live cell with 3 or more alive neighbors will die in the next generation.
    - Survival: Each live cell with between 2 included and 3 excluded alive neighbors will remain alive for the next generation.
    """

    def __init__(self, world_size: Tuple[HEXVector, HEXVector], origin: Coords = Coords(0, 0), scale: float = 10,
                 living_cells: Set[HEXVector] = None):
        self._origin = origin
        self._scale = scale
        self._alive_cells = living_cells or set()

        self._grid: Dict[HEXVector, HEXCell] = dict()
        self._size = world_size
        _from, _to = world_size
        for h in range(_from.H, _to.H + 1):
            for e in range(_from.E, _to.E + 1):
                for x in range(_from.X, _to.X + 1):
                    if h + e + x == 0:  # The addition must be vanishing to check the canonical coordinate for each hex
                        self._grid[HEXVector(h, e, x)] = self.create_cell(HEXVector(h, e, x))
    
    def get_world_size(self) -> Tuple[HEXVector, HEXVector]:
        """Returns the size of the world space"""
        return self._size
    
    def get_world_origin(self) -> Coords:
        """Returns the origin coords of the world space"""
        return self._origin
    
    def get_world_scale(self) -> float:
        """Returns the scale of the world space"""
        return self._scale
    
    def get_living_cells(self) -> Set[HEXVector]:
        """Returns the set of vectors who's coordinates are the current living cells"""
        return self._alive_cells
    
    def get_grid(self) -> Dict[HEXVector, HEXCell]:
        """Returns a dict of `HEXCell` who's keys are the `HEXVector` of the engine"""
        return self._grid

    def create_cell(self, HEX: HEXVector) -> HEXCell:
        """Returns `HEXCell` object that correspond to the given cubic (HEX) coords"""
        posx = self._origin.x + self._scale * HEX_X_RATIO * (HEX.X - HEX.E)
        posy = self._origin.y + self._scale * (HEX_Y_RATIO * (HEX.E + HEX.X) - HEX.H)
        coords = SixDoubleUplet(Coords(posx                            , posy + self._scale              ),
                                Coords(posx + HEX_X_RATIO * self._scale, posy + HEX_Y_RATIO * self._scale),
                                Coords(posx + HEX_X_RATIO * self._scale, posy - HEX_Y_RATIO * self._scale),
                                Coords(posx                            , posy - self._scale              ),
                                Coords(posx - HEX_X_RATIO * self._scale, posy - HEX_Y_RATIO * self._scale),
                                Coords(posx - HEX_X_RATIO * self._scale, posy + HEX_Y_RATIO * self._scale)
                                )
        return HEXCell(HEX, coords, HEX in self._alive_cells)

    def flip_cell(self, cell: Union[HEXVector, HEXCell]) -> None:
        """Flips a cell's living state"""
        if isinstance(cell, HEXCell):
            cell = cell.HEX
        
        if cell in self._grid.keys():
            if self[cell].alive:
                self._alive_cells.remove(cell)
                self._grid[cell].kill()
            else:
                self._alive_cells.add(cell)
                self._grid[cell].revive()
        else:
            err = repr(cell)
            raise KeyError(err)

    def neighbours(self, cell: Union[HEXVector, HEXCell]):
        """Returns all neighbours of a cell, at most 12"""
        if isinstance(cell, HEXCell):
            cell = cell.HEX

        world = self._grid.keys()
        if cell in world:
            neighbours1 = list()
            neighbours2 = list()
            for h in range(cell.H - 2, cell.H + 3):
                for e in range(cell.E - 2, cell.E + 3):
                    for x in range(cell.X - 2, cell.X + 3):
                        if h + e + x == 0:
                            neighbor = HEXVector(h, e, x)
                            abs_hex = map(abs, (h - cell.H, e - cell.E, x - cell.X))  # ensuring there is 12 neighbours max
                            if cell != neighbor and neighbor in world and 1 in abs_hex:
                                if sum(abs_hex) == 4:  # 2nd tiers neighbours which forms a David star
                                    neighbours2.append(neighbor)
                                else:  # The rest is 1st tier
                                    neighbours1.append(neighbor)
            return neighbours1, neighbours2
        else:
            err = repr(cell)
            raise KeyError(err)
    
    def canditates(self) -> Set[HEXVector]:
        """Returns a set of candidate cells that will *at least* need attention when applying the ruleset"""
        # Gets all alive cells and their neighbours
        accounting_cells = set()
        for HEX in self._alive_cells:
            accounting_cells |= {*[i for row in self.neighbours(HEX) for i in row], HEX}
        return accounting_cells

    def apply_rules(self, cell: Union[HEXVector, HEXCell]):
        """Applies the rules of the Game Of Life for Hexagones to one cell and returns a `HEXCell` that represents the evolved cell"""
        if isinstance(cell, HEXVector):
            cell = self[cell]

        neighbours = self.neighbours(cell.HEX)
        alive_neighbours = 0
        for (tier, weight) in ((0, 1), (1, 0.3)):
            alive_neighbours += weight * sum(list(map(lambda neighbor: self[neighbor].alive, neighbours[tier])))
        if cell.alive:
            if 2 <= alive_neighbours < 3:
                return cell.living()  # Survival
            self._alive_cells.remove(cell.HEX)
            return cell.dead()  # Death by overpopulation or isolation
        elif 2 < alive_neighbours <= 3:
            self._alive_cells.add(cell.HEX)
            return cell.living()  # Birth
        return cell  # Dead

    def next(self):
        """Plays to the next round of the game"""
        new_grid = deepcopy(self._grid)
        for HEX in self.canditates():
            new_grid[HEX] = self.apply_rules(HEX)
        self._grid = new_grid

    def __getitem__(self, __k: HEXVector) -> HEXCell:
        return self._grid[__k]
