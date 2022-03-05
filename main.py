import tkinter as tk
import threading as th
from time import sleep

from HexaGOL import *

size = 14
fen = tk.Tk()
can = tk.Canvas(fen, width=800, height=800)
can.pack(side=tk.BOTTOM)
tiles = set(map(lambda x: HEXVector(x[0], x[1], x[2]), [(0, 0, 0), (1, 0, -1), (0, 1, -1)]))
eng = HEXEngine((HEXVector(-16, -16, -16), HEXVector(16, 16, 16)), origin=Coords(400, 400), scale=size,
                living_cells=tiles)


def next_eng():
    candidates = eng.canditates()
    eng.next()
    grid = eng.get_grid()
    for tile in map(lambda hex: grid[hex], candidates):
        if tile.alive:
            filler = "black"
            outer = "white"
        else:
            filler = "white"
            outer = "black"
        can.itemconfig(can.find_closest(*tile.mean_coords), fill=filler, outline=outer)
        # can.create_text(*tile.mean_coords, text=tile.HEX, fill=outer, font=("Helvetica", 6))


for key, tile in eng.get_grid().items():
    if tile.alive:
        filler = "black"
        outer = "white"
    else:
        filler = "white"
        outer = "black"
    can.create_polygon(*[element for couple in tile.coords for element in couple], tag = f"{key.H};{key.E};{key.X}", fill=filler, outline=outer)

butt = tk.Button(fen, text="Next", command=next_eng)
butt.pack(side=tk.TOP)
fen.bind("<Return>", lambda _: next_eng())

def write_cell(event):
    canitem = can.find_closest(event.x, event.y)
    item = can.itemcget(canitem, "tag")
    item = item.split(' ')[0][1:-1] if "current" in item else item
    hexv = HEXVector(*[int(x) for x in item.split(';')])
    eng.flip_cell(hexv)
    if eng[hexv].alive:
        filler = "black"
        outer = "white"
    else:
        filler = "white"
        outer = "black"
    can.itemconfig(canitem, fill=filler, outline=outer)
can.bind("<Button-1>", write_cell)

text = """
The neighbours are calculated with a tier system : tier 1 neighbours are contact cells, tier 2 neighbours are
cells which points forms a David star. The "number" of alive cells is the sum of alive tier 1 neighbours + the
sum of alive tier 2 neighbours weighted 0.3;
- Births: Each dead cell with between 2 excluded and 3 included alive neighbors will become live in the next generation.
- Death by isolation: Each live cell with strictly less that 2 alive neighbors will die in the next generation.
- Death by overcrowding: Each live cell with 3 or more alive neighbors will die in the next generation.
- Survival: Each live cell with between 2 included and 3 excluded alive neighbors will remain alive for the next generation.
"""
rules = tk.Label(fen, text=text)
rules.pack(side=tk.BOTTOM)
fen.mainloop()
