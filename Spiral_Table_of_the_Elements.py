import csv
from pathlib import Path as PathlibPath

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from Draw_Spiral import *
from Position_Elements import *

from config import ITERATIONS, STEP

DATA_FILE = PathlibPath(__file__).with_name("elements_list.csv")

def load_elements(csv_path: PathlibPath):
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = []
        for row in reader:
            atomic_number = int(row["AtomicNumber"])
            atomic_name = row.get("AtomicName") or row.get("Name")
            group = float(row.get("Group", "")) if str(row.get("Group", "")).strip() else None
            starter_angle = float(row.get("angle", 0.0)) #artifact of older system
            rows.append({
                "atomic_number": atomic_number,
                "atomic_name": atomic_name,
                "group": group
            })

    rows.sort(key=lambda item: item["atomic_number"])
    return rows

def main():
    elements_raw = load_elements(DATA_FILE)
    if not elements_raw:
        raise RuntimeError(f"No element rows were found in {DATA_FILE}")

    elements = sorted(elements_raw, key=lambda e: e["atomic_number"])

    # Build initial positions
    positions = build_spiral_positions(elements)
    segments = compute_segments(elements)

    # Live relaxation loop
    fig, ax = prepare_canvas()

    # Initial draw
    draw_spiral_and_beads(elements, positions, ax)
    plt.pause(STEP)

    # Relax in small increments

    for _ in range(ITERATIONS):
        # Apply a small amount of relaxation
        positions = relax_spiral(elements, positions)

        # Clear the axes for redraw
        ax.cla()
        ax.set_aspect("equal")
        ax.axis("off")

        # Redraw updated spiral
        draw_spiral_and_beads(elements, positions, ax)

        # Allow UI to update
        plt.pause(STEP)

    # Final save
    plt.tight_layout()
    plt.show()
    fig.savefig("periodic_table_spiral.png", dpi=220, bbox_inches="tight")


if __name__ == "__main__":
    main()
