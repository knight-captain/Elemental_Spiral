import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from config import ANCHOR_GROUP, ANCHOR_ANGLE, BEAD_SIZE, SPIRAL_COLOR, GROUP_LINE_WIDTH, FONT_SIZE

# 1. Canvas preparation
def prepare_canvas(figsize=(10, 10)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig, ax

# 2. Color map builder
def make_group_colormap():
    anchor_colors = [
        "#0000ff", "#00ffaa", "#00ff00",
        "#ff9900", "#ff0000", "#5b2c83"
    ] #yellow (#ffff00) is too light
    cmap = LinearSegmentedColormap.from_list("group_cmap", anchor_colors)
    norm = plt.Normalize(1, 18)
    return cmap, norm

def group_color(group, cmap, norm):
    if group is None:
        return "#7a7a7a"
    if group == 18:
        return cmap(norm(18))
    return cmap(norm(group))

# 3. Spiral + beads + group lines
def draw_spiral_and_beads(elements, positions, ax):
    ordered = sorted(elements, key=lambda e: e["atomic_number"])

    # Build colormap
    cmap, norm = make_group_colormap()

    # Extract theta_raw = radius
    thetas = {Z: positions[Z][0] for Z in positions}

    # --- FIND ANCHOR GROUP ORIENTATION ---
    anchor_members = [e for e in ordered if e["group"] == ANCHOR_GROUP]

    if anchor_members:
        anchor_Z = anchor_members[0]["atomic_number"]
        anchor_angle_raw = positions[anchor_Z][0]
    else:
        anchor_angle_raw = 0.0  # fallback

    # Rotation needed to align anchor group to ANCHOR_ANGLE
    rotation = ANCHOR_ANGLE - anchor_angle_raw

    # --- DRAW SPIRAL CURVE ---
    max_theta = max(thetas.values())
    spiral_thetas = np.linspace(0, max_theta, 2000)

    spiral_x = spiral_thetas * np.cos(spiral_thetas + rotation)
    spiral_y = spiral_thetas * np.sin(spiral_thetas + rotation)

    ax.plot(spiral_x, spiral_y, color=SPIRAL_COLOR, linewidth=2, zorder=1)

    # --- GROUP LINES ---
    grouped = {}
    for item in ordered:
        group = item["group"]
        if group is None:
            continue
        grouped.setdefault(group, []).append(item)

    for group in sorted(grouped):
        members = grouped[group]
        color = group_color(group, cmap, norm)

        for left, right in zip(members[:-1], members[1:]):
            rL, thL_raw = positions[left["atomic_number"]]
            rR, thR_raw = positions[right["atomic_number"]]

            thL = thL_raw + rotation
            thR = thR_raw + rotation

            lx = rL * np.cos(thL)
            ly = rL * np.sin(thL)
            rx = rR * np.cos(thR)
            ry = rR * np.sin(thR)

            ax.plot([lx, rx], [ly, ry],
                    color=color, lw=1.3, alpha=0.9, zorder=3)

    # --- BEADS + LABELS ---
    for item in ordered:
        Z = item["atomic_number"]
        r, th_raw = positions[Z]

        th = th_raw + rotation
        px = r * np.cos(th)
        py = r * np.sin(th)

        group = item["group"]
        color = group_color(group, cmap, norm)

        ax.scatter(px, py, s=BEAD_SIZE, color=color,
                   edgecolors="black", linewidths=GROUP_LINE_WIDTH, zorder=4)
        ax.text(px, py, item["atomic_name"],
                ha="center", va="center",
                fontsize=FONT_SIZE, color="white", zorder=5)

    ax.set_title("Spiral Table of the Elements", fontsize=14, pad=18)
