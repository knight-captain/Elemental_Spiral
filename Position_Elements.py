import numpy as np

from config import (
    SPIRAL_START_ANGLE,
    STEP,
    REPULSION_WEIGHT,
    SPRING_GRATIENT,
    ATTRACTION_WEIGHT
    )

from utils import *
from Apply_Forces import spring_force, nearest_neighbor_attraction

def relax_spiral(elements, positions):
    """
    Relax the spiral using:
    - spring_force (repulsion + attraction toward ideal spacing)
    - nearest_neighbor_attraction (group-based attraction)
    Movement is constrained to the Archimedean spiral r = theta_raw.
    Anchors (H + anchor group) do not move.
    """

    # Build segments using anchor group
    segments = compute_segments(elements)

    # Extract theta_raw = radius
    old_thetas = {Z: positions[Z][0] for Z in positions}

    # Identify anchors (H + anchor group)
    anchors = {
        e["atomic_number"]
        for e in elements
        if is_anchor(e)
    }

    # for _ in range(ITERATIONS):
    new_thetas = dict(old_thetas)

    # Loop through segments (turns)
    for seg_index, segment in enumerate(segments):

        # Compute ideal spacing for this segment
        if is_bounded_segment(elements, segment):
            #normal segments get spaced evenly around the turn
            d_ideal = find_ideal_spacing(segment)
        else:
            #if the last segment is unbounded by an anchor, then it needs to match the previous segment's spacing
            d_ideal = find_ideal_spacing(segments[seg_index - 1])

        # Loop through elements in this segment
        for Z in segment:

            # Skip anchors
            if Z in anchors:
                continue

            # --- TOTAL TANGENTIAL FORCE ---
            Ft = 0.0

            # 1. Spring force (repulsion + attraction) towards evenly spacing elements around each segment
            raw_spring = spring_force(Z, positions, d_ideal, REPULSION_WEIGHT)
            if SPRING_GRATIENT > 0:
                gradient = SPRING_GRATIENT / len(segment)
            else:
                gradient = 1
            Ft += raw_spring * gradient

            # 2. Group attraction
            Ft += nearest_neighbor_attraction(Z, positions, elements, ATTRACTION_WEIGHT)

            # Convert tangential force to Δθ_raw
            dtheta = (Ft * STEP) / max(old_thetas[Z], d_ideal / 2)

            # Prevent overtaking
            ordered_Z = sorted(positions.keys())
            idx = ordered_Z.index(Z)

            #must be larger than the eps used in spring_force
            eps_order = 1e-2

            left_bound = (
                positions[ordered_Z[idx - 1]][0] + d_ideal / 2
                if idx > 0 else None
            )
            right_bound = (
                positions[ordered_Z[idx + 1]][0] - d_ideal / 2
                if idx < len(ordered_Z) - 1 else None
            )

            theta_new = old_thetas[Z] + dtheta

            if left_bound is not None:
                theta_new = max(theta_new, left_bound)
            if right_bound is not None:
                theta_new = min(theta_new, right_bound)

            new_thetas[Z] = theta_new

    # Update positions (r = θ_raw)
    for Z in positions:
        positions[Z] = (new_thetas[Z], new_thetas[Z])

    return positions


def build_spiral_positions(elements):
    if not elements:
        return {}

    positions = {}

    segments = compute_segments(elements)

    for seg_index, segment in enumerate(segments):

        if is_bounded_segment(elements, segment):
            #normal segments get spaced evenly around the turn
            d_ideal = find_ideal_spacing(segment)
        else:
            #if the last segment is unbounded by an anchor, then it needs to match the previous segment's spacing
            d_ideal = find_ideal_spacing(segments[seg_index - 1])

        # Place elements in this segment
        for i, Z in enumerate(segment):
            if Z == 1:
                # Hydrogen explicitly at r=0, θ=0
                positions[1] = (0.0, 0.0)
                continue

            #what each element's posotion should be
            theta_spiral = i * d_ideal

            #What segment/turn are we on?
            #to get He at SPIRAL_START_ANGLE, we need this
            theta_turn = (seg_index - (SPIRAL_START_ANGLE/(2 * np.pi))) * 2 * np.pi

            theta_raw = theta_spiral + theta_turn

            # Archimedean spiral: r = θ, but we don't want huge thetas for forces
            r = theta_spiral + theta_turn  
            theta_circular = (theta_raw % (2 * np.pi))

            positions[Z] = (theta_raw, r)

    return positions
