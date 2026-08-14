import numpy as np

from utils import *

CHECK_ELEMENT = 0

def circ_diff(a, b):
    d = (a - b) % (2*np.pi)
    if d > np.pi:
        d -= 2*np.pi
    return d

def spring_force(Z, positions, d_ideal, REPULSION_WEIGHT, eps=1e-4):
    """
    Compute tangential spring force for element Z based on ideal spacing.
    Pure θ-space force: no x,y, no anchor correction.
    Uses left/right neighbors from positions (which must be ordered).
    """

    Ft = 0.0

    if Z == CHECK_ELEMENT:
        Fhe = 0
        Fb = 0

    # Extract ordered list of atomic numbers
    ordered_Z = list(positions.keys())
    ordered_Z.sort()
    idx = ordered_Z.index(Z)

    # Current element's theta
    theta_i = positions[Z][0]

    # --- LEFT NEIGHBOR ---
    if idx > 0:
        Z_left = ordered_Z[idx - 1]
        theta_left = positions[Z_left][0]

        d = abs(circ_diff(theta_i, theta_left))
        if d > eps:
            if d < d_ideal:
                # Repulsion (inverse square)
                Fmag = REPULSION_WEIGHT * (1.0 / (d * d))
            else:
                # Attraction (linear)
                Fmag = -REPULSION_WEIGHT * d * d

            Ft += Fmag * np.sign(theta_i - theta_left)
            if Z == CHECK_ELEMENT:
                Fhe = Fmag

    # --- RIGHT NEIGHBOR ---
    if idx < len(ordered_Z) - 1:
        Z_right = ordered_Z[idx + 1]
        theta_right = positions[Z_right][0]

        d = abs(circ_diff(theta_i, theta_right))
        if d > eps:                         # <-- REQUIRED
            if d < d_ideal:
                # Repulsion (inverse square)
                Fmag = REPULSION_WEIGHT * (1.0 / (d * d))
            else:
                # Attraction (linear)
                Fmag = -REPULSION_WEIGHT * d * d

            Ft += Fmag * np.sign(theta_i - theta_right)
            if Z == CHECK_ELEMENT:
                Fb = Fmag
    
    if Z == CHECK_ELEMENT:
            print(f"{Z} spring: next={Fb} prev={Fhe} -> {Fhe-Fb}CCW")
    return Ft

def nearest_neighbor_attraction(Z, positions, elements, ATTRACTION_WEIGHT):
    """
    Tangential attraction toward adjacent group neighbors.
    Z = atomic_number
    positions = dict[Z](theta, r)
    elements = list [of dicts{"atomic_number", "atomic_name", "group"}]
    """

    # Find this element's group
    element = next(e for e in elements if e["atomic_number"] == Z)
    group = element["group"]

    if Z == CHECK_ELEMENT:
        Fh = 0
        Fna = 0

    # No group → no attraction
    if group is None:
        return 0.0

    # Build list of group members in atomic-number order
    group_members = sorted([e["atomic_number"] for e in elements if e["group"] == group])

    # If only one element in group → no attraction
    if len(group_members) <= 1:
        return 0.0

    # Find index of Z within its group
    try:
        idx = group_members.index(Z)
    except ValueError:
        print("THIS SHOULD NEVER HAPPEN")
        return 0.0

    # Current element's theta
    theta_i = positions[Z][0] # was [0]

    Ft = 0.0

    # Previous group neighbor
    if idx > 0:
        Z_prev = group_members[idx - 1]
        theta_prev = positions[Z_prev][0]
        d = circ_diff(theta_prev, theta_i)
        Ft += ATTRACTION_WEIGHT * d
        if Z == CHECK_ELEMENT:
            Fh = ATTRACTION_WEIGHT * d

    # Next group neighbor
    if idx < len(group_members) - 1:
        Z_next = group_members[idx + 1]
        theta_next = positions[Z_next][0]
        d = circ_diff(theta_next, theta_i)
        Ft += ATTRACTION_WEIGHT * d
        if Z == CHECK_ELEMENT:
            Fna = ATTRACTION_WEIGHT * d

    if Z == CHECK_ELEMENT:
        print(f"{Z} group: in={Fh} out={Fna} -> {Fna+Fh}CCW")
    return Ft #tangential force to bring the element closer to its adjacent group members




