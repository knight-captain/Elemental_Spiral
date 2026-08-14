import numpy as np
from config import ANCHOR_GROUP

def is_bounded_segment(elements, segment):

    if segment[-1] == elements[-1]["atomic_number"] and not is_anchor(elements[segment[0]]):
        #if the last element in a segment is the last element in the table, then it's the last segment
        #if the last segment is just an anchor, then it is trivially bounded
        #otherwise, it's unbounded by an anchor element, and needs special treatment
        return False

    return True

def is_anchor(element):
    """
    Return True if the element is an anchor.
    - Hydrogen (Z=1) is always an anchor.
    - If ANCHOR_GROUP is None: only Hydrogen anchors.
    - Otherwise: any element whose group == ANCHOR_GROUP anchors.
    """
    Z = element["atomic_number"]
    group = element["group"]

    if Z == 1:
        return True

    if ANCHOR_GROUP is None or ANCHOR_GROUP == 0:
        return False

    return group == ANCHOR_GROUP

def compute_segments(ordered):
    """
    Build spiral segments based on anchor elements.
    - Hydrogen always starts the first segment.
    - If ANCHOR_GROUP is None: one giant segment.
    - Otherwise: break segments at each element whose group == ANCHOR_GROUP.
    """

    segments = []
    current = []

    for e in ordered:
        Z = e["atomic_number"]
        group = e["group"]

        # If this element is an anchor (group == ANCHOR_GROUP)
        if group == ANCHOR_GROUP:
            # Close previous segment
            if current:
                segments.append(current)
            # Start new segment with this anchor
            current = [Z]
        else:
            current.append(Z)

    # Append final segment
    if current:
        segments.append(current)

    return segments

def find_ideal_spacing(segment):
    """
    Compute ideal angular spacing for a segment (turn).
    - Normal turn: d_ideal = 2π / len(segment)
    """
    count = len(segment)
    
    d_ideal = (2 * np.pi) / count if count > 1 else 0.0

    # Fallback
    return d_ideal