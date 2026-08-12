import numpy as np

# The group to anchor the orientation and segments on
ANCHOR_GROUP = 18  # Default anchor group for segmentation = 18

# Spiral geometry defaults
SPIRAL_START_ANGLE = np.pi  # starting angle of the spiral, or how far He is from H
ANCHOR_ANGLE = np.pi / 2 #what to orient the anchor group to

# Force weights
REPULSION_WEIGHT = 1
ATTRACTION_WEIGHT = 50

# Drawing & Force weights
ITERATIONS = 100
STEP = 0.01

# Drawing defaults
BEAD_SIZE = 400
FONT_SIZE = 12
SPIRAL_COLOR = "lightgray"
GROUP_LINE_WIDTH = 1.3