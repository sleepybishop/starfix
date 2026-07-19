import numpy as np
from PIL import Image

width, height = 800, 600
img = np.zeros((height, width), dtype=np.uint8)
pole_x, pole_y = 400, 300

# Draw concentric arcs
for radius in [50, 100, 150, 200, 250, 300]:
    for angle in np.linspace(0, np.pi/2, 100):
        x = int(pole_x + radius * np.cos(angle))
        y = int(pole_y + radius * np.sin(angle))
        if 0 <= x < width and 0 <= y < height:
            img[y, x] = 255
            # Add some thickness
            if 0 <= x+1 < width: img[y, x+1] = 255
            if 0 <= y+1 < height: img[y+1, x] = 255

im = Image.fromarray(img)
im.save('test_trails.pgm')
