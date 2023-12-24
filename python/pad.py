import numpy as np
import sys
from PIL import Image


image_path = sys.argv[1]

img = Image.open(image_path)

w, h = img.size
print(f'Old width {w}, height {h}')


pad_left = (200 - w)//2
pad_right = 200 - w - pad_left

pad_above = (200 - h)//2
pad_below = 200 - h - pad_above

print(f'pad_left {pad_left}, pad_right {pad_right}')
print(f'pad_above {pad_above}, pad_below {pad_below}')

arr = np.array(img)

# Zero padding
# padded_arr = np.pad(arr, ((pad_above, pad_below),(pad_left, pad_right), (0,0)), 'constant')
# new_img = Image.fromarray(padded_arr)
# new_img.save(image_path + '_zero_padded.png')

# Edged padding
padded_arr = np.pad(arr, ((pad_above, pad_below),(pad_left, pad_right), (0,0)), 'edge')
new_img = Image.fromarray(padded_arr)
new_img.save(image_path + '_edge_padded.png')

print(f'final size {new_img.size}')
