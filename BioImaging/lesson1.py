from skimage import io
import matplotlib.pyplot as plt

url = 'https://github.com/vdsukhov/bioimage-analysis-course/blob/main/data/images/dapi.png?raw=true'

img = io.imread(url)

plt.imshow(img, cmap='gray')
plt.axis("off")
plt.gcf().set_facecolor('black')