import numpy as np
from ggplot import *
#Extract Data
x,y = np.loadtxt('vicsek.txt',delimiter='\t',unpack=True)

p = ggplot(aes(x='date', y='beef'), data=meat)
p+geom_point()
