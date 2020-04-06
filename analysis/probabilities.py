import numpy as np
import matplotlib.pyplot as plt
plt.style.use('style.mplstyle')

# target thicknesses
t = [1,2,3,4,5,6,7,8,9,10,20]

# probability of proton being captured
p1 = [0.847518,0.879075,0.885870,0.889255,0.890911,0.891445,0.892520,0.893185,0.893414,0.893574,0.895203]

# probability of e+e- escaping
p2 = [0.999050,0.998559,0.997982,0.997331,0.996514,0.996042,0.995251,0.994628,0.993863,0.993076,0.985916]

# p1 x p2
p1p2 = [0.846713,0.877808,0.884082,0.886882,0.887805,0.887917,0.888281,0.888386,0.887931,0.887387,0.882595]

# plot data
plt.scatter(t, p1, s=1, color='blue', label='p1 = proton being captured')
plt.scatter(t, p2, s=1, color='red', label='p2 = e+e- escaping target')
plt.scatter(t, p1p2, s=1, color = 'black', label='p1p2')
plt.xlabel('target thickness (microns)', ha='right', x=1.0)
plt.ylabel('probability', ha='right', y=1.0)
plt.legend(frameon=True, loc='best', fontsize='x-small')
plt.tight_layout()
plt.show()
