import numpy as np
import scipy.optimize as sopt    # RLAP solver
import matplotlib.pyplot as plt  # visualizatiion
import seaborn as sns            # """
np.random.seed(1)

# Example data from
# https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html
# removed a row; will be shuffled to make it more interesting!
harvest = np.array([[0.8, 2.4, 2.5, 3.9, 0.0, 4.0, 0.0],
                    [2.4, 0.0, 4.0, 1.0, 2.7, 0.0, 0.0],
                    [1.1, 2.4, 0.8, 4.3, 1.9, 4.4, 0.0],
                    [0.6, 0.0, 0.3, 0.0, 3.1, 0.0, 0.0],
                    [0.7, 1.7, 0.6, 2.6, 2.2, 6.2, 0.0],
                    [1.3, 1.2, 0.0, 0.0, 0.0, 3.2, 5.1]],)
harvest = harvest[:, np.random.permutation(harvest.shape[1])]

# scipy: linear_sum_assignment -> able to take rectangular-problem!
# assumption: minimize -> cost-matrix to profit-matrix:
#                         remove original cost from maximum-costs
#                         Kuhn, Harold W.:
#                         "Variants of the Hungarian method for assignment problems."
max_cost = np.amax(harvest)
harvest_profit = max_cost - harvest

row_ind, col_ind = sopt.linear_sum_assignment(harvest_profit)
sol_map = np.zeros(harvest.shape, dtype=bool)
sol_map[row_ind, col_ind] = True

# Visualize
f, ax = plt.subplots(2, figsize=(9, 6))
sns.heatmap(harvest, annot=True, linewidths=.5, ax=ax[0], cbar=False,
            linecolor='black', cmap="YlGnBu")
sns.heatmap(harvest, annot=True, mask=~sol_map, linewidths=.5, ax=ax[1],
            linecolor='black', cbar=False, cmap="YlGnBu")
plt.tight_layout()
plt.show()