from numpy.random import default_rng

n = 3_000
rng = default_rng()
random_i = rng.choice(mudata.obs_names, size=n, replace=False)

sample = mudata[random_i].copy()

# The same works for anndata (scanpy)
