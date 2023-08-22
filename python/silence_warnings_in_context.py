import warnings

warnings.simplefilter("ignore") # Global

with warnings.catch_warnings(): # Only within the context
    warnings.simplefilter("ignore")
    import muon as mu # example
