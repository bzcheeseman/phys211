__author__ = "Aman LaChapelle"

import numpy as np

def est_error(data, err_type=None):
    try:
        float(err_type)

        return err_type * np.ones_like(data)

    except ValueError:
        if err_type == "Counting":
            return np.sqrt(data)
        else:
            print "Not within function parameters"
            return None
