# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: test_pickle_with_tqdm.py -- Unit test of corresponding module.
#
#  VERSION HISTORY:
#  2021-09-23 0.1  : Purposed and added bu YU Zhejian
#
# ==============================================================================
"""
test_pickle_with_tqdm.py -- Unit test of corresponding module.
"""
import pickle
import random

import test_tetgs
from commonutils import pickle_with_tqdm

test_path = test_tetgs.initialize(__name__)


def test_unpickle_with_tqdm():
    random_arr = []
    for _ in range(1000):
        random_arr.append(random.random())
    pickle_fn = f'{test_path}/rd.pickle'
    with open(pickle_fn, 'wb') as writer:
        pickle.dump(random_arr, writer)
    unpickle_obj = pickle_with_tqdm.unpickle_with_tqdm(pickle_fn)
    assert unpickle_obj == random_arr
