# General utilities live here

from datetime import datetime

from matplotlib.dates import date2num
import numpy as np


class Structure(dict):
    """Basically a dictionary whose keys are attributes too."""

    def __getattr__(self, attr):
        return self[attr]

    def __setattr__(self, attr, value):
        self[attr] = value

    def set_with_dict(self, d):
        """ set attributes with a dict """
        for k in d.keys():
            self.__setattr__(k, d[k])


class CacheDict(dict):
    """A dictionary that prevents itself from growing too much."""

    def __init__(self, maxentries):
        self.maxentries = maxentries
        super(CacheDict, self).__init__(self)

    def __setitem__(self, key, value):
        # Protection against growing the cache too much
        if len(self) > self.maxentries:
            # Remove a 10% of (arbitrary) elements from the cache
            entries_to_remove = self.maxentries // 10
            entries_to_remove += 1  # to avoid removing no entries
            for k in self.keys()[:entries_to_remove]:
                super(CacheDict, self).__delitem__(k)
        super(CacheDict, self).__setitem__(key, value)


def closest(num, numlist):
    """ returns the index of the *closest* value in a list """
    # check if we're using datetimes
    if isinstance(num, datetime):
        num = date2num(num)
        assert isinstance(numlist[0], datetime), \
               "num is date, numlist must be a list of dates"
        numlist = date2num(numlist)

    return (np.abs(numlist - num)).argmin()
