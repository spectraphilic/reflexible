"""
Module to keep functions only available inside UIO.
"""

import numpy as np

from .data_structures import Header


class UIOHeader(Header):

    def add_fires(self, **kwargs):
        """ uses the :mod:`emissions` module to read the MODIS hotspot data and
        add it to the header class as a 'fires' attribute.

        **This function is only available within UIO.**

        """

        from jfb.pflexpart import emissions as em
        self.fires = None
        for day in self.available_dates_dt:
            # day = day[:8]
            firedata = em.MODIS_hotspot(day)
            daily = firedata.daily
            # pdb.set_trace()
            if daily is None:
                continue
            else:
                if self.fires is None:
                    self.fires = daily
                else:
                    self.fires = np.hstack(
                        (self.fires, daily)).view(np.recarray)
