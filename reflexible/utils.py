# Utilities live here

class Structure(dict):
    """Basically a dictionary whose keys are attributes too."""

    def __getattr__(self, attr):
        return self[attr]

    def __setattr__(self, attr, value):
        self[attr] = value

    def set_with_dict(self, D):
        """ set attributes with a dict """
        for k in D.keys():
            self.__setattr__(k, D[k])
