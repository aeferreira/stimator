import pprint

def write2file(filename, astring):
    f = open(filename, 'w')
    f.write(astring)
    f.close()

def s2HMS(seconds):
    m, s = divmod(seconds, 60.0)
    h, m = divmod(m, 60.0)
    if h == 0:
        return "%02dm %06.3fs" % (m, s)
    return "%dh %02dm %06.3fs" % (h, m, s)

#--------- dict dot-style lookup

# CREDIT: code by David Moss, Recipe 576586 from ActiveState Python recipes.
# http://code.activestate.com/recipes/576586/
class DictDotLookup(object):
    """
    Creates objects that behave much like a dictionaries (read-only), but allow nested
    key access using object '.' (dot) lookups.
    """
    def __init__(self, d):
        for k in d:
            if isinstance(d[k], dict):
                self.__dict__[k] = DictDotLookup(d[k])
            elif isinstance(d[k], (list, tuple)):
                l = []
                for v in d[k]:
                    if isinstance(v, dict):
                        l.append(DictDotLookup(v))
                    else:
                        l.append(v)
                self.__dict__[k] = l
            else:
                self.__dict__[k] = d[k]

    def __getitem__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]

    def __iter__(self):
        return iter(self.__dict__.keys())

    def __repr__(self):
        return pprint.pformat(self.__dict__)

if __name__ == '__main__':
    cfg_data = eval("""{
        'foo' : {
            'bar' : {
                'tdata' : (
                    {'baz' : 1 },
                    {'baz' : 2 },
                    {'baz' : 3 },
                ),
            },
        },
        'quux' : False,
    }""")

    cfg = DictDotLookup(cfg_data)

    #   Standard nested dictionary lookup.
    print 'normal lookup :', cfg['foo']['bar']['tdata'][0]['baz']

    #   Dot-style nested lookup.
    print 'dot lookup    :', cfg.foo.bar.tdata[0].baz
