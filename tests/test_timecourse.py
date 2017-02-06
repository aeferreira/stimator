import pytest

import os

from six import StringIO
from numpy import isnan, array
from numpy.testing import assert_array_equal

from stimator import Solution, Solutions, read_tc
from stimator.modelparser import read_model

_THIS_DIR, _ = os.path.split(os.path.abspath(__file__))
#print(_THIS_DIR)
_UPPER, _ = os.path.split(_THIS_DIR)
#print(_UPPER)
_DATADIR = os.path.join(_UPPER, 'examples')
#print(_DATADIR)

def assert_almost_equal(x, y):
    if abs(x-y) < 0.0001:
        return True
    return False

def average(x, t):
    return array((t/2.0, (x[0]+x[-1])/2.0))

demodata = """
#this is demo data with a header
t x y z
0       0.95 0         0
0.1                  0.1

  0.2 skip 0.2 skip this
nothing really usefull here
- 0.3 0.3 this line should be skipped
#0.4 0.4
0.3 0.4 0.5 0.55
0.4 0.5 0.6 0.7
0.5 0.6 0.8 0.9
0.55 0.7 0.85 0.95
0.6  - 0.5 - -

"""
@pytest.fixture
def tc_1():
    return StringIO(demodata)

demodata_noheader = """
#this is demo data without a header
#t x y z
0       1 0         0
0.1                  0.1

  0.2 skip 0.2 skip this
nothing really usefull here
- 0.3 0.3 this line should be skipped
#0.4 0.4
0.5  - 0.5 - -
0.6 0.6 0.8 0.9

"""

demodata2 = """
#this is demo data with a header
t x y z
0       0.95 0         0
0.1                  0.09

  0.2 skip 0.2 skip this
nothing really usefull here
- 0.3 0.3 this line should be skipped
#0.4 0.4
0.3 0.45 0.55 0.58
0.4 0.5 0.65 0.75
0.5 0.65 0.85 0.98
0.55 0.7 0.9 1.45
0.6  - 0.4 - -
"""

def test_read_from(tc_1):
    sol = Solution().read_from(tc_1)
    assert sol.names == ['x', 'y', 'z']
    assert sol.t[0] == 0.0
    assert sol.t[-1] == 0.6
    assert len(sol.t) == 8
    assert sol.data.shape == (3, 8)
    assert sol.data[0, 0] == 0.95
    assert sol.data[0, 3] == 0.4
    assert isnan(sol.data[-1, -1])
    ptc = str(sol)
    lines = [line.strip() for line in ptc.split('\n') if len(line) > 0]
    assert len(lines) == 9
    assert lines[0] == 't x y z'

def test_read_str_orderByNames():
    sol = Solution().read_str(demodata)
    sol.order_by_names("z y".split())
    assert sol.names == ['z', 'y', 'x']
    assert sol.t[0] == 0.0
    assert sol.t[-1] == 0.6
    assert len(sol.t) == 8
    assert sol.data.shape == (3, 8)
    assert sol.data[0, 0] == 0
    assert sol.data[2, 3] == 0.4
    assert isnan(sol.data[-1, -1])

def test_read_str_orderByNames2():
    sol = Solution().read_str(demodata)
    sol.order_by_names(["z"])
    assert sol.names == ['z', 'x', 'y']
    assert sol.t[0] == 0.0
    assert sol.t[-1] == 0.6
    assert len(sol.t) == 8
    assert sol.data.shape == (3, 8)
    assert sol.data[0, 0] == 0
    assert sol.data[1, 3] == 0.4
    assert isnan(sol.data[1, -1])

def test_read_str_bad_order_by_names():
    sol = Solution().read_str(demodata)
    sol.order_by_names("x bof z".split())
    assert sol.names == ['x', 'z', 'y']
    assert sol.t[0] == 0.0
    assert sol.t[-1] == 0.6
    assert len(sol.t) == 8
    assert sol.data.shape == (3, 8)

def test_read_data_without_header():
    sol = Solution().read_str(demodata_noheader)
    assert sol.names == ['x1', 'x2', 'x3']
    assert sol.t[0] == 0.0
    assert sol.t[-1] == 0.6
    assert len(sol.t) == 5
    assert sol.data.shape == (3, 5)
    assert sol.data[0, 0] == 1
    assert sol.data[1, 3] == 0.5
    assert isnan(sol.data[2, 1])

def test_read_data_without_header_giving_names():
    names = ['v1', 'v2', 'v3', 'v4', 'v5']
    sol = Solution().read_str(demodata_noheader, names=names)
    assert sol.names == ['v1', 'v2', 'v3']
    assert sol.t[0] == 0.0
    assert sol.t[-1] == 0.6
    assert len(sol.t) == 5
    assert sol.data.shape == (3, 5)
    assert sol.data[0, 0] == 1
    assert sol.data[1, 3] == 0.5
    assert isnan(sol.data[2, 1])

def test_read_data_without_header_giving_names2():
    names = ['v1', 'v2']
    sol = Solution().read_str(demodata_noheader, names=names)
    assert sol.names == ['v1', 'v2', 'x3']
    assert sol.t[0] == 0.0
    assert sol.t[-1] == 0.6
    assert len(sol.t) == 5
    assert sol.data.shape == (3, 5)
    assert sol.data[0, 0] == 1
    assert sol.data[1, 3] == 0.5
    assert isnan(sol.data[2, 1])

def test_Solution_interface():
    sol = Solution().read_str(demodata)
    
    # len(sol) and sol.ntimes
    assert(len(sol)) == 3
    assert sol.ntimes == 8
    assert sol.names == ['x', 'y', 'z']
    
    # sol.t
    assert sol.t[0] == 0.0
    assert sol.t[-1] == 0.6
    assert len(sol.t) == 8
    
    # sol.data
    assert sol.data.shape == (3, 8)
    assert sol.data[0, 0] == 0.95
    assert sol.data[0, 3] == 0.4
    assert isnan(sol.data[-1, -1])
    # sol is an array, vectorial operators apply
    y = 2.0 * sol.data[:, -1]
    assert y[1] == 1
    
    # indexing
    assert_array_equal(sol[0], sol.data[0])
    assert_array_equal(sol['x'], sol.data[0])
    with pytest.raises(ValueError):
        kseries = sol['k']
        assert kseries[0] == 0.0
    
    # state_at(), returns dictionaries
    s02 = sol.state_at(0.2)
    assert isnan(s02['x'])
    assert isnan(s02['z'])
    assert assert_almost_equal(s02['y'], 0.2)
    s045 = sol.state_at(0.45) # linear interpolation
    assert assert_almost_equal(s045['x'], 0.55)
    assert assert_almost_equal(s045['y'], 0.7)
    assert assert_almost_equal(s045['z'], 0.8)

    # init() and last(), returns dictionaries
    sinit = sol.init
    assert assert_almost_equal(sinit['x'], 0.95)
    assert assert_almost_equal(sinit['y'], 0.0)
    assert assert_almost_equal(sinit['z'], 0.0)
    slast = sol.last
    assert isnan(slast['x'])
    assert isnan(slast['z'])
    assert assert_almost_equal(slast['y'], 0.5)
    
    # iteration
    for series, sdata in zip(sol, sol.data):
        assert_array_equal(series, sdata)
    
    # writing to file
    outfile = StringIO()
    sol.write_to(outfile)
    outfile.seek(0)
    sol.read_from(outfile)
    assert(len(sol)) == 3
    assert sol.ntimes == 8
    assert sol.names == ['x', 'y', 'z']
    assert sol.t[0] == 0.0
    assert sol.t[-1] == 0.6
    assert len(sol.t) == 8
    assert sol.data.shape == (3, 8)
    assert sol.data[0, 0] == 0.95
    assert sol.data[0, 3] == 0.4
    assert isnan(sol.data[-1, -1])

def test_Solution_transformation():
    sol = Solution(title='original time course').read_str(demodata)
    # before transformation
    assert sol.names == ['x', 'y', 'z']
    assert sol.t[0] == 0.0
    assert sol.t[-1] == 0.6
    assert len(sol.t) == 8
    assert sol.data.shape == (3, 8)
    assert sol.data[0, 0] == 0.95
    assert sol.data[0, 3] == 0.4
    assert isnan(sol.data[-1, -1])
    sol = sol.transform(average,
                        newnames=['t/2', 'mid point'],
                        new_title='after transformation')

    # after transformation
    assert sol.names == ['t/2', 'mid point']
    assert sol.t[0] == 0.0
    assert sol.t[-1] == 0.6
    assert len(sol.t) == 8
    assert sol.data.shape == (2, 8)
    assert sol.data[0, 0] == sol.t[0] / 2.0
    assert sol.data[0, -1] == sol.t[-1] / 2.0
    assert sol.data[1, 0] == (0.95 + 0.0)/ 2.0
    assert isnan(sol.data[1, -1])

def test_Solution_clone_copy():
    sol = Solution(title='original time course').read_str(demodata)
    # original
    assert sol.names == ['x', 'y', 'z']
    assert sol.t[0] == 0.0
    assert sol.t[-1] == 0.6
    assert sol.ntimes == 8
    assert sol.data.shape == (3, 8)
    assert sol.init['x'] == 0.95
    assert sol.init['z'] == 0.0
    assert isnan(sol.last['z'])

    sol2 = sol.clone()
    sol3 = sol.copy()

    # cloned solutions
    assert sol2.names == ['x', 'y', 'z']
    assert sol2.t[0] == 0.0
    assert sol2.t[-1] == 0.6
    assert sol2.ntimes == 8
    assert sol2.data.shape == (3, 8)
    assert sol2.init['x'] == 0.95
    assert sol2.init['z'] == 0.0
    assert isnan(sol2.last['z'])

    assert sol3.names == ['x', 'y', 'z']
    assert sol3.t[0] == 0.0
    assert sol3.t[-1] == 0.6
    assert sol3.ntimes == 8
    assert sol3.data.shape == (3, 8)
    assert sol3.init['x'] == 0.95
    assert sol3.init['z'] == 0.0
    assert isnan(sol3.last['z'])

    sol4 = sol.copy('y')

    assert sol4.names == ['y']
    assert sol4.t[0] == 0.0
    assert sol4.t[-1] == 0.6
    assert sol4.ntimes == 8
    assert sol4.data.shape == (1, 8)
    assert sol4.init['y'] == 0.0
    assert sol4.last['y'] == 0.5


def test_Solutions_construction_and_iadd():
    sols = Solutions(title='all time courses')
    bcontext = True if sols else False
    assert not bcontext
    s = Solution(title='1st time course').read_str(demodata)
    sols += s
    s = Solution(title='2nd time course').read_str(demodata_noheader)
    sols += s
    bcontext = True if sols else False
    assert bcontext
    
    print_1st_line = 't x y z'
    ssols = str(sols)
    ssols = [line.strip() for line in ssols.split('\n')]
    assert print_1st_line == ssols[0]

def test_readTCs():
    tcs = read_tc(['TSH2b.txt', 'TSH2a.txt'], _DATADIR, verbose=False)
    assert len(tcs) == 2

    assert tcs[0].shape == (2, 347)
    assert tcs[0].names == ['SDLTSH', 'HTA']
    assert assert_almost_equal(tcs[0].init['SDLTSH'], 0.001246154)
    assert assert_almost_equal(tcs[0].init['HTA'], 0.2688)
    assert assert_almost_equal(tcs[0].last['SDLTSH'], 0.042815385)
    assert isnan(tcs[0].last['HTA'])
    assert tcs[0].shortname == 'TSH2b.txt'

    assert tcs[1].shape == (1, 244)
    assert tcs[1].names == ['x1']
    assert assert_almost_equal(tcs[1].init['x1'], 7.69231E-05)
    assert assert_almost_equal(tcs[1].last['x1'], 0.022615385)
    assert tcs[1].shortname == 'TSH2a.txt'

def test_readTCs_default_names():
    tcs = read_tc(['TSH2b.txt', 'TSH2a.txt'], _DATADIR,
                  names="SDLTSH HTA".split(),
                  verbose=False)
    assert len(tcs) == 2

    assert tcs[0].shape == (2, 347)
    assert tcs[0].names == ['SDLTSH', 'HTA']
    assert assert_almost_equal(tcs[0].init['SDLTSH'], 0.001246154)
    assert assert_almost_equal(tcs[0].init['HTA'], 0.2688)
    assert assert_almost_equal(tcs[0].last['SDLTSH'], 0.042815385)
    assert isnan(tcs[0].last['HTA'])
    assert tcs[0].shortname == 'TSH2b.txt'

    assert tcs[1].shape == (1, 244)
    assert tcs[1].names == ['SDLTSH']
    assert assert_almost_equal(tcs[1].init['SDLTSH'], 7.69231E-05)
    assert assert_almost_equal(tcs[1].last['SDLTSH'], 0.022615385)
    assert tcs[1].shortname == 'TSH2a.txt'
    assert tcs.get_common_full_vars() == ['SDLTSH']

def test_readTCs_and_change_order():
    tcs = read_tc(['TSH2b.txt', 'TSH2a.txt'], _DATADIR, verbose=False)
    tcs.order_by_names('HTA SDLTSH'.split())
    assert len(tcs) == 2

    assert tcs[0].shape == (2, 347)
    assert tcs[0].names == ['HTA', 'SDLTSH']
    assert assert_almost_equal(tcs[0].init['SDLTSH'], 0.001246154)
    assert assert_almost_equal(tcs[0].init['HTA'], 0.2688)
    assert assert_almost_equal(tcs[0].last['SDLTSH'], 0.042815385)
    assert isnan(tcs[0].last['HTA'])
    assert tcs[0].shortname == 'TSH2b.txt'

    assert tcs[1].shape == (1, 244)
    assert tcs[1].names == ['x1']
    assert assert_almost_equal(tcs[1].init['x1'], 7.69231E-05)
    assert assert_almost_equal(tcs[1].last['x1'], 0.022615385)
    assert tcs[1].shortname == 'TSH2a.txt'

def test_write_to():
    tcs = read_tc(['TSH2b.txt', 'TSH2a.txt'], _DATADIR, verbose=False)
    tcs.write_to(['TSH2b_2.txt', 'TSH2a_2.txt'], filedir='../', verbose=False)
    tcs = read_tc(['TSH2b_2.txt', 'TSH2a_2.txt'], '../', verbose=False)
    assert len(tcs) == 2

    assert tcs[0].shape == (2, 347)
    assert tcs[0].names == ['SDLTSH', 'HTA']
    assert assert_almost_equal(tcs[0].init['SDLTSH'], 0.001246154)
    assert assert_almost_equal(tcs[0].init['HTA'], 0.2688)
    assert assert_almost_equal(tcs[0].last['SDLTSH'], 0.042815385)
    assert isnan(tcs[0].last['HTA'])
    assert tcs[0].shortname == 'TSH2b_2.txt'

    assert tcs[1].shape == (1, 244)
    assert tcs[1].names == ['x1']
    assert assert_almost_equal(tcs[1].init['x1'], 7.69231E-05)
    assert assert_almost_equal(tcs[1].last['x1'], 0.022615385)
    assert tcs[1].shortname == 'TSH2a_2.txt'
    assert os.path.isfile('../TSH2b_2.txt')
    assert os.path.isfile('../TSH2a_2.txt')
    os.remove('../TSH2b_2.txt')
    os.remove('../TSH2a_2.txt')
    assert not os.path.isfile('../TSH2b_2.txt')
    assert not os.path.isfile('../TSH2a_2.txt')

def test_read_tc_declared_in_model():
    m = read_model("""
    v1: HTA -> SDLTSH, rate = 1 ..
    v2: SDLTSH -> ,    rate = 2 ..
    timecourse TSH2b.txt
    timecourse TSH2a.txt
    variables SDLTSH HTA
    """)
    
    tcs = read_tc(m, _DATADIR, verbose=False)
    
    assert len(tcs) == 2

    assert tcs[0].shape == (2, 347)
    assert tcs[0].names == ['SDLTSH', 'HTA']
    assert assert_almost_equal(tcs[0].init['SDLTSH'], 0.001246154)
    assert assert_almost_equal(tcs[0].init['HTA'], 0.2688)
    assert assert_almost_equal(tcs[0].last['SDLTSH'], 0.042815385)
    assert isnan(tcs[0].last['HTA'])
    assert tcs[0].shortname == 'TSH2b.txt'

    assert tcs[1].shape == (1, 244)
    assert tcs[1].names == ['SDLTSH']
    assert assert_almost_equal(tcs[1].init['SDLTSH'], 7.69231E-05)
    assert assert_almost_equal(tcs[1].last['SDLTSH'], 0.022615385)
    assert tcs[1].shortname == 'TSH2a.txt'
    
    tcs.order_by_modelvars(m)

    assert len(tcs) == 2

    assert tcs[0].shape == (2, 347)
    assert tcs[0].names == ['HTA', 'SDLTSH']
    assert assert_almost_equal(tcs[0].init['SDLTSH'], 0.001246154)
    assert assert_almost_equal(tcs[0].init['HTA'], 0.2688)
    assert assert_almost_equal(tcs[0].last['SDLTSH'], 0.042815385)
    assert isnan(tcs[0].last['HTA'])
    assert tcs[0].shortname == 'TSH2b.txt'

    assert tcs[1].shape == (1, 244)
    assert tcs[1].names == ['SDLTSH']
    assert assert_almost_equal(tcs[1].init['SDLTSH'], 7.69231E-05)
    assert assert_almost_equal(tcs[1].last['SDLTSH'], 0.022615385)
    assert tcs[1].shortname == 'TSH2a.txt'

if __name__ == '__main__':
    pytest.main()