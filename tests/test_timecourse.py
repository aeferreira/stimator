import pytest

from six import StringIO
from numpy import isnan
from numpy.testing import assert_array_equal

from stimator import Solution
from stimator.modelparser import read_model
from stimator.timecourse import StimatorTCError

def assert_almost_equal(x, y):
    if abs(x-y) < 0.0001:
        return True
    return False

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
    sol.orderByNames("z y".split())
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
    sol.orderByNames(["z"])
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
    sol.orderByNames("x bof z".split())
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

if __name__ == '__main__':
    pytest.main()