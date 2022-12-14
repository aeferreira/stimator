import pytest

from stimator.modelparser import read_model, StimatorParserError

model_text = """
#This is an example of a valid model:
title: A model to test parsing.
variables: X1 X2 X3

React1 : X2  + X3 -> X1, rate = Vmax1*X2*X3 / ((KmX3+X3)*(KmX2+X2))
leak : X3 -> 4.2 X3out, 10 ..
reaction React2 : X1 ->  2  OutVar,  \\
    step(t, 2.0, Vmax2*X1 / (Km2 + X1)) #reaction 2
kout_global = 3.14
export: OutVar ->, kout * OutVar, kout = sqrt(4.0)/2.0 * kout_global

in i1 = 20 - X2
-> i2 = i1 * 15
input i3 = i1 + i2

~ totX = X2 + X1
~ OutVarmult = mult * OutVar,      mult = (kout_global/export.kout) * 2
pi   = 3.1416
pi2  = 2*pi
pypi = pi**2  #this is pi square
KmX3 = sqrt(1e-2)
Vmax1 = 0.0001
find Vmax1 in [1e-9, 1e-3]
find   KmX3  in [1e-5, 1]
find KmX2 in [1e-5, pi/pi]

find Km2   in [1e-5, 1]
find Vmax2 in (1e-9, 1e-3)

find export.kout in (3,4)

@ 3.4 pi = 2*pi
x' = X3/2
#init  = state(X2 = 0.1, X3 = 0.63655, X1 = 0.0, x = 0)
init: X2 = 0.1, X3 = 0.63655, X1 = 0.0, x = 0

genomesize = 50 #should be enough
generations = 400
popsize = 20

timecourse my file.txt  # this is a timecourse filename
timecourse anotherfile.txt
#timecourse stillanotherfile.txt
tf: 10
!! X1 > X2 -> ~ ..

"""


def _get_error_loc(st_parse_error):
    sp = st_parse_error.value.physloc
    return (sp.nstartline, sp.nendline, sp.startlinepos, sp.endlinepos)


def _insert_line_and_string(textlines, i, line, replace=False):
    newlines = textlines[:]
    if not replace:
        newlines.insert(i, line)
    else:
        newlines[i] = line
    return '\n'.join(newlines)


@pytest.fixture
def textlines():
    return model_text.splitlines()


def test_valid_model():
    m = read_model(model_text)
    assert isinstance(m.info(), str)
    assert m.metadata['title'] == 'A model to test parsing.'
    assert len(m.varnames) == 6
    assert len(m.extvariables) == 0
    assert len(m.reactions) == 5
    assert len(m.reactions.React1.reagents) == 2
    assert len(m.reactions.d_x_dt.products) == 1
    assert len(m.transformations) == 2
    assert len(m.input_variables) == 3
    assert len(m.init) == 4
    assert len(m.parameters) == 11
    assert len(m.with_bounds) == 6
    assert m.parameters.KmX3 == 0.1
    assert m.parameters.Km2.bounds.upper == 1
    assert m.metadata['optSettings']['max_generations'] == 400
    assert len(m.metadata['timecourses']['filenames']) == 2
    assert m.metadata['timecourses']['filenames'][-1] == 'anotherfile.txt'
    assert len(m.metadata['timecourses']['defaultnames']) == 3
    assert m.metadata['timecourses']['defaultnames'][-1] == 'X3'


def test_name_undef(textlines):
    modelText = _insert_line_and_string(textlines,
                12, 'pypi = pip  #this is an error')
                #    0....v....1....v....2....v....3
    with pytest.raises(StimatorParserError) as spe:
        _ = read_model(modelText)

    sl, el, slp, elp = _get_error_loc(spe)
    assert 'NameError' in spe.value.value
    assert sl == 12 and el == 12 and slp == 7 and elp == 10


def test_name_undef_in_find(textlines):
    modelText = _insert_line_and_string(textlines,
                6, 'find pypi in [1e-5, 2 + kkk]  #this is an error')
                #   0....v....1....v....2....v....3....v....4....v....5
    with pytest.raises(StimatorParserError) as spe:
        _ = read_model(modelText)

    sl, el, slp, elp = _get_error_loc(spe)
    assert 'NameError' in spe.value.value
    assert sl == 6 and el == 6 and slp == 24 and elp == 27


def test_overflow(textlines):
    modelText = _insert_line_and_string(textlines,
                6, 'pypipip = pi*1e100**10000  #this is an overflow')
                #   0....v....1....v....2....v....3....v....4....v....5
    with pytest.raises(StimatorParserError) as spe:
        _ = read_model(modelText)

    sl, el, slp, elp = _get_error_loc(spe)
    assert 'OverflowError' in spe.value.value
    assert sl == 6 and el == 6 and slp == 10 and elp == 25


def test_repeated_decl(textlines):
    modelText = _insert_line_and_string(textlines,
                12, 'React1 : X2  + X3 -> X1, rate = 2 * X2 * X3')
                #    0....v....1....v....2....v....3....v....4....v
    with pytest.raises(StimatorParserError) as spe:
        _ = read_model(modelText)

    sl, el, slp, elp = _get_error_loc(spe)
    assert 'Repeated declaration' in spe.value.value
    assert sl == 12 and el == 12 and slp == 0 and elp == 43


def test_bad_name_in_rate(textlines):
    modelText = _insert_line_and_string(textlines,
                5, 'React1 : X2  + X3 -> X1, rate = X2*X3 / (KmX3+X32)', True)
                #   0....v....1....v....2....v....3....v....4....v....5
    with pytest.raises(StimatorParserError) as spe:
        _ = read_model(modelText)

    sl, el, slp, elp = _get_error_loc(spe)
    assert 'NameError' in spe.value.value
    assert sl == 5 and el == 5 and slp == 46 and elp == 49


def test_bad_rate(textlines):
    modelText = _insert_line_and_string(textlines,
                8, '    Vmax2*X1 / (Km2 + X1)) #reaction 2', True)
                #   0....v....1....v....2....v....3....v....4
    with pytest.raises(StimatorParserError) as spe:
        _ = read_model(modelText)

    sl, el, slp, elp = _get_error_loc(spe)
    assert 'Syntax Error' in spe.value.value
    assert sl == 8 and el == 8 and slp == 4 and elp == 27


def test_bad_syntax(textlines):
    modelText = _insert_line_and_string(textlines,
                6, 'OK!! not good...')
                #   0....v....1....v..
    with pytest.raises(StimatorParserError) as spe:
        _ = read_model(modelText)

    sl, el, slp, elp = _get_error_loc(spe)
    assert 'Invalid syntax' in spe.value.value
    assert sl == 6 and el == 6 and slp == 0 and elp == 16
