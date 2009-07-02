import sys
import os.path

#append parent directory to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir))

import nose
args = ['nose', '-v']

result = nose.run(argv=args)
