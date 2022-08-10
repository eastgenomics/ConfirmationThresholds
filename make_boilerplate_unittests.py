from inspect import getmembers, isfunction
import FPTP 

out = 'test_boilerplate.py'
 
with open(out, 'w') as fh:
    output = 'from FPTP import *\nimport unittest\n\n\nclass TestModule(unittest.TestCase):\n'
    for i in getmembers(FPTP, isfunction):
            f = i[0]
            test = f'\tdef test_{f}(self):\n\t\tpass\n\n'
            output = f'{output}{test}'
    output = f'{output}\n\nif __name__ == \'__main__\':\n\tunittest.main()\n'
    fh.write(output)

