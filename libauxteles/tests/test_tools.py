'''
Created on 6 nov. 2012

@author: colley
'''


import unittest
from tools import *

    

class Test(unittest.TestCase):


    def test_01(self):
        ret = readTextFileColumn('file_tools_01.txt',';')       
        print ret
        
        
    def test_02(self):
        # mon comm
        ret = readTextFileColumn('/home/colley/temp/lsst/modtran/TemplateT04.01_1.txt')        
        print ret


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()