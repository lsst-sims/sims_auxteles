'''
Created on 17 mars 2014

@author: colley
'''

from modtranTools import *

S_FileTP7 = "/home/colley/temp/TempUS76.tp7"


def test_getKeywordTP7():
    getKeywordTP7()
    
    
def test_readTP7_multiCard():
    aRet, aCard = readTP7_multiCard(S_FileTP7, ['O3_TRANS', 'FREQ_CM-1'])
    print aRet.shape
    print aRet[0, :, 0]
    print aRet[10, :, 0]
    printCardLine(aCard[0])
    printCardLine(aCard[1])
    
    
def test_extractTemplate():
    extractTemplateFromTP7(S_FileTP7, '../../output/')



   
#test_getKeywordTP7()
test_readTP7_multiCard()
test_extractTemplate()

pl.show()