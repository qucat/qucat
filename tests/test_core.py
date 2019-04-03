import unittest
import Qcircuits.core as core

class TestTesting(unittest.TestCase):
    def test_test(self):
        self.assertEqual(0,0)
        
class TestNetworkAnalysis(unittest.TestCase):

    '''
    Trivial cases
    '''
    def test_transfer_left_right_port_identical(self):
        net = core.Network([])
        self.assertEqual(net.transfer(0,1,0,1),1)

    def test_transfer_left_right_port_indentical_inverted(self):
        net = core.Network([])
        self.assertEqual(net.transfer(0,1,1,0),-1)

    '''
    Voltage divider
    '''
    def test_transfer_voltage_divider(self):
        