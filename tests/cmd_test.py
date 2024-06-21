#! /usr/bin/python
import os
import sys
import unittest
import json
#from distutils.spawn import find_executable
import subprocess


class Util_test(unittest.TestCase):
    maxDiff = None
    def test_example(self):
        self.assertEqual(1,1) 

    def test_os_popen_example(self):
        cmd = 'echo example'
        with os.popen(cmd) as process:
            result = process.read()
        self.assertEqual(result, 'example\n')

    def test_help(self):        
        cmd = 'python src/cluster.py -h'
        with os.popen(cmd) as process:
            result = process.read()
        expected_output = (
            " _______________________________________________________________________________________________________________________ \n"
            "|***********************************************************************************************************************|\n"
            "| * List all available commands.                                                                                        |\n"
            "| python3 src/cluster.py -h                                                                                             |\n"
            "|_______________________________________________________________________________________________________________________|\n"
            "| * Make prediction with all input parameters in JSON file                                                              |\n"
            "| python3 src/cluster.py -j [input_json_file]                                                                           |\n"
            "| Example: python3 src/cluster.py -j examples/cluster.json                                                              |\n"
            "|_______________________________________________________________________________________________________________________|\n"
            "| * Make prediction with all input parameters in JSON file and write result to a JSON file                              |\n"
            "| python3 src/cluster.py -j [input_json_file] -o [output_prefix] -f json                                                |\n"
            "| Example: python3 src/cluster.py -j examples/cluster.json -o output -f json                                            |\n"
            "|_______________________________________________________________________________________________________________________|\n"
        )
        self.assertEqual(result, expected_output)
        # if not options are given, help info will be returned as well
        cmd = 'python3 src/cluster.py'
        with os.popen(cmd) as process:
            result = process.read()
        self.assertEqual(result, expected_output)



class TestCommonPrediction(unittest.TestCase):
    maxDiff = None

    def test_cluster_break(self):

        cmd="python src/cluster.py  -j test_data/cluster_break/input.json -o test_data/cluster_break/output -f json"
        with os.popen(cmd) as process:
            printed_result = process.read()
        with open('test_data/cluster_break/output.json','r') as f:
            result = json.load(f)
        with open('test_data/cluster_break/expected_output.json','r') as f:
            expected_output = json.load(f)

        self.assertEqual(result, expected_output)


def test():    
    unittest.main()

def main():
    print( 'run Cluster Analysis tests...')
    test()

if __name__ == '__main__':   
    main()

