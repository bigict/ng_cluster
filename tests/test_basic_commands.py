import os
import re
import json
import subprocess
import unittest
from pathlib import Path


class TestCluster(unittest.TestCase):
    HOME_DIR = '/src/cluster'
    CWD = os.getcwd()

    def get_parent_dir(self):
        pdir = Path(self.CWD).parent.absolute()
        if pdir != self.HOME_DIR: pdir = '.'
        return pdir
    
    def standardize_format(self, text):
        ''' Take the input and replace all spaces with tabs, and remove empty strings '''
        text = text.replace(' ', '\t')
        final_text = re.sub(r"[\t]+", "\t", text)
        tmp_list = final_text.split('\n')
        tmp_list = list(filter(None, tmp_list))
        return ''.join(tmp_list)

    def test_basic_command(self):
        pdir = self.get_parent_dir()
        cmd = 'python %s/src/cluster.py -j %s/examples/cluster.json' %(pdir, pdir)
        formatted_cmd = cmd.split(' ')
        result = subprocess.run(formatted_cmd, capture_output=True)
        result = result.stdout.decode()

        result = self.standardize_format(result)

        eresult = '''1
cluster_number  peptide_number  alignment       position        sequence_number peptide cluster_consensus
1.1     Singleton       LEQIHVLENSLVL   -       Mus Pep1        LEQIHVLENSLVL   LEQIHVLENSLVL
2.1     Singleton       FVEHIHVLENSLAFK -       Mus Pep2        FVEHIHVLENSLAFK FVEHIHVLENSLAFK
3.1     Singleton       GLYGREPDLSSDIKERFA      -       Mus Pep3        GLYGREPDLSSDIKERFA      GLYGREPDLSSDIKERFA
4.1     Singleton       EWFSILLASDKREKI -       Mus Pep4        EWFSILLASDKREKI EWFSILLASDKREKI'''

        eresult = self.standardize_format(eresult)

        self.assertEqual(result, eresult)

    
    def test_job_descriptionn_and_its_output(self):
        '''
        This test case is specifically for the following command :
        > src/cluster.py -j examples/cluster.json --split --split-dir=examples/job/parameter_units

        This is the third example that is shown in the README.

        The above command will produce job description output path defaulted to '/src/cluster/examples/job/job_descriptions.json'.
        In the 'job_description.json' file, it will show two shell commands.

        This unittest will run both shell commands, and make sure their outputs matches to the expected output files located
        in the 'tests' directory.
        '''
        pdir = self.get_parent_dir()
        cmd = 'python %s/src/cluster.py -j %s/examples/cluster.json --split --split-dir=examples/job/parameter_units' %(pdir, pdir)
        formatted_cmd = cmd.split(' ')
        descr_result = subprocess.run(formatted_cmd, capture_output=True)
        descr_result = descr_result.stdout.decode()
        descr_result = list(filter(None, descr_result.split('\n')))
        job_description_path = descr_result[-1].split(':')[1].strip()
        with open(job_description_path, 'r') as f:
            content = json.load(f)
            first_content = content[0]
            second_content = content[1]


        ''' Testing the first job ID '''
        job_cmd = first_content['shell_cmd']
        
        # This particular example has only one output path
        job_result_file = first_content['expected_outputs'][0]        

        # Run the job command
        job_cmd = job_cmd.split(' ')
        subprocess.run(job_cmd)

        # Check expected output
        with open(job_result_file, 'r') as f:
            resulting_output = f.readlines()
            resulting_output = ''.join(resulting_output)
            resulting_output = resulting_output.replace(' ', '')
        

        expected_job_result_file = '%s/tests/expected_0.json' %(pdir)
        with open(expected_job_result_file, 'r') as f:
            expected_resulting_output = f.readlines()            
            expected_resulting_output = ''.join(expected_resulting_output)
            expected_resulting_output = expected_resulting_output.replace(' ', '')

        self.assertEqual(resulting_output, expected_resulting_output)


        ''' Testing the second job ID - aggregation '''
        job_cmd = second_content['shell_cmd']
        
        # This particular example has only one output path
        job_result_file = second_content['expected_outputs'][0]


        # Run the job command
        job_cmd = job_cmd.split(' ')
        subprocess.run(job_cmd)

        # Check expected output
        with open(job_result_file, 'r') as f:
            resulting_output = f.readlines()
            resulting_output = ''.join(resulting_output)
            resulting_output = resulting_output.replace(' ', '')
        

        expected_job_result_file = '%s/tests/expected_aggregated_result.json' %(pdir)
        with open(expected_job_result_file, 'r') as f:
            expected_resulting_output = f.readlines()            
            expected_resulting_output = ''.join(expected_resulting_output)
            expected_resulting_output = expected_resulting_output.replace(' ', '')

        self.assertEqual(resulting_output, expected_resulting_output)

    
if __name__=='__main__':
    unittest.main(warnings='ignore')