'''
Home page form
'''
import os
from django import forms
from django.conf import settings
from django.forms.fields import EMPTY_VALUES
import tempfile
import re
import logging
from Bio import SeqIO

# from common.util import * #@UnusedWildImport
from common.setupinfo_parser import *  # @UnusedWildImport
from django.conf import settings

# ***TODO: MEDIA_ROOT is not an appropriate place for temporary data.
MEDIA_ROOT = settings.MEDIA_ROOT
tmpdir = os.path.join(MEDIA_ROOT, 'tmp')
cluster2_path = os.path.join(tmpdir, 'cluster2')
if not os.path.exists(cluster2_path):
    os.makedirs(cluster2_path)


class DynamicChoiceField(forms.ChoiceField):
    def clean(self, value):
        if value in EMPTY_VALUES:
            return None
        return value


class InputError(Exception):
    """Exception raised for errors in the input."""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value


class UserInputForm(forms.Form):
    sequence_text = forms.CharField(widget=forms.Textarea(
        attrs={'rows': 10, 'cols': 70, 'onKeyDown': 'countChar(document.formular.sequence_text,10485760)',
               'onKeyUp': 'countChar(document.formular.sequence_text,10485760)'}), required=False)
    sequence_file = forms.Field(widget=forms.FileInput, required=False)

    def clean(self):

        cleaned_data = super(UserInputForm, self).clean()
        textarea = cleaned_data['sequence_text']
        sequence_file = cleaned_data.get('sequence_file', None)

        # get input_sequences
        if sequence_file:
            input_sequences = sequence_file.read().decode()
        else:
            input_sequences = textarea

        # the max limit of the input peptides are 3000:
        peptides_num = len([line for line in input_sequences.splitlines() if not line.startswith('>')])
        if peptides_num > 3000:
            err = "please limit the number of peptides to 3,000 or less."
            raise forms.ValidationError(err)
        if peptides_num < 2:
            err = "Please provide input sequences. Minimum two peptides are required to perform clustering."
            raise forms.ValidationError(err)

            # create tmpdir/jobid
        cluster2_tmpdir = tempfile.mkdtemp(dir=cluster2_path)
        jobid = cluster2_tmpdir.split('/')[-1]
        tmpfile = open(cluster2_tmpdir + "/raw.fa", "wt")


        if input_sequences.startswith('>'):
            tmpfile.write(input_sequences)
        else:
            counting = 1
            for each in input_sequences.splitlines():
                tmpfile.write('>seq' + str(counting) + '\n' + str(each).upper() + '\n');
                counting += 1
        tmpfile.close()

        # check invalid amino acid
        counting = 0
        for ids, entry in enumerate(SeqIO.parse(cluster2_tmpdir + '/raw.fa', 'fasta')):
            for letter in str(entry.seq).upper():
                index = "ACDEFGHIKLMNPQRSTVWY".find(letter)
                if index == -1:
                    err = str(entry.seq) + " Sequence contains invalid amino acid."
                    raise forms.ValidationError(err)
            counting += 1
        if counting < 2:
            logging.info('err: less than two peptides input')
            err = "Please provide input sequences. Minimum two peptides are required to perform clustering."
            raise forms.ValidationError(err)

        logging.debug('cleaned_data: %s' % cleaned_data)
        cleaned_data['page2_data'] = {'jobid': jobid, }
        logging.debug('cleaned_data: %s' % cleaned_data)
        """
        cleaned_data = super(UserInputForm, self).clean()
        sequence_text = cleaned_data['sequence_text']
        sequence_file = cleaned_data['sequence_file']
        
        # raise error if no input sequence is entered
        if not sequence_text and not sequence_file:
            raise forms.ValidationError("You must enter a sequence(s).")
        
        sequence_list = []
        if '>' in sequence_text:
            input_sequences = sequence_text.split(">")
            for i in input_sequences[1:]:
                if(len(i) > 0):
                    end_of_name = i.find("\n")
                seq = i[end_of_name:].split()
                sequence_list.append(''.join(seq))
        else: sequence_list = sequence_text.split('\n')
        
        for sequence in sequence_list:
            for amino_acid in sequence.strip():
                if not amino_acid.upper() in "ACDEFGHIKLMNPQRSTVWY":
                    raise forms.ValidationError("Sequence: '%s' contains an invalid character: '%c' at position %d." %(sequence, amino_acid, sequence.find(amino_acid)))
        if len(sequence_text) > settings.MAX_UPLOAD_SIZE:
            raise forms.ValidationError("Please keep file size under %s. Current file size: %s" %(settings.MAX_UPLOAD_SIZE, len(sequence_text)))
        
        #if minimum_length > maximum_length:
        #    raise forms.ValidationError("Minimum peptide length selected %s. should not be higher than maximum peptide length %s" %(minimum_length, maximum_length))

        """

        return cleaned_data

    def process(self, request, interface='web'):

        if request.method == 'POST':
            textarea = request.POST.get('sequence_text', '')
            upload = ''

            if interface == 'api':
                upload = request.POST.get('sequence_file')
                upload = open(upload).read()
                input_source = 'local file'
            else:
                if request.FILES.has_key('sequence_file') or textarea is None:
                    input_source = 'local file'
                    upload = request.FILES['sequence_file'].read()
                else:
                    input_source = 'textarea'

            threshold = request.POST.get('threshold', '0.7')
            input_sequences = upload + textarea

            max_epitope_sequences = 50000

            if len(input_sequences) > 0:
                pass

            session = []
            session.append({'cluster_sequence': input_sequences.encode('ascii')})
            session.append(interface)
            request.session['session_data'] = session


class Step2Form(forms.Form):
    jobid = forms.CharField(required=False)
    threshold = forms.CharField(required=False)
    minimum_length = forms.CharField(required=False)
    maximum_length = forms.CharField(required=False)

    def clean(self):
        cleaned_data = super(Step2Form, self).clean()
        logging.debug('cleaned_data: %s' % cleaned_data)
        jobid = cleaned_data['jobid']
        threshold = float(cleaned_data.get('threshold', '0.7'))
        minimum_len = int(cleaned_data.get('minimum_length', '0'))
        maximum_len = int(cleaned_data.get('maximum_length', '1000'))

        cluster2_tmpdir = os.path.join(cluster2_path, jobid)

        # check Minimum two peptides in selected length range
        # SeqIO.parse(cluster2_tmpdir + '/raw.fa','fasta') is a generator
        if sum(1 for s in SeqIO.parse(cluster2_tmpdir + '/raw.fa', 'fasta') if
               len(s) >= minimum_len and len(s) <= maximum_len) < 2:
            logging.info('err: less than two peptides in selected length range')
            err = "Please provide proper input sequences and length range. Minimum two peptides in selected length range are required to perform clustering."
            raise forms.ValidationError(err)

        logging.debug('cleaned_data: %s' % cleaned_data)

        cleaned_data['page3_data'] = {'threshold': threshold, 'minimum_length': minimum_len,
                                      'maximum_length': maximum_len, 'jobid': jobid}
        logging.debug('cleaned_data: %s' % cleaned_data)

        return cleaned_data
