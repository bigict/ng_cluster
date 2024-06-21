'''The validation should find all the errors and part of the warnings'''
from logging import getLogger

logger = getLogger(__name__)

def cluster_validate(data):
    errors = []
    warnings = []
    input_sequences = data.get('input_sequence_text', '')
    input_sequence_format = data.get('input_sequence_format', 'auto')
    if not input_sequences:
        errors.append('no input sequences')

    threshold = data.get('threshold', '0.7')
    threshold = float(threshold)
    if threshold < 0 or threshold > 1:
        errors.append('threshold "%s" is out of range 0 to 1' % threshold)
    return dict(errors=errors, warnings=warnings)