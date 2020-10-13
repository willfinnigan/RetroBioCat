from unidecode import unidecode

INVALID_NAME_CHARS = [".", "(", ")", "'", "/"]

def sanitise_string(text):
    for char in text:
        if char in INVALID_NAME_CHARS:
            text = text.replace(char, '')

    return unidecode(str(text))


def sanitise_sequence(seq_string):
    seq_string = seq_string.replace('\n', '')
    seq_string = seq_string.replace(' ', '')
    return seq_string


def sequence_check(sequence):
    amino_acids_list = ['X', 'V', 'G', 'F', 'E', 'N', 'P', 'Q', 'M', 'K', 'T', 'S', 'W', 'A', 'R', 'D', 'L', 'Y', 'H',
                        'I', 'C', '*']

    bad_chars = []
    for letter in sequence:
        if letter.upper() not in amino_acids_list:
            bad_chars.append(letter)

    return bad_chars