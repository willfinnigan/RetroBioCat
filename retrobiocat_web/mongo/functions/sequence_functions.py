from unidecode import unidecode

INVALID_NAME_CHARS = [".", "(", ")", "'", "/"]

def sanitise_string(text):
    for char in text:
        if char in INVALID_NAME_CHARS:
            text = text.replace(char, '')

    return unidecode(str(text))
