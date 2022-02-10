"""
A LL(1) parser to parse GTF-style attributes.
"""
from typing import Union, Dict, List, Tuple

# gene_id "STRG.23"; transcript_id "STRG.23.1"; reference_id "XM_017000355.2"; ref_gene_id "MIB2"; ref_gene_name "MIB2"; cov "1.073021"; FPKM "196.106232"; TPM "586.075928";

from libc.string cimport strlen, memset, strcmp, strcpy

parse_result_t = Dict[str, Union[str, float, int, bool, None]]

def parse(attributes:bytes) -> parse_result_t:
    retd = {}
    tokenizer = Tokenizer(attributes)
    tokenizer.run()
    current_keyname = ""

    for i in range(tokenizer.processed_number_of_tokens):
        token = tokenizer.get_token(i)
        print(token)
        if token[0] == _token_type.KEY:
            current_keyname = token[1]
        elif token[0] == _token_type.BOOL:
            retd[current_keyname] = _is_true(token[1])
        elif token[0] == _token_type.FLOAT:
            retd[current_keyname] = float(token[1])
        elif  token[0] == _token_type.INTEGER:
            retd[current_keyname] = int(token[1])
        elif token[0] == _token_type.NONE:
            retd[current_keyname] = None
        elif token[0] == _token_type.STRING:
            retd[current_keyname] = str(token[1], encoding='UTF-8')
    return retd

cdef struct _token:
    int TYPE
    char[256] VALUE
ctypedef _token token_t

cdef enum _token_state_type:
    WAITING_FOR_KEY = 0
    EXTENDING_KEY = 1
    WAITING_FOR_VALUE = 2
    EXTENDING_VALUE = 3

cdef enum _token_type:
    KEY = 0
    INTEGER = 1
    FLOAT = 2
    BOOL = 3
    STRING = 4
    NONE = 5

cdef bint _is_blank(char getv):
    return getv == '\t' or getv == '\n' or getv == ' '

cdef bint _is_quote(char getv):
    return getv == '\'' or getv == '\"'

cdef bint _is_true(char * token):
    return strcmp(token, b'TRUE') == 0 or \
           strcmp(token, b'true') == 0 or \
           strcmp(token, b'True') == 0 or \
           strcmp(token, b'YES') == 0 or \
           strcmp(token, b'yes') == 0 or \
           strcmp(token, b'Yes') == 0

cdef bint _is_false(char * token):
    return strcmp(token, b'FALSE') == 0 or \
           strcmp(token, b'false') == 0 or \
           strcmp(token, b'False') == 0 or \
           strcmp(token, b'NO') == 0 or \
           strcmp(token, b'no') == 0 or \
           strcmp(token, b'No') == 0



cdef class Tokenizer:
    cdef int attributes_len
    """Length of input attributes bytes"""

    cdef public token_t tokens[256]
    """Generated tokens"""

    cdef int current_ptr
    """Pointer of insertion locus of current attributes bytes"""

    cdef char* attributes
    """Input bytes of attributes"""

    cdef char current_token[256]
    """Value of the token being processed"""

    cdef int current_token_ptr
    """Pointer of insertion locus of current token"""

    cdef public int processed_number_of_tokens
    """Number of tokens processed"""

    def get_token(self, i: int) -> Tuple[int, bytes]:
        return (self.tokens[i].TYPE, self.tokens[i].VALUE)

    def __cinit__(self, char* attributes):
        self.attributes=attributes
        self.attributes_len = strlen(attributes)

    cdef bint peekable(self):
        """Whether have another char to get"""
        return self.current_ptr + 1 < self.attributes_len

    cdef char peek(self):
        """Get next char"""
        return self.attributes[self.current_ptr + 1]

    cdef char get(self):
        """Get char at current locus of attributes string"""
        return self.attributes[self.current_ptr]

    cdef void record_key(self):
        """Record key in self.tokens"""
        self.record_value(_token_type.KEY)

    cdef void record_value(self, int current_token_type):
        """Record value in self.tokens"""
        self.processed_number_of_tokens += 1
        self.tokens[self.processed_number_of_tokens].TYPE = current_token_type
        # self.tokens[self.processed_number_of_tokens].VALUE =
        strcpy(self.tokens[self.processed_number_of_tokens].VALUE, self.current_token)
        # print(current_token_type, self.current_token)

    cdef void init_new_token(self):
        memset(self.current_token, 0, 256)
        self.current_token_ptr = 0

    cdef void append_token(self, char getv):
        self.current_token[self.current_token_ptr] = getv
        self.current_token_ptr += 1


    cdef int update_token_type(self, char getv, int current_token_type):
        cdef int new_token_type = 0
        if current_token_type != _token_type.STRING:
            if 48 <= getv <= 57:  # Is numeric
                new_token_type = current_token_type
            elif getv == '.':
                if current_token_type == _token_type.FLOAT:
                    new_token_type =  _token_type.STRING
                elif current_token_type == _token_type.INTEGER:
                    new_token_type =  _token_type.FLOAT
            else:
                new_token_type =  _token_type.STRING
        else:
            new_token_type =  current_token_type
        # print(chr(getv), current_token_type, new_token_type)
        return new_token_type

    cdef int finalize_token_type(self, char* token, int current_token_type):
        if current_token_type == _token_type.FLOAT or current_token_type == _token_type.STRING:
            if strlen(token) == 0 or strlen(token) == 1 and token[0] == '.':
                return _token_type.NONE
            if _is_true(token) or _is_false(token):
                return _token_type.BOOL
        return current_token_type

    cdef void run(self) except *:
        cdef char getv
        cdef int current_status = _token_state_type.WAITING_FOR_KEY
        self.current_ptr = 0
        self.processed_number_of_tokens= -1
        cdef int current_token_type = _token_type.INTEGER
        cdef int quotation_mark = 0
        while self.current_ptr < self.attributes_len:
            getv = self.get()
            # print(chr(getv), current_status)
            if current_status == _token_state_type.WAITING_FOR_KEY:
                if _is_blank(getv) or getv == ';':
                    pass
                elif _is_quote(getv):
                    quotation_mark = getv
                    current_status = _token_state_type.EXTENDING_KEY
                    self.init_new_token()
                else:
                    current_status = _token_state_type.EXTENDING_KEY
                    self.init_new_token()
                    self.append_token(getv)
            elif current_status == _token_state_type.EXTENDING_KEY:
                if quotation_mark != 0 and _is_quote(getv) and getv == quotation_mark:
                    if self.peekable():
                        if _is_blank(self.peek()):
                            current_status = _token_state_type.WAITING_FOR_VALUE
                            quotation_mark = 0
                            self.record_key()
                        else:
                            self.append_token(getv)
                    else:
                        raise ValueError(f"Waiting for value at {self.current_ptr}")
                elif _is_blank(getv):
                    current_status = _token_state_type.WAITING_FOR_VALUE
                    self.record_key()
                else:
                    self.append_token(getv)
            elif current_status == _token_state_type.WAITING_FOR_VALUE:
                if _is_blank(getv):
                    pass
                elif _is_quote(getv):
                    quotation_mark = getv
                    current_status = _token_state_type.EXTENDING_VALUE
                    self.init_new_token()
                    current_token_type = _token_type.INTEGER
                else:
                    current_status = _token_state_type.EXTENDING_VALUE
                    self.init_new_token()
                    current_token_type = _token_type.INTEGER
                    self.append_token(getv)
            elif current_status == _token_state_type.EXTENDING_VALUE:
                if quotation_mark != 0 and _is_quote(getv) and getv == quotation_mark:
                    if self.peekable():
                        if _is_blank(self.peek()) or self.peek() == ';':
                            current_status = _token_state_type.WAITING_FOR_KEY
                            quotation_mark = 0
                            self.record_value(current_token_type)
                        else:
                            current_token_type = self.update_token_type(getv, current_token_type)
                            self.append_token(getv)
                elif _is_blank(getv):
                    current_status = _token_state_type.WAITING_FOR_KEY
                    self.record_value(current_token_type)
                else:
                    current_token_type = self.update_token_type(getv, current_token_type)
                    self.append_token(getv)
            self.current_ptr += 1
        if current_status == _token_state_type.WAITING_FOR_KEY:
            pass
        elif current_status == _token_state_type.EXTENDING_KEY:
            if quotation_mark == 0:
                raise ValueError(f"Waiting for value at {self.current_ptr}")
            else:
                raise ValueError(f"Waiting for end quote at {self.current_ptr}")
        elif current_status == _token_state_type.WAITING_FOR_VALUE:
            raise ValueError(f"Waiting for value at {self.current_ptr}")
        elif current_status == _token_state_type.EXTENDING_VALUE:
            if quotation_mark == 0:
                current_status = _token_state_type.WAITING_FOR_KEY
                self.record_value(current_token_type)
            else:
                raise ValueError(f"Waiting for end quote at {self.current_ptr}")


