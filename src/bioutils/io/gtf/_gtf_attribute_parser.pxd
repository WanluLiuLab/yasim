"""
A LL(1) parser to parse GTF-style attributes.
"""
from typing import Union, Dict

parse_result_t = Dict[str, Union[str, float, int, bool, None]]



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

cdef bint _is_blank(char getv)

cdef bint _is_quote(char getv)

cdef bint _is_true(char * token)

cdef bint _is_false(char * token)


cdef class GtfAttributeTokenizer:
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

    cdef bint peekable(self)

    cdef char peek(self)

    cdef char get(self)

    cdef void record_key(self)

    cdef void record_value(self, int current_token_type)

    cdef void init_new_token(self)

    cdef void append_token(self, char getv)


    cdef int update_token_type(self, char getv, int current_token_type)

    cdef int finalize_token_type(self, char* token, int current_token_type)

    cdef void run(self) except *
