#!/usr/bin/env python3

from .valueparsers import ValueParsers


class ParserDict(dict):

    def __init__(self):

        if(input_parser is None):
            val_parser_list = self.get_method_names(ValueParsers)
        else:
            val_parser_list = self.get_method_names(input_parser)

        for parser in val_parser_list:
            if(parser.startswith("parse_")):
                parse_key = parser[6:]
                self[parse_key] = getattr(ValueParsers, parser)
            else:
                raise SyntaxError(("A function in the ValueParsers class did "
                                   "not start with 'parse_'. Function is "
                                   "named: {}".format(parser)))

    @staticmethod
    def get_method_names(cls):
        return [func for func in dir(cls) if(callable(getattr(cls, func))
                                             and not func.startswith("__"))]
