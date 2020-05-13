#!/usr/bin/env python3

import json
import os.path

from .parserdict import ParserDict
from .exceptions import CGECoreOutTypeError, CGECoreOutInputError


class Result(dict):

    BEONE_JSON_FILE = "beone.json"

    beone_json_path = os.path.join(os.path.dirname(__file__), BEONE_JSON_FILE)
    beone_defs = {}
    with open(beone_json_path, "r") as fh:
        beone_defs = json.load(fh)

    val_parsers = ParserDict()

    def __init__(self, result_type=None, **kwargs):
        type = self._get_type(result_type, **kwargs)
        self._set_type(type)
        self._parser = ResultParser(result_def=self.beone_defs[type])
        for d in self._parser.arrays:
            self[d] = []
        for d in self._parser.dicts:
            self[d] = {}

        self.add(**kwargs)

    def _set_type(self, type):
        if(type in self.beone_defs):
            self["type"] = type
        else:
            raise CGECoreOutTypeError(
                "Unknown result type given. Type given: {}. Type must be one "
                "of:\n{}".format(type, list(self.beone_defs.keys())))

    def _get_type(self, result_type=None, **kwargs):
        type = None
        if(result_type is not None):
            type = result_type
        if(kwargs):
            kw_type = kwargs.get("type", None)
            if(type is not None and kw_type is not None and type != kw_type):
                raise CGECoreOutTypeError(
                    "Type was given as argument to method call and as an "
                    "attribute in the given dictionary, but they did not "
                    "match. {} (method) != {} (dict)".format(type, kw_type))
            elif(kw_type is not None):
                type = kw_type
        if(type is None):
            raise CGECoreOutTypeError(
                "The class format requires a 'type' attribute. The given "
                "dictionary contained the following attributes: {}"
                .format(kwargs.keys()))
        return type

    def add(self, **kwargs):
        for key, val in kwargs.items():
            if(val is None):
                continue
            self[key] = val

    def add_class(self, cl, result_type=None, **kwargs):
        type = self._get_type(result_type, **kwargs)
        res = Result(result_type=type, **kwargs)
        if(cl in self._parser.arrays):
            self[cl].append(res)
        elif(cl in self._parser.dicts):
            self[cl][res["key"]] = res
        else:
            self[cl] = res

    def check_results(self, errors=None):
        self.errors = {}

        for key, val in self.items():
            if(key == "type"):
                continue
            self._check_result(key, val, self.errors)

        # errors is not None if called recursively
        if(errors is not None):
            errors[self["key"]] = self.errors
            return None
        # errors is None if it is the first/root call
        elif(errors is None and self._no_errors(self.errors)):
            return None
        else:
            raise CGECoreOutInputError(
                "Some input data did not pass validation, please consult the "
                "Dictionary of ERRORS:{}".format(self.errors),
                self.errors)

    def _check_result(self, key, val, errors, index=None):
        # Remember Result is a dict object and therefore this test should
        # be before the dict test.
        if(isinstance(val, Result)):
            val.check_results(errors)
        elif(isinstance(val, dict)):
            self._check_result_dict(key, val, errors)
        elif(isinstance(val, list)):
            self._check_result_list(key, val, errors)
        else:
            self._check_result_val(key, val, errors, index)

    def _no_errors(self, errors):
        no_errors = True

        for key, val in errors.items():

            if(isinstance(val, dict)):
                no_errors = self._no_errors(val)
                if(no_errors is False):
                    return False

            elif(val is not None):
                return False

        return no_errors

    def _check_result_val(self, key, val, errors, index=None):
        val_type = self._parser[key]

        if(val_type.endswith("*")):
            val_type = val_type[:-1]

        val_error = self.val_parsers[val_type](val)

        if(val_error):
            if(index is not None):
                val_error = "{}:{} ".format(index, val_error)
            errors[key] = val_error

    def _check_result_dict(self, result_key, result_dict, errors):
        errors[result_key] = {}
        for key, val in result_dict.items():
            self._check_result(key, val, errors[result_key])

    def _check_result_list(self, result_key, result_list, errors):
        errors[result_key] = {}
        for i, val in enumerate(result_list):
            self._check_result(result_key, val, errors[result_key], index=i)


class ResultParser(dict):
    """"""
    def __init__(self, result_def):
        self.classes = set()
        self.arrays = {}
        self.dicts = {}

        for key, val_def_str in result_def.items():
            val_def, *sub_def = val_def_str.split(" ")
            if(sub_def and val_def == "dict"):
                self.dicts[key] = sub_def
                self[key] = sub_def
            elif(sub_def and val_def == "array"):
                self.arrays[key] = sub_def
                self[key] = sub_def
            else:
                self[key] = val_def
