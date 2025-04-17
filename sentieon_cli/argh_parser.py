import sys

import argh


class CustomArghParser(argh.ArghParser):
    """
    Modify the argh parser to accept hyphens or underscores
    """

    def _parse_optional(self, arg_string):
        # if it's an empty string, it was meant to be a positional
        if not arg_string:
            return None

        # if it doesn't start with a prefix, it was meant to be positional
        if not arg_string[0] in self.prefix_chars:
            return None

        # if the option string is present in the parser, return the action
        if arg_string in self._option_string_actions:
            action = self._option_string_actions[arg_string]
            if sys.version_info.minor == 12:
                return [(action, arg_string, None, None)]
            elif sys.version_info.minor == 11:
                return action, arg_string, None, None
            elif sys.version_info.minor >= 8:
                return action, arg_string, None

        # if a long argument, convert hypens to underscores and match
        if len(arg_string) > 2 and arg_string[1] in self.prefix_chars:
            arg_string_normalized = "--" + arg_string[2:].replace("-", "_")
            if arg_string_normalized in self._option_string_actions:
                action = self._option_string_actions[arg_string_normalized]
                if sys.version_info.minor == 12:
                    return [(action, arg_string, None, None)]
                elif sys.version_info.minor == 11:
                    return action, arg_string, None, None
                elif sys.version_info.minor >= 8:
                    return action, arg_string, None

        return super()._parse_optional(arg_string)
