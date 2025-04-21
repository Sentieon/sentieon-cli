import argh


class CustomArghParser(argh.ArghParser):
    """
    Modify the argh parser to accept hyphens or underscores
    """

    def _parse_optional(self, arg_string):
        # normalize the string - replace hypens with underscores
        if (
            arg_string
            and len(arg_string) > 2
            and arg_string[0] in self.prefix_chars
            and arg_string[1] in self.prefix_chars
        ):
            arg_string = "--" + arg_string[2:].replace("-", "_")
        return super()._parse_optional(arg_string)
