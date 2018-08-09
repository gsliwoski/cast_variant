# Custom exceptions used while casting a variant

class ParseException(Exception):
    def __init__(self, source, msg):
        self.args = (source, msg)
        self.source = source
        self.msg = msg
        self.fullmsg = "Unable to parse {}: {}".format(source,msg)

class ParseWarning(Exception):
    def __init__(self, source, msg):
        self.args = (source, msg)
        self.source = source
        self.msg =  msg
        self.fullmsg = "Warning: parsefail in {}: {}".format(source,msg)
        