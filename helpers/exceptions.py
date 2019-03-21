# Custom exceptions used while casting a variant

class ParseException(Exception):
    '''
    Critical failure parsing input files
    '''
    def __init__(self, source, msg):
        self.args = (source, msg)
        self.source = source
        self.msg = msg
        self.fullmsg = "Unable to parse {}: {}".format(source,msg)

class ParseWarning(Exception):
    '''
    Noncritical failure parsing input files
    '''
    def __init__(self, source, msg):
        self.args = (source, msg)
        self.source = source
        self.msg =  msg
        self.fullmsg = "Warning: parsefail in {}: {}".format(source,msg)

class AlignException(Exception):
    '''
    Critical failure during alignment, full worker fail
    '''
    def __init__(self, source, msg):
        self.args = (source, msg)
        self.source = source
        self.msg = msg
        self.fullmsg = "{} failed align: {}".format(source, msg)

class WorkerException(Exception):
    '''
    Isolated worker failure to PDB or model lookup
    '''
    def __init__(self, source, msg):
        self.args = (source, msg)
        self.source = source
        self.msg = msg
        self.fullmsg = "worker: {} failed: {}".format(source, msg)

class DescriptorException(Exception):
    '''
    Failure during descriptor attachment
    '''
    def __init__(self, source, msg):
        self.args = (source,msg)
        self.source = source
        self.msg = msg
        self.fullmsg = "descriptor gen failure: {} {}".format(source,msg)

class NoModelException(Exception):
    '''
    No custom models match (avoids writing as error since expected to be common
    '''
    def __init__(self):
        self.msg = 'no custom model for transcript'
        self.fullmsg = self.msg

class CustomModelException(Exception):
    '''
    Failed to load a custom model, report this because it is unusual
    '''
    def __init__(self, source, msg):
        self.args = (source, msg)
        self.source = source
        self.msg = msg
        self.fullmsg = "Custom model loading gen failure: {} {}".format(source,msg)
                
