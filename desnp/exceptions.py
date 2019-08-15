# -*- coding: utf-8 -*-

#
# Collection of module errors
#


class DeSNPError(Exception):
    """
    Simple base exception, root of all DeSNP exceptions
    """
    def __init__(self, msg=None):
        self.msg = msg
        super(DeSNPError, self).__init__(self.msg)


class KeyboardInterruptError(Exception):
    """
    Keyboard Interrupt errors
    """
    pass
