# -*- coding: utf-8 -*-

#
# Collection of desnp utility functions and classes
#

import bz2
import gzip
import os
import urllib

from .exceptions import DeSNPError
from . import compat


def check_file(filename, mode='r'):
    if mode == 'r':

        if filename and os.path.exists(filename):
            return os.path.abspath(filename)

        raise DeSNPError(
            "The following file does not exist: {0}".format(filename))

    elif mode == 'w':
        file_dir = '.'

        if filename:
            file_name = os.path.abspath(filename)
            file_dir = os.path.dirname(file_name)

            if not os.access(file_dir, os.W_OK | os.X_OK):
                raise DeSNPError("Cannot generate file: {0}".format(filename))

            return file_name

    raise DeSNPError("Unspecified mode to open file, '{}'".format(mode))


def s(value):
    if compat.is_py2:
        if isinstance(value, unicode):
            return value.decode('ascii', 'ignore')

        if isinstance(value, str):
            return value

    else:
        if isinstance(value, bytes):
            return n(value)

        if isinstance(value, str):
            return value

    return value


def open_resource(resource, mode='rb'):
    """
    Open different types of files and return the handle.

    :param resource: a file located locally or on the internet.
                     Gzip'd and zip'd files are handled.
    :param mode: standard file open modes
    :return: the resource (file) handle
    """
    if not resource:
        return None

    if compat.is_py2:
        if not isinstance(resource, basestring):
            return resource
    else:
        if not isinstance(resource, str):
            return resource

    resource = s(resource)

    if resource.endswith(('.gz', '.Z', '.z')):
        return gzip.open(resource, mode)
    elif resource.endswith(('.bz', '.bz2', '.bzip2')):
        return bz2.BZ2File(resource, mode)
    elif resource.startswith(('http://', 'https://', 'ftp://')):
        return urllib.urlopen(resource)
    else:
        return open(resource, mode)
