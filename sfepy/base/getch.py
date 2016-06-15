"""
getch()-like unbuffered character reading from stdin on both Windows and Unix

_Getch classes inspired by Danny Yoo, iskeydown() based on code by Zachary
Pincus.
"""
from __future__ import absolute_import
from __future__ import print_function
import os, sys

class _Getch:
    """Gets a single character from standard input. Does not echo to the
    screen."""
    def __init__(self):
        if os.name == 'nt':
            self.impl = _GetchWindows()
        elif os.name == 'posix':
            self.impl = _GetchUnix()
        else:
            self.impl = _GetchDefault()

    def __call__(self): return self.impl()

class _GetchWindows:
    def __init__(self):
        import msvcrt
        self.msvcrt = msvcrt

    def __call__(self):
        msvcrt = self.msvcrt
        
        c = msvcrt.getch()
        if c == '\x00' or c == '\xE0':    #functions keys
            msvcrt.getch()
        return c

    def iskeydown(self):
        msvcrt = self.msvcrt
        return msvcrt.kbhit()

class _GetchUnix:
    def __init__(self):
        import tty, sys, select
        self.mods = (tty, sys, select)

    def __call__(self):
        tty, sys, select = self.mods

        fd = sys.stdin.fileno()
        old = tty.tcgetattr(fd)
        tty.setcbreak(fd, tty.TCSANOW)
        try:
            return sys.stdin.read(1)
        finally:
            tty.tcsetattr(fd, tty.TCSAFLUSH, old)

    def iskeydown(self):
        tty, sys, select = self.mods

        fd = sys.stdin.fileno()
        old = tty.tcgetattr(fd)
        tty.setcbreak(fd, tty.TCSANOW)
        try:
            if select.select([sys.stdin], [], [], 0) == ([sys.stdin], [], []):
                return sys.stdin.read(1)
            else:
                return False
        finally:
            tty.tcsetattr(fd, tty.TCSAFLUSH, old)

class _GetchDefault:
    def __call__(self):
        return raw_input()[0]

    def iskeydown(self):
        raise NotImplementedError

getch = _Getch()

if __name__ == '__main__':
    from .base import pause, spause
    
    pause('press a key anytime the script stops!')
    pause()
    spause('last time...')
    print('done.')
