from functools import wraps
import platform


class MetaStandbyLock(type):
    """  pong """

    SYSTEM = platform.system()

    def __new__(cls, name: str, bases: tuple, attrs: dict) -> type:
        if not ('inhibit' in attrs and 'release' in attrs):
            raise TypeError(
                "Missing implementations for classmethods 'inhibit(cls)' and 'release(cls)'."
            )
        else:
            if name == 'StandbyLock':
                cls._superclass = super().__new__(cls, name, bases, attrs)
                return cls._superclass
            if cls.SYSTEM.upper() in name.upper():
                if not hasattr(cls, '_superclass'):
                    raise ValueError(
                        "Class 'StandbyLock' must be implemented.")
                cls._superclass._subclass = super().__new__(
                    cls, name, bases, attrs)
                return cls._superclass._subclass
            else:
                return super().__new__(cls, name, bases, attrs)


class StandbyLock(metaclass=MetaStandbyLock):
    """
    """

    _subclass = None

    @classmethod
    def inhibit(cls):
        if cls._subclass is None:
            raise OSError(
                f"There is no 'StandbyLock' implementation for OS '{platform.system()}'."
            )
        else:
            return cls._subclass.inhibit()

    @classmethod
    def release(cls):
        if cls._subclass is None:
            raise OSError(
                f"There is no 'StandbyLock' implementation for OS '{platform.system()}'."
            )
        else:
            return cls._subclass.release()

    def __enter__(self, *args, **kwargs):
        self.inhibit()
        return self

    def __exit__(self, *args, **kwargs):
        self.release()


class WindowsStandbyLock(StandbyLock):
    """
    """

    ES_CONTINUOUS = 0x80000000
    ES_SYSTEM_REQUIRED = 0x00000001

    INHIBIT = ES_CONTINUOUS | ES_SYSTEM_REQUIRED
    RELEASE = ES_CONTINUOUS

    @classmethod
    def inhibit(cls):
        import ctypes
        ctypes.windll.kernel32.SetThreadExecutionState(cls.INHIBIT)

    @classmethod
    def release(cls):
        import ctypes
        ctypes.windll.kernel32.SetThreadExecutionState(cls.RELEASE)


class LinuxStandbyLock(metaclass=MetaStandbyLock):
    """
    """

    COMMAND = 'systemctl'
    ARGS = [
        'sleep.target', 'suspend.target', 'hibernate.target',
        'hybrid-sleep.target'
    ]

    @classmethod
    def inhibit(cls):
        import subprocess
        subprocess.run([cls.COMMAND, 'mask', *cls.ARGS])

    @classmethod
    def release(cls):
        import subprocess
        subprocess.run([cls.COMMAND, 'unmask', *cls.ARGS])


class DarwinStandbyLock(metaclass=MetaStandbyLock):
    """
    """

    COMMAND = 'caffeinate'
    BREAK = b'\003'

    _process = None

    @classmethod
    def inhibit(cls):
        from subprocess import Popen, PIPE
        cls._process = Popen([cls.COMMAND], stdin=PIPE, stdout=PIPE)

    @classmethod
    def release(cls):
        cls._process.stdin.write(cls.BREAK)
        cls._process.stdin.flush()
        cls._process.stdin.close()
        cls._process.wait()


def standby_lock(callback):
    """ standby_lock(callable) -> callable
        This decorator guarantees that the system will not enter standby mode while 'callable' is running.
    """
    @wraps(callback)
    def new_callback(*args, **kwargs):
        with StandbyLock():
            return callback(*args, **kwargs)

    return new_callback


@standby_lock
def test_function():

    for _ in range(2000):
        import time
        time.sleep(5)
