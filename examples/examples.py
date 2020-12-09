

# accessing a class' name:
import math
print(math.pi.__class__.__name__)
print(math.pi.__class__)
print(type(math.pi).__name__)

# another way to import something:
mod = __import__('math')
print(mod.pi)

import src.feedingorcas
from pprint import pprint


print(src.feedingorcas.__doc__)


# things which have to work properly:
pprint(dir(src.feedingorcas))            # presented namespace

pprint(src.feedingorcas.__dict__)   # dict of namespace ()

print(src.feedingorcas.__doc__)   # grabs a single module/package/function/class/variable



# pprint(globals())  # global namespace
# print('location: examples.py, package:', __package__)
# print('location: examples.py, file:', __file__)
# print('location: examples.py, name:', __name__)

# # built-in namespace:
# dir(__builtins__)

# import functools
#
# # default cache maxsize = 1
# @functools.lru_cache(maxsize=3)
# def square(x: float) -> float:
#     print(f'running: {x}')
#     return x * x
