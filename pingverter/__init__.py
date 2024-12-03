
from .humminbird_class import hum
from .lowrance_class import low
from .converter import hum2pingmapper, low2pingmapper, low2hum
from . import _version
__version__ = _version.get_versions()['version']
