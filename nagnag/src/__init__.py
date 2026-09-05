''' App design components
* `ui` (`shiny.Tag`)
* `server` (`function`)
'''

__all__ = ["ui", "server"]

from .app_ui import app_ui as ui
from .app_server import server
