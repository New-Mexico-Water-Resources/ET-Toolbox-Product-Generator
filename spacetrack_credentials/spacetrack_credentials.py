from typing import Dict
from os.path import join, abspath, dirname
from credentials import get_credentials

FILENAME = join(abspath(dirname(__file__)), "spacetrack_credentials.txt")

def get_spacetrack_credentials(filename: str = FILENAME) -> Dict[str, str]:
    return get_credentials(
        filename=filename,
        displayed=["username"],
        hidden=["password"],
        prompt="credentials for Spacetrack https://www.space-track.org/auth/createAccount"
    )
