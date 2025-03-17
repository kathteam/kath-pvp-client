# TODO Fix the following block to enable macOS builds

import os
import shutil
from distutils.core import setup


def tree(src):
    return [
        (root, map(lambda f: os.path.join(root, f), files))
        for (root, dirs, files) in os.walk(os.path.normpath(src))
    ]


if os.path.exists("build"):
    shutil.rmtree("build")

if os.path.exists("dist/Kath.app"):
    shutil.rmtree("dist/Kath.app")

ENTRY_POINT = ["backend/main.py"]

DATA_FILES = tree("gui")
OPTIONS = {
    "argv_emulation": True,
    "strip": False,
    "iconfile": "backend/assets/logo.icns",
    "packages": ["WebKit", "Foundation", "webview"],
    "plist": {"NSRequiresAquaSystemAppearance": False},
    "resources": DATA_FILES,
}

setup(
    name="Kath",
    app=ENTRY_POINT,
    options={"py2app": OPTIONS},
    setup_requires=["py2app"],
)
