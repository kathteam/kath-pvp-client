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

if os.path.exists("dist/main.app"):
    shutil.rmtree("dist/main.app")

ENTRY_POINT = ["backend/main.py"]

DATA_FILES = tree("gui")
OPTIONS = {
    "argv_emulation": True,
    "strip": False,
    "iconfile": "backend/assets/logo.icns",
    "packages": ["WebKit", "Foundation", "webview"],
    "plist": {
        "CFBundleName": "MainApp",
        "CFBundleShortVersionString": "0.1.0",
        "CFBundleVersion": "0.1.0",
        "CFBundleIdentifier": "com.example.mainapp",
    },
    # "codesign_identity": "Developer ID Application: Your Name (Team ID)",
    # "entitlements_file": "path/to/your/entitlements.plist",
    "resources": DATA_FILES,
}

setup(
    app=ENTRY_POINT,
    options={"py2app": OPTIONS},
    setup_requires=["py2app"],
)
