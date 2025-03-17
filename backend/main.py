import os
import webview


class Api:
    def fullscreen(self) -> None:
        webview.windows[0].toggle_fullscreen()


def get_entrypoint() -> str:
    def exists(path):
        return os.path.exists(os.path.join(os.path.dirname(__file__), path))

    if exists("../gui/index.html"):
        return "../gui/index.html"

    if exists("../Resources/gui/index.html"):
        return "../Resources/gui/index.html"

    if exists("./gui/index.html"):
        return "./gui/index.html"

    raise Exception("No index.html found")


if __name__ == "__main__":
    window = webview.create_window(title="", url=get_entrypoint(), js_api=Api())
    webview.start()
