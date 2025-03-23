import uvicorn
import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    uvicorn.run("app.main:app", host="0.0.0.0", port=8001, reload=True)