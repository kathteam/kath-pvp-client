# kath-pvp-client

The desktop application frontend is developed using Vite, TypeScript and React for user interface. On the backend, Python is used in combination with the PyWebView library, which enables embedding a web-based UI inside a native desktop window. PyWebView acts as a bridge between the frontend and backend, facilitating seamless communication for executing logic, handling user actions and processing data.

**Note:** _All commands must be run from the project root directory._

## Prerequisites
- [Python 3](https://www.python.org/downloads/) (recommended version 3.12)
- [Node.js](https://nodejs.org/en/download/) (recommended version 22.x)

## Initialization

### Windows

```shell
npm run init
python -m venv .venv 
.venv\Scripts\activate
pip install -r requirements.txt
pip install -r requirements_dev.txt
```

### Linux

```shell
npm run init
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -r requirements_dev.txt
```

### Mac OS

Currently not supported (good luck).

## Development

### Frontend

```shell
npm run dev
```
This starts a development server for the frontend at `http://localhost:5173/`.

### Frontend & Backend

```shell
npm run start
```

This builds the frontend and runs the backend in development mode.

### Linting

```shell
npm run lint
```

For frontend directory this will check and automatically fix the Typescript code using ESLint. In the backend directory it will check the Python code using pylint.

## Build

```shell
npm run build
```

This generates a production-ready build, placing the output in the `./dist` directory.