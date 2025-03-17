# kath-pvp-client

The desktop application frontend is developed using Vite, TypeScript and React for user interface. On the backend, Python is used in combination with the PyWebView library, which enables embedding a web-based UI inside a native desktop window. PyWebView acts as a bridge between the frontend and backend, facilitating seamless communication for executing logic, handling user actions and processing data.

## Prerequisites
- [Python 3](https://www.python.org/downloads/)
- [Node.js](https://nodejs.org/en/download)

## Initialization

### Windows

```shell
npm run init
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

### Linux

```shell
npm run init
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Mac OS

```shell
npm run init
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Development

### Only the frontend

```shell
npm run dev-frontend
```

### Build the frontend and run the backend in dev mode

```shell
npm run start
```

## Build

```shell
npm run build
```

This will create an executable file in `./dist` directory.