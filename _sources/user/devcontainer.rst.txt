.. _devcontainer:

*********************
Dev Container
*********************

The edrixs repository includes a `Dev Container
<https://containers.readthedocs.io/>`_ configuration, which gives you a
fully working edrixs development environment without any manual installation.

Supported environments
----------------------

* **VS Code** with the `Dev Containers
  <https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers>`_
  extension
* **GitHub Codespaces** — runs the container in the cloud directly from the
  GitHub repository page, with no local setup at all
* **JetBrains IDEs** (PyCharm, IntelliJ, etc.) via their Dev Containers support
* **The devcontainer CLI** (``@devcontainers/cli``) for terminal-only use

Using VS Code
-------------

1. Install `Docker <https://www.docker.com/>`_ and start the Docker daemon.
2. Install the `Dev Containers
   <https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers>`_
   extension in VS Code.
3. Clone the edrixs repository and open it in VS Code.
4. When prompted *"Reopen in Container"*, click it — or open the Command
   Palette (``Ctrl+Shift+P``) and run **Dev Containers: Reopen in Container**.

VS Code will build the image on the first launch (this takes a few minutes)
and then reopen with edrixs fully installed and ready to use.

Using GitHub Codespaces
-----------------------

1. Go to the edrixs repository on GitHub.
2. Click **Code → Codespaces → Create codespace on master**.

The codespace opens a browser-based VS Code session with edrixs already
installed — no local software required.
